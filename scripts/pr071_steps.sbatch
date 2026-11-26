#!/usr/bin/env bash
#SBATCH --job-name=kmi_pr071steps
#SBATCH --partition=gpu
#SBATCH --nodelist=genogpu2
#SBATCH --exclusive
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=08:00:00
# [2026-11-25] STEPS run for the 19 Nov deck: every optimization step on the PR code (v0.7.1 + pr1..pr4) plus,
# for measurement only, names read only for the hits (lazy_measurement_only.patch; not in the PRs). bin_ref from
# job 17744756 is reused. Predictions: PREDICTIONS_steps.txt.
# Previous header:
# [2026-11-24] The PRs as proposed for Friday's review (kmindex v0.7.1 + pr1..pr4: early stop removed from the fast
# count, macOS guard for posix_fadvise), built next to v0.7.1 untouched, then timed on the full-width index with the
# protocol of 17707000 / 17713887. Predictions: PREDICTIONS_pr071.txt.
set -uo pipefail
R=/projects/logan_compression/rsalles/prefetch_campaign
W=/projects/logan_compression/rsalles/span_tuning_paper
D=$R/pr071
OUT=$D/results_${SLURM_JOB_ID}; mkdir -p $OUT
# ---- builds (system g++, cmake-only conda env, as every earlier build)
export PATH="/projects/logan_compression/rsalles/zst_native_patch_build/buildenv/bin:/usr/bin:/bin"; export CC=gcc CXX=g++
cmake --version | head -1 | tee $OUT/build.txt; g++ --version | head -1 | tee -a $OUT/build.txt
rm -rf $D/src_x && cp -a $D/kmindex_v071 $D/src_x && cd $D/src_x && for p in $D/pr1_whole_row_reads.patch $D/pr2_fast_count.patch $D/pr3_streamed_rows.patch $D/pr4_one_copy_of_names.patch $D/lazy_measurement_only.patch; do git apply $p || { echo "patch $p failed"; exit 2; }; done
git diff --stat | tail -1 | tee -a $OUT/build.txt
rm -rf $D/build_x && mkdir -p $D/build_x && cd $D/build_x
cmake $D/src_x -DCMAKE_BUILD_TYPE=Release -DWITH_SERVER=OFF -DWITH_COMPRESSION=OFF -DPORTABLE_BUILD=ON -DCMAKE_IGNORE_PREFIX_PATH=/local/miniconda3 > cmake.log 2>&1 || { tail -30 cmake.log; exit 3; }
make -j16 > make.log 2>&1 || { grep -n "error" make.log | head -30; exit 4; }
B=$(find $D/build_x -name kmindex -type f -perm -u+x | head -1); mkdir -p $D/bin_x; cp $B $D/bin_x/kmindex
echo "x md5=$(md5sum < $D/bin_x/kmindex | cut -c1-12)  ref md5=$(md5sum < $D/bin_ref/kmindex | cut -c1-12)" | tee -a $OUT/build.txt
# ---- timing session (as fullwidth_v5_lean.sbatch)
source /projects/logan_compression/rsalles/kmhelpers_logan_paper_pfactor/scripts/activate_env.sh   # numpy for gen_index
REAL_G=$W/pfactor_builds/pf_b1p1_g1_p32_t8/kmindex_data/eval/pf_b1p1_g1_p32_t8_g0
MODES=(ref pr1 pr12 pr123 pr1234 pr12345)
LOCAL=/tmp/pr071s_${SLURM_JOB_ID}; mkdir -p $LOCAL; trap 'rm -rf $LOCAL' EXIT
{ echo "host=$(hostname) job=$SLURM_JOB_ID $(date -Iseconds) fstype=$(df --output=fstype $LOCAL | tail -1)"; lsblk -d -o NAME,ROTA,SIZE,MODEL | grep -v loop; } | tee $OUT/host.txt
mkdir -p $LOCAL/q/31; for L in 250 500 1000; do mkdir -p $LOCAL/q/$L; cp $W/queries/$L/q00[0-4]_$L.fa $LOCAL/q/$L/; done
python3 - $LOCAL/q/250/q000_250.fa $LOCAL/q/31/q000_31.fa <<'PY'
import sys
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
open(sys.argv[2], "w").write(">q000_31\n" + seq[:31] + "\n")
PY
evict() { sync; python3 - "$LOCAL/idx" <<'PY'
import os, sys
for root, _, files in os.walk(sys.argv[1]):
    for n in files:
        p = os.path.join(root, n)
        if os.path.islink(p): continue
        fd = os.open(p, os.O_RDONLY)
        try: os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        finally: os.close(fd)
PY
}
resident() { find "$LOCAL/idx" -type f -print0 | xargs -0 -r fincore -b -n -o RES 2>/dev/null | awk '{s+=$1} END{print s+0}'; }
setmode() {
  unset KMINDEX_PREFETCH KMINDEX_FASTCOUNT KMINDEX_COUNT_THREADS KMINDEX_STREAM KMINDEX_LEAN KMINDEX_LAZY_NAMES
  BIN=$D/bin_x/kmindex
  case "$1" in
    ref) BIN=$D/bin_ref/kmindex ;;
    pr1) export KMINDEX_PREFETCH=16 ;;
    pr12) export KMINDEX_PREFETCH=16 KMINDEX_FASTCOUNT=1 KMINDEX_COUNT_THREADS=8 ;;
    pr123) export KMINDEX_PREFETCH=16 KMINDEX_FASTCOUNT=1 KMINDEX_COUNT_THREADS=8 KMINDEX_STREAM=32 ;;
    pr1234) export KMINDEX_PREFETCH=16 KMINDEX_FASTCOUNT=1 KMINDEX_COUNT_THREADS=8 KMINDEX_STREAM=32 KMINDEX_LEAN=1 ;;
    pr12345) export KMINDEX_PREFETCH=16 KMINDEX_FASTCOUNT=1 KMINDEX_COUNT_THREADS=8 KMINDEX_STREAM=32 KMINDEX_LEAN=1 KMINDEX_LAZY_NAMES=1 ;;
  esac
}
CSV=$OUT/timings.csv
echo "mode,query_length,query,rep,cache_state,resident_bytes,elapsed_s,user_s,sys_s,maxrss_kb,out_md5,out_bytes" > $CSV
run1() {  # mode L Q rep state threshold resident
  local M=$1 L=$2 Q=$3 REP=$4 ST=$5 THR=$6 RES=$7 QB=$(basename $3 .fa) O=$LOCAL/out
  rm -rf $O
  /usr/bin/time -f '%e,%U,%S,%M' -o $LOCAL/t $BIN query -i $LOCAL/idx -q $Q -o $O -t 4 -z 6 -r $THR > $OUT/${QB}_${M}_${REP}_${THR}.stdout 2> $OUT/${QB}_${M}_${REP}_${THR}.stderr
  echo "$M,$L,$QB,$REP,$ST,$RES,$(cat $LOCAL/t),$(find $O -name '*.json' -exec cat {} + | md5sum | cut -c1-12),$(find $O -name '*.json' -exec cat {} + | wc -c)" | tee -a $CSV
}
echo "$(date -Iseconds) generating index"
python3 $R/fullwidth/gen_index.py $LOCAL/idx --rows 128 --procs 12 --density 0.3125 --hit-frac 0.001 \
  --header-from $REAL_G/matrices/matrix_0.cmbf --repart-from $REAL_G/repartition_gatb/repartition.minimRepart 2>&1 | tail -3 | tee $OUT/gen.log
sync
for M in "${MODES[@]}"; do setmode $M; for REP in 1 2 3; do run1 $M 31 $LOCAL/q/31/q000_31.fa $REP floor 0.99 $(resident); done; done
for M in "${MODES[@]}"; do setmode $M; run1 $M 250 $LOCAL/q/250/q000_250.fa 1 fullout 0 $(resident); done
ROT=0
for L in 1000 250 500; do
  for Q in $LOCAL/q/$L/*.fa; do
    ROT=$((ROT + 1))
    for i in 0 1 2 3 4 5; do
      M=${MODES[$(( (i + ROT) % 6 ))]}; setmode $M
      evict; RES=$(resident); (( RES <= 16777216 )) && ST=cold || ST=cold_unverified
      run1 $M $L $Q 1 $ST 0.99 $RES
      run1 $M $L $Q 2 warm 0.99 $(resident)
    done
  done
done
echo "$(date -Iseconds) pr071 steps done"
