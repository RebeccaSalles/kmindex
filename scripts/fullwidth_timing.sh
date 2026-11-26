#!/usr/bin/env bash
# Full-width timing on your own machine or node: generates a synthetic one-group index as wide as full Logan
# (gen_index.py), then times v0.7.1 and the PR build, every mode, cold, on the same local disk.
# Requirements and disk space: README §3. Runs as `bash` or as a SLURM job, e.g.
#   sbatch -p <partition> --nodelist=<node> --exclusive --cpus-per-task=16 --mem=32G --time=03:00:00 \
#          --export=ALL,KM_REF=...,KM_PR=...,OUT=...,SRC=... scripts/fullwidth_timing.sh
# Settings (environment): those of common.sh, plus
#   SRC      your kmindex v0.7.1 checkout (its tests/data gives the 49-byte file header and the repartition format)
#   SAMPLES  samples per row (default 26857911, full Logan)
#   PARTS    partition files (default 256)        ROWS  rows per file, a multiple of 64 (default 128)
#   PROCS    processes for the generator (default 8)
#   KEEP=1   keep the generated index after the run (default: deleted)
set -euo pipefail
source "$(dirname "$0")/common.sh"
SRC=${SRC:?set SRC to your kmindex v0.7.1 checkout}
SAMPLES=${SAMPLES:-26857911}; PARTS=${PARTS:-256}; ROWS=${ROWS:-128}; PROCS=${PROCS:-8}
IDX=$LOCAL/kmi_ev_fullwidth_$$
BYTES=$(( PARTS * (49 + ROWS * ((SAMPLES + 7) / 8)) ))
FREE=$(( $(df --output=avail -B1 "$LOCAL" | tail -1) ))
echo "index: $SAMPLES samples, $PARTS files x $ROWS rows = $(( BYTES / 1000000000 )) GB in $LOCAL (free: $(( FREE / 1000000000 )) GB)"
(( FREE > BYTES + 2000000000 )) || { echo "not enough space in $LOCAL: lower ROWS or PARTS (README §3)" >&2; exit 6; }
[[ ${KEEP:-0} == 1 ]] || trap 'rm -rf "$IDX"' EXIT
{ echo "host=$(hostname) $(date -Iseconds) samples=$SAMPLES parts=$PARTS rows=$ROWS fstype=$(df --output=fstype "$LOCAL" | tail -1)"
  lsblk -d -o NAME,ROTA,SIZE,MODEL 2>/dev/null | grep -v loop || true
  echo "ref md5=$(md5sum < "$KM_REF" | cut -c1-12)  pr md5=$(md5sum < "$KM_PR" | cut -c1-12)"; } | tee "$OUT/host.txt"
python3 "$(dirname "$0")/gen_index.py" "$IDX" --samples "$SAMPLES" --partitions "$PARTS" --rows "$ROWS" --procs "$PROCS" \
  --density 0.3125 --hit-frac 0.001 \
  --header-from "$SRC/tests/data/indexes/pa_index/matrices/matrix_0.cmbf" \
  --repart-from "$SRC/tests/data/indexes/pa_index/repartition_gatb/repartition.minimRepart" | tail -2 | tee -a "$OUT/host.txt"
CSV=$OUT/fullwidth_timings.csv
echo "mode,query_length,query,state,resident_bytes,elapsed_s,user_s,sys_s,maxrss_kb,out_md5" > "$CSV"
read -ra M <<< "$MODES"; ROT=0
for L in $LENGTHS; do
  for Q in "$QUERIES/$L"/*.fa; do
    ROT=$((ROT + 1))
    for i in "${!M[@]}"; do                      # the order of the modes rotates from query to query
      mode=${M[$(( (i + ROT) % ${#M[@]} ))]}; setmode "$mode"
      evict "$IDX"; RES=$(resident "$IDX")
      run1 "$IDX" "$Q" "$mode,$L,$(basename "$Q" .fa),$(state_of "$RES"),$RES" "$CSV"
    done
  done
done
echo "done: $CSV"
python3 - "$CSV" <<'PY'
import csv, sys, statistics as st, collections as c
d = c.defaultdict(list); md5 = c.defaultdict(set)
for r in csv.DictReader(open(sys.argv[1])):
    d[(r['mode'], int(r['query_length']))].append(float(r['elapsed_s'])); md5[(r['query'])].add(r['out_md5'])
for k in sorted(d, key=lambda k: (k[1], k[0])): print(f'{k[1]:>5} bp  {k[0]:7s} median {st.median(d[k]):7.2f} s  ({len(d[k])} queries)')
print('outputs identical across modes for every query:', all(len(v) == 1 for v in md5.values()))
PY
