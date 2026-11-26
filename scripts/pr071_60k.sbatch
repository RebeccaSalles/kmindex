#!/usr/bin/env bash
#SBATCH --job-name=kmi_pr071_60k
#SBATCH --partition=gpu
#SBATCH --nodelist=genogpu2
#SBATCH --exclusive
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=14:00:00
# [2026-11-25] 60k / 90k campaign with the PR code (kmindex v0.7.1 + pr1..pr4, bin_x; switches off = untouched code
# path, identity-tested): modes none (today), in16 (whole-row reads), in16fcT8 (+ fast count, no early stop), st32
# (streamed rows). Same designs, protocol and anchors as opt_st. For the 19 Nov deck's slide 21 and the model refit.
# Predictions: pr071/PREDICTIONS_60k.txt.
# Previous header:
# [2026-11-16] Same session as 17705724 (opt_all), two modes, no other change:
#   in16fcT8 : record build (prefetch 16 + fast count, 8 counting threads, ONE reader) -> fit the
#              model-with-fixes constants on exactly the configuration that model describes;
#   st32     : stream build (KMINDEX_STREAM=32: rows read in query order into a ring of 32 + z
#              rows, counters updated as they arrive; same prefetch and counting) -> memory and
#              time on real narrow-row indexes. Outputs compared with opt_all's "off" afterwards.
# Original header:
# Prefetch campaign (uncommitted patched kmindex, prefetch_campaign/bin): the SSD session
# (genogpu2, Dell BOSS) with every query run once per prefetch mode, back to back on the same
# local copy, cold each time, rotating order. Env: PF_NAME (results), PREFETCH_MODES.
# Whole node: 60k anchors,
# the 16 designs on both samples (60k and 90k alternated per design, the
# order alternated between designs), 60k anchors.
# Same verified-cold local protocol (requery_cold_local_pf.sh) throughout.
set -euo pipefail
P=/projects/logan_compression/rsalles/kmhelpers_logan_paper_pfactor/scripts
N=/projects/logan_compression/rsalles/kmhelpers_logan_paper_90k/scripts
W=/projects/logan_compression/rsalles/span_tuning_paper
R=/projects/logan_compression/rsalles/prefetch_campaign
export KMINDEX_BIN_DIR_OVERRIDE=$R/pr071/bin_x STREAM_BIN_DIR=$R/pr071/bin_x PREFETCH_MODES="${PREFETCH_MODES:-none in16 in16fcT8 st32}"
PF_NAME=${PF_NAME:-opt_pr071}
O=$W/$PF_NAME
mkdir -p $O
echo "$(date -Iseconds) session on $(hostname)"
ANCH="pf_b1p1_g1_p4_t8:4 pf_b1p1_g1_p32_t8:32"
run() {  # results_name phase build_root config_entry
  RESULTS_NAME_OVERRIDE="$1" PHASE_OVERRIDE="$2" BUILD_ROOT_OVERRIDE="$3" CONFIGS_OVERRIDE="$4" bash $R/requery_cold_local_pf2.sh
}
run ${PF_NAME}/anchors_before ${PF_NAME}_anchors_before $W/pfactor_builds "$ANCH"
i=0
while read -r id base g e60 m60 e90 m90; do
  E60="$id:4:$base:$g:$e60:$m60"; E90="$id:4:$base:$g:$e90:$m90"
  if (( i % 2 == 0 )); then
    run ${PF_NAME}/s60k ${PF_NAME}_60k $W/training_rerun_builds "$E60"
    run ${PF_NAME}/s90k ${PF_NAME}_90k $W/training90k_builds "$E90"
  else
    run ${PF_NAME}/s90k ${PF_NAME}_90k $W/training90k_builds "$E90"
    run ${PF_NAME}/s60k ${PF_NAME}_60k $W/training_rerun_builds "$E60"
  fi
  i=$((i + 1))
done <<'CFG'
b1p1_g02 1.1 2 17700000000 8870000000 26550000000 13305000000
b1p05_g249 1.05 249 6550000000 336000000 9825000000 504000000
b1p1_g20 1.1 20 7110000000 628000000 10665000000 942000000
b1p5_g32 1.5 32 7810000000 2470000000 11715000000 3705000000
b1p2_g70 1.2 70 6970000000 1260000000 10455000000 1890000000
b1p1_g131 1.1 131 6670000000 628000000 10005000000 942000000
b1p3_g04 1.3 4 10200000000 3750000000 15300000000 5625000000
b1p75_g07 1.75 7 8680000000 3520000000 13020000000 5280000000
b1p05_g03 1.05 3 12500000000 4400000000 18750000000 6600000000
b1p4_g02 1.4 2 21400000000 12900000000 32100000000 19350000000
b1p5_g02 1.5 2 17500000000 12500000000 26250000000 18750000000
b1p5_g03 1.5 3 13100000000 5000000000 19650000000 7500000000
b1p75_g02 1.75 2 17300000000 13500000000 25950000000 20250000000
b1p75_g03 1.75 3 13000000000 5750000000 19500000000 8625000000
b2p0_g03 2.0 3 15200000000 9090000000 22800000000 13635000000
b2p0_g04 2.0 4 10800000000 4380000000 16200000000 6570000000
CFG
run ${PF_NAME}/anchors_after ${PF_NAME}_anchors_after $W/pfactor_builds "$ANCH"
echo "$(date -Iseconds) session done"
