#!/usr/bin/env bash
# Real-index timing: copies already-built kmindex indexes (registries) to node-local disk, rewrites the absolute paths
# kmindex stores in them, then times v0.7.1 and the PR build, every mode, cold, on that copy.
# Usage:  scripts/realindex_timing.sh <index dir> [<index dir> ...]
# Our 16 designs on Genouest (readable by the group prj_logan_compression):
#   /projects/logan_compression/rsalles/span_tuning_paper/training_rerun_builds/<design>   (60k samples)
#   /projects/logan_compression/rsalles/span_tuning_paper/training90k_builds/<design>      (90k samples)
# As a SLURM job, e.g.
#   sbatch -p <partition> --nodelist=<node> --exclusive --cpus-per-task=16 --mem=32G --time=04:00:00 \
#          --export=ALL,KM_REF=...,KM_PR=...,OUT=... scripts/realindex_timing.sh <index dir> ...
# Settings: those of common.sh (MODES default: ref pr1 pr12 pr123 pr1234). Requirements: README §3.
set -euo pipefail
source "$(dirname "$0")/common.sh"
(( $# > 0 )) || { echo "usage: $0 <index dir> [...]" >&2; exit 1; }
CSV=$OUT/realindex_timings.csv
[[ -f $CSV ]] || echo "design,mode,query_length,query,state,resident_bytes,elapsed_s,user_s,sys_s,maxrss_kb,out_md5" > "$CSV"
echo "host=$(hostname) $(date -Iseconds) ref md5=$(md5sum < "$KM_REF" | cut -c1-12) pr md5=$(md5sum < "$KM_PR" | cut -c1-12)" >> "$OUT/host.txt"
read -ra M <<< "$MODES"; ROT=0
for SRCIDX in "$@"; do
  SRCIDX=$(cd "$SRCIDX" && pwd); NAME=$(basename "$SRCIDX")
  [[ -f $SRCIDX/index.json ]] || { echo "$SRCIDX has no index.json: not a kmindex registry" >&2; exit 7; }
  NEED=$(du -sb "$SRCIDX" | awk '{print $1}'); FREE=$(( $(df --output=avail -B1 "$LOCAL" | tail -1) ))
  (( FREE > NEED + 2000000000 )) || { echo "not enough space in $LOCAL for $NAME ($(( NEED / 1000000000 )) GB)" >&2; exit 6; }
  IDX=$LOCAL/kmi_ev_${NAME}_$$; rm -rf "$IDX"; trap 'rm -rf "$IDX"' EXIT
  cp -a "$SRCIDX" "$IDX"
  # kmindex registries store absolute paths (index.json, symlinks): point them at the local copy
  while IFS= read -r -d '' l; do t=$(readlink "$l"); if [[ $t == "$SRCIDX"* ]]; then ln -sfn "$IDX${t#"$SRCIDX"}" "$l"; fi; done < <(find "$IDX" -type l -print0)
  { grep -rlIZ --fixed-strings "$SRCIDX" "$IDX" 2>/dev/null || true; } | xargs -0 -r sed -i "s#$SRCIDX#$IDX#g"
  if grep -rqI --fixed-strings "$SRCIDX" "$IDX" || [[ -n $(find "$IDX" -type l -lname "$SRCIDX*" -print -quit) ]]; then
    echo "the local copy of $NAME still points to $SRCIDX" >&2; exit 8; fi
  for L in $LENGTHS; do
    for Q in "$QUERIES/$L"/*.fa; do
      ROT=$((ROT + 1))
      for i in "${!M[@]}"; do
        mode=${M[$(( (i + ROT) % ${#M[@]} ))]}; setmode "$mode"
        evict "$IDX"; RES=$(resident "$IDX")
        run1 "$IDX" "$Q" "$NAME,$mode,$L,$(basename "$Q" .fa),$(state_of "$RES"),$RES" "$CSV"
      done
    done
  done
  rm -rf "$IDX"; echo "$(date -Iseconds) $NAME done"
done
python3 - "$CSV" <<'PY'
import csv, sys, statistics as st, collections as c
d = c.defaultdict(list); md5 = c.defaultdict(set)
for r in csv.DictReader(open(sys.argv[1])):
    d[(r['design'], int(r['query_length']), r['mode'])].append(float(r['elapsed_s'])); md5[(r['design'], r['query'])].add(r['out_md5'])
for k in sorted(d): print(f'{k[0]:12s} {k[1]:>5} bp  {k[2]:7s} median {st.median(d[k]):7.2f} s')
print('outputs identical across modes for every design and query:', all(len(v) == 1 for v in md5.values()))
PY
