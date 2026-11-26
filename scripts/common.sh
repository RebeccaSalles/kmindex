# Shared settings of the two timing scripts (fullwidth_timing.sh, realindex_timing.sh). Every value can be set from
# the environment before running them. Sourced, not run.

# --- your two kmindex builds (README §1): v0.7.1 untouched, and v0.7.1 + PR 1 to 4
KM_REF=${KM_REF:?set KM_REF to the kmindex binary built from v0.7.1, e.g. KM_REF=$PWD/build_ref/app/kmindex/kmindex}
KM_PR=${KM_PR:?set KM_PR to the kmindex binary built from branch pr4-one-copy-of-names}
# --- where results go (a folder you can write; created) and node-local scratch for the index copy (not tmpfs/NFS)
OUT=${OUT:?set OUT to a results folder you can write, e.g. OUT=\$HOME/kmindex_evidence_results}
LOCAL=${LOCAL:-${TMPDIR:-/tmp}}
# --- what to run
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
QUERIES=${QUERIES:-$HERE/../queries}               # queries/<length>/*.fa, 5 real sequences per length
LENGTHS=${LENGTHS:-"250 500 1000"}
MODES=${MODES:-"ref pr1 pr12 pr123 pr1234"}       # ref = v0.7.1; prN = PR build with switches of PRs 1..N
THREADS=${THREADS:-4}; Z=${Z:-6}; R=${R:-0.99}    # the query settings of every run in the decks
COLD_MAX=$(( 16 * 1024 * 1024 ))                  # a cold run counts as verified at <= 16 MiB still resident

mkdir -p "$OUT"
for b in "$KM_REF" "$KM_PR"; do [[ -x $b ]] || { echo "not an executable kmindex: $b" >&2; exit 2; }; done
case "$(df --output=fstype "$LOCAL" | tail -1 | tr -d ' ')" in
  tmpfs|nfs|nfs4|ramfs) echo "LOCAL=$LOCAL is not a local disk: evicting its pages would not make queries cold" >&2; exit 3 ;;
esac

setmode() {   # sets BIN and the switches for one mode; every switch is reset first
  unset KMINDEX_PREFETCH KMINDEX_FASTCOUNT KMINDEX_COUNT_THREADS KMINDEX_STREAM KMINDEX_LEAN
  BIN=$KM_PR
  case "$1" in
    ref)    BIN=$KM_REF ;;
    pr1)    export KMINDEX_PREFETCH=16 ;;
    pr12)   export KMINDEX_PREFETCH=16 KMINDEX_FASTCOUNT=1 KMINDEX_COUNT_THREADS=8 ;;
    pr123)  export KMINDEX_PREFETCH=16 KMINDEX_FASTCOUNT=1 KMINDEX_COUNT_THREADS=8 KMINDEX_STREAM=32 ;;
    pr1234) export KMINDEX_PREFETCH=16 KMINDEX_FASTCOUNT=1 KMINDEX_COUNT_THREADS=8 KMINDEX_STREAM=32 KMINDEX_LEAN=1 ;;
    *) echo "unknown mode $1" >&2; exit 4 ;;
  esac
}
evict() {     # drop this tree's pages from the page cache: no root, no node-wide setting
  sync
  python3 - "$1" <<'PY'
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
resident() {  # bytes of the tree still in the page cache (fincore, util-linux); NA if fincore is missing
  command -v fincore >/dev/null || { echo NA; return; }
  find "$1" -type f -print0 | xargs -0 -r fincore -b -n -o RES 2>/dev/null | awk '{s+=$1} END{print s+0}'
}
state_of() {  # cold (verified) | cold_unverified
  [[ $1 != NA && $1 -le $COLD_MAX ]] && echo cold || echo cold_unverified
}
run1() {      # run1 <index dir> <query> <label for the CSV> <csv>: one query, timed, output md5 recorded
  local idx=$1 q=$2 lab=$3 csv=$4 o=$LOCAL/kmi_ev_out_$$ t=$LOCAL/kmi_ev_time_$$
  rm -rf "$o"
  /usr/bin/time -f '%e,%U,%S,%M' -o "$t" "$BIN" query -i "$idx" -q "$q" -o "$o" -t "$THREADS" -z "$Z" -r "$R" > /dev/null 2> "$OUT/last.stderr" \
    || { echo "query failed ($lab), see $OUT/last.stderr" >&2; exit 5; }
  echo "$lab,$(cat "$t"),$(find "$o" -name '*.json' -exec cat {} + | md5sum | cut -c1-12)" >> "$csv"
  rm -rf "$o" "$t"
}
