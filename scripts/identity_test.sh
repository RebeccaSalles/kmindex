#!/usr/bin/env bash
# Identity test: kmindex v0.7.1 untouched (REF) vs v0.7.1 + PR1-4 (NEW), every switch combination, several indexes,
# formats, z values and thread counts. Outputs compared as parsed JSON (json/jsonl) or as a multiset of result
# lines (matrix: upstream itself splits batches differently between runs, one header per batch file).
#
# usage: identity_test.sh REF_KMINDEX NEW_KMINDEX WORKDIR [EXTRA_DIR]
#   REF_KMINDEX  kmindex binary built from v0.7.1 (b4a07f6)
#   NEW_KMINDEX  kmindex binary built from v0.7.1 + PR 1-4
#   WORKDIR      scratch folder (created; its data/ copy of tests/data is rewritten)
#   EXTRA_DIR    optional: folder with indexes idx/ and idxP4096/ (4096 partitions) and queries q/q1_1000.fa,
#                q/multi.fa, q/multi20.fa (several sequences per file); built with kmtricks, not in the repo
# Needs: python3; the tests/data folder of the v0.7.1 checkout (SRC below).
# Result on 25 Nov 2026: 480 comparisons on the bundled indexes (+120 with EXTRA_DIR), 0 different.
set -u
REF=$(realpath "$1"); NEW=$(realpath "$2"); T=$(realpath -m "$3"); EXTRA=${4:-}
for b in "$1" "$2"; do [[ -x $b ]] || { echo "not an executable kmindex: $b (build it first)" >&2; exit 2; }; done
SRC=${SRC:?set SRC to your kmindex v0.7.1 checkout, e.g. SRC=$PWD/ref}
mkdir -p $T
rm -rf $T/data && cp -a $SRC/tests/data $T/data
python3 - $T <<'PY'
import sys, glob, os
T = sys.argv[1]
def recs(f):
    out, name, seq = [], None, []
    for l in open(f):
        l = l.strip()
        if l.startswith('>'):
            if name: out.append((name, ''.join(seq)))
            name, seq = l[1:].split()[0], []
        else: seq.append(l)
    if name: out.append((name, ''.join(seq)))
    return out
for kind in ('pa', 'abs'):
    fs = sorted(glob.glob(f'{T}/data/datasets/{kind}_dataset/*.fasta'))[:4]
    with open(f'{T}/q_{kind}.fa', 'w') as o:
        for i, f in enumerate(fs):
            for n, s in recs(f)[:2]: o.write(f'>{kind}{i}_{n}\n{s[:1000]}\n')
PY
cat > $T/cmp.py <<'PY'
import glob, json, os, sys
def load(d):
    m = {}
    for p in sorted(glob.glob(os.path.join(d, '**', '*'), recursive=True)):
        if os.path.isdir(p): continue
        rel = os.path.relpath(p, d); b = os.path.basename(p)
        if b.endswith('.json'):
            for idx, qs in json.load(open(p)).items():
                for q, v in qs.items(): m[('json', idx, q)] = v
        elif b.endswith('.jsonl'):
            for line in open(p):
                if line.strip(): d2 = json.loads(line); m[('jsonl', json.dumps(d2, sort_keys=True))] = 1
        else:   # matrix: lines of every batch file with this name, as a multiset
            m.setdefault(('text', b), []).extend(l for l in open(p).read().splitlines() if ':' in l.split('\t')[0])  # result lines only: one header per batch file
    for k in m:
        if k[0] == 'text': m[k] = sorted(m[k])
    return m
a, b = load(sys.argv[1]), load(sys.argv[2])
print('same' if a and a == b else f'DIFF({len(a)},{len(b)})')
PY
MODES="off: pr1:KMINDEX_PREFETCH=16 pr2t1:KMINDEX_FASTCOUNT=1,KMINDEX_COUNT_THREADS=1 pr2t8:KMINDEX_FASTCOUNT=1,KMINDEX_COUNT_THREADS=8 pr12:KMINDEX_PREFETCH=16,KMINDEX_FASTCOUNT=1,KMINDEX_COUNT_THREADS=8 pr3:KMINDEX_STREAM=32,KMINDEX_PREFETCH=16,KMINDEX_FASTCOUNT=1,KMINDEX_COUNT_THREADS=8 pr3c1:KMINDEX_STREAM=1 pr4:KMINDEX_LEAN=1 pr34:KMINDEX_STREAM=32,KMINDEX_PREFETCH=16,KMINDEX_FASTCOUNT=1,KMINDEX_COUNT_THREADS=8,KMINDEX_LEAN=1 pr124:KMINDEX_PREFETCH=16,KMINDEX_FASTCOUNT=1,KMINDEX_COUNT_THREADS=8,KMINDEX_LEAN=1"
tot=0; bad=0
run() {  # label idxdir indexpath query fmt z threads
  local lab=$1 dir=$2 idx=$3 q=$4 fmt=$5 z=$6 th=$7
  (cd $dir && rm -rf $T/o_ref && $REF query -i $idx -q $q -o $T/o_ref -z $z -r 0 -f $fmt -t $th > $T/ref.log 2>&1) || { echo "REF FAILED $lab $fmt z$z t$th"; tail -3 $T/ref.log; return; }
  for m in $MODES; do
    local name=${m%%:*} envs=${m#*:}
    (cd $dir && rm -rf $T/o_new && env $(echo $envs | tr ',' ' ') $NEW query -i $idx -q $q -o $T/o_new -z $z -r 0 -f $fmt -t $th > $T/new.log 2>&1) || { echo "NEW FAILED $lab $name $fmt z$z t$th: $(tail -1 $T/new.log)"; bad=$((bad+1)); tot=$((tot+1)); continue; }
    r=$(python3 $T/cmp.py $T/o_ref $T/o_new); tot=$((tot+1)); [[ $r == same ]] || { bad=$((bad+1)); echo "DIFF $lab $name $fmt z$z t$th: $r"; }
  done
}
for fmt in json jsonl matrix json_vec; do for z in 0 3; do for th in 1 4; do
  run pa_global $T/data ./indexes/pa_global $T/q_pa.fa $fmt $z $th
  run abs_global $T/data ./indexes/abs_global $T/q_abs.fa $fmt $z $th
  run registry2 $T/data ./indexes/index $T/q_pa.fa $fmt $z $th
done; done; done
if [[ -n $EXTRA ]]; then
for idx in idx idxP4096; do for q in q1_1000.fa multi.fa multi20.fa; do for z in 0 6; do
  run $idx $EXTRA $EXTRA/$idx $EXTRA/q/$q json $z 4
done; done; done
fi
echo "TOTAL $tot comparisons, $bad different or failed"
