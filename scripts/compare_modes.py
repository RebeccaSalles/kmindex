#!/usr/bin/env python3
"""[2026-11-25] Identity of outputs across modes in a requery session: for every design / length / query, the parsed
JSON of mode none (q..._L) vs each other mode (q..._L_<mode>). Prints counts."""
import glob, json, os, sys
def load(d):
    m = {}
    for p in glob.glob(os.path.join(d, '**', '*.json'), recursive=True):
        for idx, qs in json.load(open(p)).items():
            for q, v in qs.items(): m[(idx, q)] = v
    return m
tot = same = 0; bad = []
for root in sys.argv[1:]:
    for cfg in sorted(os.listdir(root)):
        for L in sorted(os.listdir(os.path.join(root, cfg))):
            base = os.path.join(root, cfg, L)
            for q in sorted(d for d in os.listdir(base) if d.count('_') == 1):
                ref = load(os.path.join(base, q))
                for other in sorted(d for d in os.listdir(base) if d.startswith(q + '_')):
                    tot += 1; ok = bool(ref) and load(os.path.join(base, other)) == ref; same += ok
                    if not ok: bad.append(f'{cfg}/{L}/{other}')
print(f'comparisons {tot}, identical {same}, different {len(bad)}'); print('\n'.join(bad[:20]))
