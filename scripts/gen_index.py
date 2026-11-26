#!/usr/bin/env python3
"""Synthetic full-width kmindex index (one group), real kmindex v0.6.1 file format.

Purpose: time the real kmindex binary end to end on rows as wide as full Logan 1.2
(N samples per row, N/8 bytes), which no built index of ours reaches. Only the row COUNT is
small (R rows per partition); the row WIDTH, the number of partitions, the file layout, the
index.json (with N sample names) and every query code path are real.

Layout written under ROOT:
  index.json                      {"index": {NAME: {...}}, "path": ROOT}
  NAME -> kmindex_data/NAME       (symlink, as kmindex register makes)
  kmindex_data/NAME/hash.info     kmtricks HashWindow: u64 bloom, u64 P, u64 R, u64 R/8, u32 minim
  kmindex_data/NAME/repartition_gatb/repartition.minimRepart
                                  u16 P, u64 4^m, u16 nb_pass, u16[4^m], bool has_freq, u32 magic
  kmindex_data/NAME/matrices/matrix_p.cmbf
                                  49-byte header (copied from a real matrix; kmindex skips it)
                                  + R rows of ceil(N/8) bytes, bits set with probability DENSITY
Row bits: random, P(bit) = density. kmhelpers sizes a group's Bloom filters for the largest
possible sample of its top span (m = n_max x -ln(FP) / ln(2)^2, FP = 0.25), so no sample can fill
more than 1 - exp(-ln(2)^2 / -ln(0.25)) = 0.293 of its bits; density 0.3125 (a & (b | (c & d)))
puts EVERY sample above that bound (pessimistic: real rows are much sparser for the many small
samples of a wide group). --hit-frac plants true hits: that fraction of samples (random, fixed
seed) has every bit set, so every query matches them (pessimistic case: a very common sequence;
exercises the output path and defeats the fast loop's early exit). The original counting loop does
not depend on density.
"""
import argparse, json, os, struct, sys
from multiprocessing import Pool
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("root"); ap.add_argument("--name", default="fw_g0")
ap.add_argument("--samples", type=int, default=26857911)
ap.add_argument("--partitions", type=int, default=256)
ap.add_argument("--rows", type=int, default=256, help="rows per partition (multiple of 64)")
ap.add_argument("--density", choices=["0.25", "0.3125", "0.5"], default="0.3125")
ap.add_argument("--hit-frac", type=float, default=0.0)
ap.add_argument("--minim", type=int, default=10)
ap.add_argument("--smer", type=int, default=25)
ap.add_argument("--header-from", required=True, help="a real matrix_*.cmbf (first 49 bytes copied)")
ap.add_argument("--repart-from", required=True, help="a real repartition.minimRepart (nb_pass copied)")
ap.add_argument("--procs", type=int, default=8)
ap.add_argument("--seed", type=int, default=20261113)
a = ap.parse_args()

assert a.rows % 64 == 0
N, P, R = a.samples, a.partitions, a.rows
B = (N + 7) // 8
G = os.path.join(a.root, "kmindex_data", a.name)
os.makedirs(os.path.join(G, "matrices"), exist_ok=True)
os.makedirs(os.path.join(G, "repartition_gatb"), exist_ok=True)

with open(a.header_from, "rb") as f:
    HDR = f.read(49)
assert len(HDR) == 49

# hash.info (kmtricks HashWindow): bloom = P * R, window bits = R
with open(os.path.join(G, "hash.info"), "wb") as f:
    f.write(struct.pack("<QQQQI", P * R, P, R, R // 8, a.minim))

# repartition: same field layout as kmtricks; minimizer -> partition uniform at random
with open(a.repart_from, "rb") as f:
    rp, rn = struct.unpack("<HQ", f.read(10)); (npass,) = struct.unpack("<H", f.read(2))
nm = 4 ** a.minim
assert rn == nm, (rn, nm)
rng = np.random.default_rng(a.seed)
table = rng.integers(0, P, size=nm, dtype=np.uint16)
with open(os.path.join(G, "repartition_gatb", "repartition.minimRepart"), "wb") as f:
    f.write(struct.pack("<HQH", P, nm, npass)); f.write(table.tobytes())
    f.write(struct.pack("<?I", False, 0x12345678))


HITMASK = None
if a.hit_frac > 0:
    hits = np.random.default_rng((a.seed, 999)).random(N) < a.hit_frac
    HITMASK = np.packbits(np.concatenate([hits, np.zeros(B * 8 - N, bool)]), bitorder="little")
    print(f"planted hits: {int(hits.sum())} samples ({a.hit_frac:.4f})", flush=True)


def write_partition(p):
    path = os.path.join(G, "matrices", f"matrix_{p}.cmbf")
    r = np.random.default_rng((a.seed, p))
    chunk = 16
    with open(path, "wb") as f:
        f.write(HDR)
        for r0 in range(0, R, chunk):
            n = min(chunk, R - r0) * B
            x = np.frombuffer(r.bytes(n), dtype=np.uint8)
            u = lambda: np.frombuffer(r.bytes(n), dtype=np.uint8)
            if a.density == "0.25":
                x = x & u()
            elif a.density == "0.3125":
                x = x & (u() | (u() & u()))
            if HITMASK is not None:
                x = (x.reshape(-1, B) | HITMASK).reshape(-1)
            f.write(x.tobytes())
    return p


with Pool(a.procs) as pool:
    for i, _ in enumerate(pool.imap_unordered(write_partition, range(P))):
        if i % 32 == 31:
            print(f"partitions written: {i + 1}/{P}", flush=True)

# index.json with N sample names, as kmindex register writes it
names = [f"SRR{i:09d}" for i in range(N)]  # 12 chars, at least as long as real accessions
meta = {"bloom_size": P * R, "bw": 1, "index_size": P * R * B, "kmindex_version": "0.6.1",
        "kmtricks_version": "1.6.0", "minim_size": a.minim, "nb_partitions": P, "nb_samples": N,
        "samples": names, "sha1": "0" * 40, "smer_size": a.smer}
with open(os.path.join(a.root, "index.json"), "w") as f:
    json.dump({"index": {a.name: meta}, "path": os.path.abspath(a.root)}, f)
link = os.path.join(a.root, a.name)
if not os.path.islink(link):
    os.symlink(os.path.abspath(G), link)
print(f"done: N={N} P={P} R={R} row_bytes={B} file_bytes={49 + R * B} total_GB={P * (49 + R * B) / 1e9:.1f}")
