# Evidence for the four kmindex PRs

This branch holds what backs the numbers in the review deck, and how to rerun it.
The PRs are RebeccaSalles/kmindex #1 to #4, on top of kmindex v0.7.1 (b4a07f6).

```
git clone -b evidence https://github.com/RebeccaSalles/kmindex kmindex-evidence && cd kmindex-evidence
```

Run every command below from that folder.
Each optimization sits behind its own environment switch; with no switch set, the code path is the one in v0.7.1.

| switch | PR | what it turns on |
|---|---|---|
| `KMINDEX_PREFETCH=16` | 1 | whole-row reads: `posix_madvise` WILLNEED on each row, 16 rows ahead, in 1 MiB pieces |
| `KMINDEX_FASTCOUNT=1`, `KMINDEX_COUNT_THREADS=8` | 2 | fast count: AND on 64-bit words, visit only the set bits, threads over word ranges |
| `KMINDEX_STREAM=32` | 3 | streamed rows: `pread` in query order into a ring of 32 + z rows |
| `KMINDEX_LEAN=1` | 4 | one copy of the sample names, ratios computed when written |

There are three checks, from quickest to heaviest.

## 1. Identity test (any machine, about 30 s)

This compares the outputs of v0.7.1 with v0.7.1 + PR 1 to 4, under 10 switch combinations:
- 4 output formats (json, jsonl, matrix, json_vec);
- z = 0 and 3;
- 1 and 4 threads;
- the indexes bundled in `tests/data`.

json and jsonl outputs are compared as parsed records. Matrix outputs are compared as a multiset of result lines, because v0.7.1 itself splits the batches differently from one run to the next (one header per batch file).

```
git clone https://github.com/tlemane/kmindex ref && git -C ref checkout b4a07f6 && git -C ref submodule update --init --recursive
git clone --recursive -b pr4-one-copy-of-names https://github.com/RebeccaSalles/kmindex new   # stacked branches: PR 4 contains PR 1 to 3
for d in ref new; do cmake -S $d -B build_$d -DCMAKE_BUILD_TYPE=Release -DWITH_SERVER=OFF -DWITH_COMPRESSION=OFF -DPORTABLE_BUILD=ON -DWITH_TESTS=OFF && make -C build_$d -j4; done
SRC=$PWD/ref bash scripts/identity_test.sh build_ref/app/kmindex/kmindex build_new/app/kmindex/kmindex /tmp/idtest
```

Result on 25 Nov 2026 (Ubuntu under WSL2, g++, local builds): `TOTAL 480 comparisons, 0 different or failed`.

An optional fourth argument adds two indexes built with kmtricks. They are not in this folder:
- the same data with 4096 partitions;
- query files with 1 and 20 sequences of 1000 bp, z = 0 and 6.

That brings the total to 600 comparisons; it also gave 0 different.

## 2. Counting benchmark (one C++ file)

`scripts/count_bench2.cpp` contains:
- **the original loop:** the v0.7.1 `compute_ratios` loop, copied verbatim;
- **the word versions:** the 64-bit word versions, on 1 and 8 threads.

It runs them on synthetic rows as wide as full Logan (26,857,911 samples, 976 k-mers, z = 6, 31% of bits set, 0.1% of samples with every bit set) and on rows of a real 60k-sample partition file. The total count must be the same for every variant.

```
g++ -O3 -DNDEBUG -std=c++17 -pthread -o count_bench2 scripts/count_bench2.cpp
./count_bench2 syn 26857911 976 0.3125 0.001 base      # also: words, wordsT8, sparse, sparseT8, stats
```

It needs about 3.3 GB of RAM for the rows. Results on genogpu2 are in `results/count_bench_17721442.csv`.

| variant | ns per sample and k-mer |
|---|---|
| base (the original loop) | 1.08 |
| words | 0.19 |
| wordsT8 | 0.025 |
| sparse | 0.32 |
| sparseT8 | 0.044 |

`sparse` is the early stop: skip a word's remaining rows once it is all zeros. It was slower, so it is not in PR 2.

## 3. Timing campaigns (Genouest, genogpu2, local SSD)

What we ran is below, with our job scripts (`scripts/pr071_*.sbatch`: our account and paths, kept as the record). To rerun the same protocol with your own account, on Genouest or any Linux machine, use the two portable scripts of §4.

- **The node and its disk:** genogpu2, an exclusive node. The index sits on its local SSD: a Dell BOSS card, two M.2 SATA SSDs mirrored. `lsblk` shows ROTA=1 because the card is a RAID controller.
- **Cold queries:** before each one, only the index's own files are evicted (`posix_fadvise` DONTNEED on each file). No disk-wide or node-wide setting is touched. Residency is checked with `fincore`; a run counts as cold only below 16 MiB resident.
- **Order of the modes:** the modes rotate from query to query, so no mode always runs first.
- **Predictions first:** the expected values were written down before each job was submitted (`results/PREDICTIONS_*.txt`).

**Full width** (`scripts/pr071_steps.sbatch`, job 17746813; portable version: §4). It runs on a synthetic one-group index (`scripts/gen_index.py`) with real kmindex files:
- 256 partitions, 128 rows each, every row 3.4 MB (26.9 M samples);
- an index.json with 26.9 M names;
- 31% of the bits set;
- 0.1% planted samples (26,728) with every bit set, so every query has exactly those answers.

It needs about 110 GB of local disk. The queries are 5 per length (250, 500 and 1000 bp), cold, `-t 4 -z 6 -r 0.99`. One extra run at 250 bp uses `-r 0`, so all 26.9 M samples are written out.

Cold medians at 1000 bp (`results/fullwidth_17746813_timings.csv`):

| mode | switches | time | max RSS |
|---|---|---|---|
| ref | v0.7.1 | 54.04 s | 6.30 GB |
| pr1 | 1 | 36.54 s | 6.30 GB |
| pr12 | 1, 2 | 9.02 s | 6.19 GB |
| pr123 | 1, 2, 3 | 6.77 s | 2.91 GB |
| pr1234 | 1, 2, 3, 4 | 5.91 s | 1.84 GB |

The output md5 is the same in every mode for every query, including the run with all samples written out.

A note on the binary:
- **What was timed:** the PR modes ran a build of PR 1 to 4 plus one measurement-only patch, which reads the names only for the hits (switch `KMINDEX_LAZY_NAMES`, mode pr12345 in the CSV).
- **Why that is fine for the PR modes:** that patch is not in the PRs, and it does nothing unless its switch is set.
- **Where the PR-only build is checked:** the identity test above runs it.

**Real indexes** (`scripts/pr071_60k.sbatch`, job 17746814; portable version: §4). These are 16 designs built with kmhelpers on a 60k-sample and a 90k-sample subset of Logan. The queries are 5 per length, cold.
- **The modes:** today; PR 1; PR 1+2; PR 1+2+3.
- **Speedup at 1000 bp (PR 1+2 over today):** median ×13.4 on the 60k sample and ×12.1 on the 90k sample; over all designs, from ×4.0 to ×15.6.
- **Outputs:** `scripts/compare_modes.py` compared every mode against today, as parsed JSON: 1,440 comparisons, 1,440 identical.

## 4. Rerun the timings yourself

Two portable scripts take your own folders and builds; every setting is an environment variable.

- `scripts/fullwidth_timing.sh`: generates the full-width index on local disk, then times every mode (§3, full width).
- `scripts/realindex_timing.sh <index dir> ...`: copies already-built indexes to local disk, rewrites the absolute paths kmindex stores in them, then times every mode (§3, real indexes).

Both rotate the modes from query to query, evict only the index's own pages before each query, record the residency (`fincore`) and the output md5, and print the medians and whether the outputs are identical across modes.

**Requirements**

| | full width | real indexes |
|---|---|---|
| local disk (not NFS, not tmpfs) | PARTS × ROWS × 3.36 MB: **110 GB** by default (see below) | the size of each design (60k: 7 to 21 GB; 90k: 10 to 32 GB) |
| memory | 32 GB is enough (one query peaks at 6.3 GB; the generator holds the 26.9 M names) | 32 GB |
| time on genogpu2 | about 45 min (index 5 min, then 5 modes × 15 cold queries) | about 3 h for the 32 designs of §3 |
| software | Linux, bash, python3 with numpy, GNU time (`/usr/bin/time`); `fincore` (util-linux) to verify the queries are cold, optional | the same |
| data | none: generated; its 49-byte file header and repartition format come from `tests/data` of v0.7.1 (`SRC`) | read access to the indexes (on Genouest, ours are readable by the group `prj_logan_compression`) |

**Less disk:** `ROWS` (a multiple of 64) and `PARTS` set the index size. `ROWS=64` halves it (55 GB); `ROWS=64 PARTS=128` takes 27.5 GB. Rows stay 3.36 MB wide, so each row still costs what it costs at full width. With fewer rows in all, a query hits the same row more often, which then comes from memory, so the times are a little lower. Keep the defaults to compare with our numbers.

**Settings shared by both scripts** (`scripts/common.sh`):

| variable | meaning | default |
|---|---|---|
| `KM_REF` | kmindex built from v0.7.1 (§1) | required |
| `KM_PR` | kmindex built from `pr4-one-copy-of-names` (§1) | required |
| `OUT` | results folder you can write | required |
| `LOCAL` | node-local scratch for the index | `$TMPDIR` or `/tmp` |
| `MODES` | `ref` = v0.7.1; `prN` = PR build with the switches of PRs 1..N | `ref pr1 pr12 pr123 pr1234` |
| `LENGTHS`, `QUERIES` | query lengths, and the folder of `queries/<length>/*.fa` | `250 500 1000`, `queries/` (5 real sequences per length, from Logan unitigs) |
| `THREADS`, `Z`, `R` | `kmindex query -t -z -r` | `4`, `6`, `0.99` |
| `SRC`, `SAMPLES`, `PARTS`, `ROWS`, `PROCS` | full width only: v0.7.1 checkout, index size, generator processes | required, 26857911, 256, 128, 8 |

**On Genouest** (on a compute node: the login node does not run builds). Build the two binaries as in §1 inside an `srun` session, then, for example:

```
export KM_REF=$PWD/build_ref/app/kmindex/kmindex KM_PR=$PWD/build_new/app/kmindex/kmindex SRC=$PWD/ref
export OUT=/projects/<your project>/<you>/kmindex_evidence   # or any folder you can write
sbatch -p <partition> --nodelist=<node> --exclusive --cpus-per-task=16 --mem=32G --time=03:00:00 --export=ALL \
       scripts/fullwidth_timing.sh
B=/projects/logan_compression/rsalles/span_tuning_paper
sbatch -p <partition> --nodelist=<node> --exclusive --cpus-per-task=16 --mem=32G --time=06:00:00 --export=ALL \
       scripts/realindex_timing.sh $B/training_rerun_builds/b1p1_g02 $B/training_rerun_builds/b1p05_g249
```

We used `-p gpu --nodelist=genogpu2 --exclusive` (its local disk is the SSD of §3). Both scripts were checked on a laptop on 27 Nov at small size (`SAMPLES=100000 PARTS=8 ROWS=64`): outputs identical across the five modes, the 101 planted samples found by v0.7.1 and by the PRs.

## What this does not show

- **Other machines:** times are from one machine (genogpu2) and one kind of disk (a local SATA SSD). The target, TACC Horizon, has network flash storage and Arm CPUs, so the times there must be measured again.
- **Index builds:** none of the PRs changes the index format or the build. We have not built a production-size index with the PR code: the full-width index is synthetic, and it has one group.
- **Real full-width data:** at full width the speedups come from synthetic rows. The 31% density was chosen above what kmhelpers' sizing allows, so it is meant as a pessimistic case. The fast count's cost depends on how many bits survive the AND, while the original loop costs the same whatever the density. We have not measured real rows as wide as full Logan.
