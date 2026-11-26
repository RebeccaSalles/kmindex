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

These are SLURM scripts for our cluster. The paths are ours, but the protocol is in the scripts.

- **The node and its disk:** genogpu2, an exclusive node. The index sits on its local SSD: a Dell BOSS card, two M.2 SATA SSDs mirrored. `lsblk` shows ROTA=1 because the card is a RAID controller.
- **Cold queries:** before each one, only the index's own files are evicted (`posix_fadvise` DONTNEED on each file). No disk-wide or node-wide setting is touched. Residency is checked with `fincore`; a run counts as cold only below 16 MiB resident.
- **Order of the modes:** the modes rotate from query to query, so no mode always runs first.
- **Predictions first:** the expected values were written down before each job was submitted (`results/PREDICTIONS_*.txt`).

**Full width** (`scripts/pr071_steps.sbatch`, job 17746813). It runs on a synthetic one-group index (`scripts/gen_index.py`) with real kmindex files:
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

**Real indexes** (`scripts/pr071_60k.sbatch`, job 17746814). These are 16 designs built with kmhelpers on a 60k-sample and a 90k-sample subset of Logan. The queries are 5 per length, cold.
- **The modes:** today; PR 1; PR 1+2; PR 1+2+3.
- **Speedup at 1000 bp (PR 1+2 over today):** median ×13.4 on the 60k sample and ×12.1 on the 90k sample; over all designs, from ×4.0 to ×15.6.
- **Outputs:** `scripts/compare_modes.py` compared every mode against today, as parsed JSON: 1,440 comparisons, 1,440 identical.

## What this does not show

- **Other machines:** times are from one machine (genogpu2) and one kind of disk (a local SATA SSD). The target, TACC Horizon, has network flash storage and Arm CPUs, so the times there must be measured again.
- **Index builds:** none of the PRs changes the index format or the build. We have not built a production-size index with the PR code: the full-width index is synthetic, and it has one group.
- **Real full-width data:** at full width the speedups come from synthetic rows. The 31% density was chosen above what kmhelpers' sizing allows, so it is meant as a pessimistic case. The fast count's cost depends on how many bits survive the AND, while the original loop costs the same whatever the density. We have not measured real rows as wide as full Logan.
