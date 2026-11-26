// [2026-11-20] Counting loop, step by step (deck: how much each part of optimization 2 gives), on genogpu2.
// base   : kmindex v0.6.1 query_result::compute_ratios, verbatim (byte AND of the z+1 rows, then one BITCHECK per sample)
// words  : AND of the z+1 rows 64 samples (one 64-bit word) at a time, ALL z+1 rows, then visit only the set bits
// sparse : same + early stop: skip the remaining rows of a word once it is all zeros (= kmindex patch, 1 thread)
// wordsT8 / sparseT8 : the same on 8 threads (samples split into 8 word ranges)
// Rows: synthetic (each bit 1 with probability p, plus a fraction h of planted samples with every bit 1, as the
// full-width index) or real (K random rows of a real kmindex partition file). Counts checked equal to base.
// Also prints, for sparse, how many of the z+1 rows a word needs on average and the share of words empty after r rows.
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <thread>
#include <vector>
#define BITSLOT(b) ((b) / 8)
#define BITMASK(b) (1 << ((b) % 8))
#define BITCHECK(a, b) ((a)[BITSLOT(b)] & BITMASK(b))
using clk = std::chrono::steady_clock;

static void base(const std::uint8_t* data, std::size_t bs, std::size_t nbk, std::size_t z, std::size_t N, std::vector<std::uint32_t>& c)
{
  std::vector<uint8_t> kres(bs, 255);
  const std::size_t bsz = bs * z;
  for (std::size_t i = 0; i < nbk * bs; i += bs)
  {
    for (std::size_t j = i; j <= i + bsz; j += bs)
      for (std::size_t k = 0; k < bs; ++k) kres[k] &= data[j + k];
    for (std::size_t s = 0; s < N; ++s) c[s] += static_cast<bool>(BITCHECK(kres, s));
    std::fill(kres.begin(), kres.end(), 255);
  }
}

template <bool EARLY>
static void words(const std::uint8_t* data, std::size_t bs, std::size_t nbk, std::size_t z, std::size_t N,
                  std::uint32_t* c, std::size_t w0, std::size_t w1)
{
  for (std::size_t i = 0; i < nbk; ++i)
  {
    const std::uint8_t* row = data + i * bs;
    for (std::size_t w = w0; w < w1; ++w)
    {
      const std::size_t off = w * 8, n = std::min<std::size_t>(8, bs - off);
      std::uint64_t acc = ~0ULL;
      for (std::size_t r = 0; r <= z && (!EARLY || acc); ++r)
      {
        std::uint64_t v = 0; std::memcpy(&v, row + r * bs + off, n); acc &= v;
      }
      while (acc)
      {
        const std::size_t s = off * 8 + static_cast<std::size_t>(__builtin_ctzll(acc));
        if (s < N) ++c[s];
        acc &= acc - 1;
      }
    }
  }
}

int main(int argc, char** argv)
{
  // synthetic: count_bench2 syn N K p h variant      real: count_bench2 real <matrix_file> N K seed variant
  const std::string src = argv[1]; const std::size_t z = 6;
  std::size_t N, K; double p = 0, h = 0; std::string variant, label;
  std::vector<std::uint8_t> data;
  if (src == "syn")
  {
    N = std::strtoull(argv[2], nullptr, 10); K = std::strtoull(argv[3], nullptr, 10); p = std::atof(argv[4]); h = std::atof(argv[5]); variant = argv[6];
    const std::size_t bs = (N + 7) / 8; data.resize(bs * K);
    std::mt19937_64 g(1); std::bernoulli_distribution bit(p), hit(h);
    for (auto& b : data) { std::uint8_t v = 0; for (int k = 0; k < 8; ++k) v |= std::uint8_t(bit(g)) << k; b = v; }
    for (std::size_t s = 0; s < N; ++s) if (hit(g)) for (std::size_t r = 0; r < K; ++r) data[r * bs + s / 8] |= std::uint8_t(1u << (s % 8));
    label = "syn_p" + std::to_string(p).substr(0, 6) + "_h" + std::to_string(h).substr(0, 6);
  }
  else
  {
    const std::string f = argv[2]; N = std::strtoull(argv[3], nullptr, 10); K = std::strtoull(argv[4], nullptr, 10);
    const unsigned seed = std::strtoul(argv[5], nullptr, 10); variant = argv[6];
    const std::size_t bs = (N + 7) / 8; data.resize(bs * K);
    std::ifstream in(f, std::ios::binary); in.seekg(0, std::ios::end); const std::size_t rows = (static_cast<std::size_t>(in.tellg()) - 49) / bs;
    std::mt19937_64 g(seed); std::uniform_int_distribution<std::size_t> pick(0, rows - 1);
    for (std::size_t r = 0; r < K; ++r) { in.seekg(49 + pick(g) * bs); in.read(reinterpret_cast<char*>(data.data() + r * bs), bs); }
    label = "real_" + f.substr(f.find_last_of('/') + 1);
  }
  const std::size_t bs = (N + 7) / 8, nbk = K - z, W = (bs + 7) / 8;
  double ones = 0; for (auto b : data) ones += __builtin_popcount(b); const double fill = ones / (8.0 * bs * K);
  std::vector<std::uint32_t> ref(N, 0), cnt(N, 0);
  double tb = -1;
  if (variant == "base" || N <= 1000000)
  { auto a = clk::now(); base(data.data(), bs, nbk, z, N, ref); tb = std::chrono::duration<double>(clk::now() - a).count(); }
  double t = tb;
  if (variant != "base" && variant != "stats")
  {
    const bool early = variant.rfind("sparse", 0) == 0;
    const std::size_t T = variant.find("T8") != std::string::npos ? 8 : 1;
    auto a = clk::now(); std::vector<std::thread> th;
    for (std::size_t k = 0; k < T; ++k)
    {
      const std::size_t w0 = W * k / T, w1 = W * (k + 1) / T;
      if (early) th.emplace_back(words<true>, data.data(), bs, nbk, z, N, cnt.data(), w0, w1);
      else th.emplace_back(words<false>, data.data(), bs, nbk, z, N, cnt.data(), w0, w1);
    }
    for (auto& x : th) x.join();
    t = std::chrono::duration<double>(clk::now() - a).count();
  }
  if (variant == "stats")   // rows needed per word with the early stop, and share of words empty after r rows
  {
    std::vector<double> empty_after(z + 2, 0); double rows_read = 0, total = 0;
    for (std::size_t i = 0; i < nbk; ++i)
      for (std::size_t w = 0; w < W; ++w)
      {
        const std::size_t off = w * 8, n = std::min<std::size_t>(8, bs - off);
        std::uint64_t acc = ~0ULL; std::size_t r = 0;
        for (; r <= z && acc; ++r) { std::uint64_t v = 0; std::memcpy(&v, data.data() + (i + r) * bs + off, n); acc &= v; if (!acc) empty_after[r + 1] += 1; }
        rows_read += r; total += 1;
      }
    std::printf("stats,%s,N=%zu,K=%zu,fill=%.4f,rows_per_word=%.3f,of=%zu", label.c_str(), N, K, fill, rows_read / total, z + 1);
    double cum = 0; for (std::size_t r = 1; r <= z + 1; ++r) { cum += empty_after[r] / total; std::printf(",empty_after_%zu=%.4f", r, cum); }
    std::printf("\n"); return 0;
  }
  const bool checked = (variant == "base") || (N <= 1000000);
  const bool same = variant == "base" || !checked || ref == cnt;
  unsigned long long hits = 0; for (auto c : (variant == "base" ? ref : cnt)) hits += c;
  std::printf("time,%s,%s,%zu,%zu,%.4f,%.6f,%.4f,%s,%llu\n", label.c_str(), variant.c_str(), N, K, fill, t, t / (double(nbk) * N) * 1e9,
              checked ? (same ? "identical" : "DIFFERENT") : "not_checked", hits);
}
