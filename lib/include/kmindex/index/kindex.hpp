#ifndef INDEX_HPP_FJYOTLJN
#define INDEX_HPP_FJYOTLJN

#include <memory>
#include <kmindex/query/query_results.hpp>
#include <kmindex/index/index_infos.hpp>
#include <kmindex/spinlock.hpp>
#include <mio/mmap.hpp>

#ifdef KMINDEX_WITH_COMPRESSION
#include <block_compressor/block_decompressor.hpp>
#include <block_compressor/config.hpp>
#include <block_compressor/decompressor.hpp>
#include <block_compressor/int_container.hpp>
#endif

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <thread>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <kmindex/exceptions.hpp>

namespace kmq {

  // KMINDEX_PREFETCH=D: while row k of a partition is copied, ask the kernel for row k + D
  // (unset or 0 = original behaviour, no extra calls)
  inline std::size_t prefetch_depth()
  {
    static const std::size_t d = [] { const char* e = std::getenv("KMINDEX_PREFETCH"); return e ? std::strtoul(e, nullptr, 10) : 0UL; }();
    return d;
  }

  // counting threads for the streamed path (same switch as the fast count)
  inline std::size_t stream_count_threads()
  {
    static const std::size_t t = [] { const char* e = std::getenv("KMINDEX_COUNT_THREADS"); std::size_t v = e ? std::strtoul(e, nullptr, 10) : 1UL; return v ? v : 1UL; }();
    return t;
  }

  // POSIX_FADV_WILLNEED on one row of an open (unmapped) partition file, in 1 MiB pieces
  // (posix_fadvise does not exist on macOS: there the read-ahead is skipped, the rows are read the same way)
  inline void prefetch_file_row(int fd, std::size_t bytes, std::uint64_t pos)
  {
#ifdef POSIX_FADV_WILLNEED
    const std::size_t off = (bytes * pos) + 49, end = off + bytes, piece = 1UL << 20;
    for (std::size_t c = off; c < end; c += piece)
      posix_fadvise(fd, c, std::min(piece, end - c), POSIX_FADV_WILLNEED);
#else
    (void)fd; (void)bytes; (void)pos;
#endif
  }

  class partition_interface
  {
    public:
      virtual ~partition_interface() = default;
      virtual void query(std::uint64_t pos, std::uint8_t* dest) = 0;
      virtual void prefetch(std::uint64_t pos) { (void)pos; }
  };

  class partition : public partition_interface
  {
    public:
      partition(const std::string& matrix_path, std::size_t nb_samples, std::size_t width);

      ~partition();

      virtual void query(std::uint64_t pos, std::uint8_t* dest);
      virtual void prefetch(std::uint64_t pos);

    private:
      int m_fd {0};
      mio::mmap_source m_mapped;
      std::size_t m_nb_samples {0};
      std::size_t m_bytes {0};
  };

#ifdef KMINDEX_WITH_COMPRESSION
  class compressed_partition : public partition_interface
  {
    public:
      compressed_partition(const std::string& matrix_path, const std::string& config_path, std::size_t nb_samples, std::size_t width);

      ~compressed_partition();

      virtual void query(std::uint64_t pos, std::uint8_t* dest);

    private:
      std::unique_ptr<block_compressor::BlockDecompressor> m_ptr_bd;
      std::unique_ptr<block_compressor::Decompressor> m_ptr_decompressor;
      std::unique_ptr<block_compressor::IntContainer<std::uint64_t>> m_ptr_int_container;
      std::size_t m_nb_samples {0};
      std::size_t m_bytes {0};
  };
#endif

  class kindex
  {
    public:

      kindex();
      ~kindex();
      kindex(const index_infos& i, bool cache = false);

      void init(std::size_t p);
      void unmap(std::size_t p);

    public:
      std::string name() const;
      std::string directory() const;

      void solve(batch_query& bq)
      {
        std::vector<std::size_t> order; order.reserve(m_infos.nb_partitions());
        for (std::size_t p = 0; p < m_infos.nb_partitions(); p++)
          order.push_back(p);
#ifndef __APPLE__
        static std::mt19937 g{std::random_device{}()};
        std::shuffle(std::begin(order), std::end(order), g);
#endif
        for (auto const& p : order)
          solve_one(bq, p);
      }

      void solve_one(batch_query& bq, std::size_t p)
      {
        auto& smers = bq.partition(p);
        auto& responses = bq.response();

        std::sort(std::begin(smers), std::end(smers));

        std::unique_lock<spinlock> lock(m_mutexes[p]);
        init(p);
        const std::size_t D = prefetch_depth();
        if (!D)
        {
          for (auto& [mer, qid] : smers)
          {
            m_partitions[p]->query(mer.h, responses[qid]->get(mer.i));
          }
          m_partitions[p] = nullptr;
          return;
        }
        // rows are requested D ahead of the one being copied
        const std::size_t n = smers.size();
        for (std::size_t t = 0; t < D && t < n; ++t) m_partitions[p]->prefetch(smers[t].first.h);
        for (std::size_t k = 0; k < n; ++k)
        {
          if (k + D < n) m_partitions[p]->prefetch(smers[k + D].first.h);
          m_partitions[p]->query(smers[k].first.h, responses[smers[k].second]->get(smers[k].first.i));
        }
        m_partitions[p] = nullptr;
      }

      void solve_cache(batch_query& bq)
      {
        std::vector<std::size_t> order; order.reserve(m_infos.nb_partitions());
        for (std::size_t p = 0; p < m_infos.nb_partitions(); p++)
          order.push_back(p);
#ifndef __APPLE__
        static std::mt19937 g{std::random_device{}()};
        std::shuffle(std::begin(order), std::end(order), g);
#endif
        for (auto const& p : order)
          solve_one_cache(bq, p);
      }

      void solve_one_cache(batch_query& bq, std::size_t p)
      {
        auto& smers = bq.partition(p);
        auto& responses = bq.response();

        std::sort(std::begin(smers), std::end(smers));

        std::unique_lock<spinlock> lock(m_mutexes[p]);
        for (auto& [mer, qid] : smers)
        {
          m_partitions[p]->query(mer.h, responses[qid]->get(mer.i));
        }
      }

      // KMINDEX_STREAM=C: each query's rows are read in query order with pread into a ring of C + z rows
      // (C k-mers per step), D = KMINDEX_PREFETCH rows ahead (fadvise WILLNEED, any partition), and each
      // step's k-mers are counted at once (AND of the z + 1 rows, 64 samples per word, set bits only,
      // KMINDEX_COUNT_THREADS word ranges). Counts are the same as the original loop's; only the order
      // of the reads and the buffer change. A partition file is opened at its first row in a query
      // sequence and closed right after its last one.
      void solve_stream(batch_query& bq)
      {
        if (m_infos.is_compressed_index())
          throw kmq_error("KMINDEX_STREAM does not support compressed indexes");
        auto& responses = bq.response();
        const std::size_t np = m_infos.nb_partitions(), z = bq.z(), N = m_infos.nb_samples();
        const std::uint32_t none = static_cast<std::uint32_t>(-1);
        std::vector<std::vector<std::pair<std::uint32_t, std::uint64_t>>> where(responses.size());
        for (std::size_t q = 0; q < responses.size(); ++q)
          where[q].assign(responses[q]->nbk(), {none, 0});
        for (std::size_t p = 0; p < np; ++p)
          for (auto& [mer, qid] : bq.partition(p))
            where[qid][mer.i] = {static_cast<std::uint32_t>(p), mer.h};
        std::vector<std::size_t> last(np, static_cast<std::size_t>(-1));   // last row of each file in this sequence
        std::vector<int> fds(np, -1);
        auto shut = [&](std::uint32_t p) { if (fds[p] >= 0) { close(fds[p]); fds[p] = -1; } };
        auto fd = [&](std::uint32_t p) {
          if (fds[p] < 0 && (fds[p] = open(m_infos.get_partition(p).c_str(), O_RDONLY)) < 0)
            throw kmq_error("cannot open partition " + m_infos.get_partition(p));
          return fds[p];
        };
        const std::size_t C0 = stream_chunk(), D = prefetch_depth(), T0 = stream_count_threads();
        for (std::size_t q = 0; q < responses.size(); ++q)
        {
          auto& wq = where[q];
          const std::size_t n = wq.size(), nbk = n > z ? n - z : 0, bs = responses[q]->block_size();
          std::fill(last.begin(), last.end(), static_cast<std::size_t>(-1));
          for (std::size_t t = 0; t < n; ++t)
            if (wq[t].first != none) last[wq[t].first] = t;
          std::vector<std::uint32_t> counts(N, 0);
          if (nbk == 0) { responses[q]->set_counts(std::move(counts)); continue; }
          const std::size_t C = std::min(C0, nbk), R = C + z;   // ring never larger than the query needs
          std::vector<std::uint8_t> ring(R * bs);
          auto ahead = [&](std::size_t t) {
            if (D && t < n && wq[t].first != none) prefetch_file_row(fd(wq[t].first), bs, wq[t].second);
          };
          auto fetch = [&](std::size_t t) {
            std::uint8_t* dst = ring.data() + (t % R) * bs;
            if (wq[t].first == none) { std::memset(dst, 0, bs); return; }
            const int f = fd(wq[t].first);
            const off_t off = static_cast<off_t>(bs * wq[t].second + 49);
            for (std::size_t got = 0; got < bs; )
            {
              const ssize_t r = pread(f, dst + got, bs - got, off + static_cast<off_t>(got));
              if (r <= 0) throw kmq_error("short read in partition " + m_infos.get_partition(wq[t].first));
              got += static_cast<std::size_t>(r);
            }
            if (last[wq[t].first] == t) shut(wq[t].first);
          };
          // threads only for wide rows (>= 4096 words, 262k samples, per thread): the counting step runs
          // once per C k-mers, so starting threads for a narrow group would cost more than it saves
          const std::size_t words = (bs + 7) / 8, T = std::max<std::size_t>(1, std::min(T0, words / 4096));
          std::uint32_t* cnt = counts.data(); const std::uint8_t* rg = ring.data();
          auto count = [&](std::size_t a, std::size_t b) {
            auto work = [=](std::size_t w0, std::size_t w1) {
              std::vector<const std::uint8_t*> rp(z + 1);   // the z + 1 rows of k-mer i, located once
              for (std::size_t i = a; i < b; ++i)
              {
                for (std::size_t r = 0; r <= z; ++r) rp[r] = rg + ((i + r) % R) * bs;
                for (std::size_t w = w0; w < w1; ++w)
                {
                  const std::size_t o = w * 8, nb = std::min<std::size_t>(8, bs - o);
                  std::uint64_t acc = ~0ULL;
                  for (std::size_t r = 0; r <= z; ++r)
                  {
                    std::uint64_t v = 0; std::memcpy(&v, rp[r] + o, nb); acc &= v;
                  }
                  while (acc)
                  {
                    const std::size_t s = o * 8 + static_cast<std::size_t>(__builtin_ctzll(acc));
                    if (s < N) ++cnt[s];
                    acc &= acc - 1;
                  }
                }
              }
            };
            if (T <= 1) { work(0, words); return; }
            std::vector<std::thread> th;
            for (std::size_t k = 0; k < T; ++k) th.emplace_back(work, words * k / T, words * (k + 1) / T);
            for (auto& x : th) x.join();
          };
          for (std::size_t t = 0; t < D; ++t) ahead(t);
          std::size_t next = 0;
          for (std::size_t a = 0; a < nbk; a += C)
          {
            const std::size_t b = std::min(a + C, nbk);   // k-mers [a, b) need s-mers [a, b + z)
            for (; next < b + z; ++next) { ahead(next + D); fetch(next); }
            count(a, b);
          }
          responses[q]->set_counts(std::move(counts));
          for (std::uint32_t p = 0; p < np; ++p) shut(p);   // nothing left open between sequences
        }
      }

      void solve_batch(batch_query& bq)
      {
        if (!bq.response().empty() && bq.response()[0]->streamed())
          solve_stream(bq);
        else if (m_cache)
          solve_cache(bq);
        else
          solve(bq);
      }

      index_infos& infos();
    private:
      index_infos m_infos;
      std::vector<std::unique_ptr<partition_interface>> m_partitions;
      std::vector<spinlock> m_mutexes;
      bool m_cache {false};
  };
}



#endif /* end of include guard: INDEX_HPP_FJYOTLJN */
