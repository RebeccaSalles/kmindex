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

namespace kmq {

  // KMINDEX_PREFETCH=D: while row k of a partition is copied, ask the kernel for row k + D
  // (unset or 0 = original behaviour, no extra calls)
  inline std::size_t prefetch_depth()
  {
    static const std::size_t d = [] { const char* e = std::getenv("KMINDEX_PREFETCH"); return e ? std::strtoul(e, nullptr, 10) : 0UL; }();
    return d;
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

      void solve_batch(batch_query& bq)
      {
        if (m_cache)
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
