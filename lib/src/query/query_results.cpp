#include <kmindex/query/query_results.hpp>
#include <kmindex/index/index_infos.hpp>
#include <kmindex/query/format.hpp>
#include <kmindex/exceptions.hpp>

#include <bitpacker/bitpacker.hpp>
#include <nonstd/span.hpp>

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <vector>

namespace kmq {

  // KMINDEX_FASTCOUNT=1: count per-sample hits by ANDing the z+1 rows 64 samples (one 64-bit word) at a time
  // and visiting only the set bits of the result,
  // KMINDEX_COUNT_THREADS=n: split the samples into n ranges, one thread each.
  // Unset = original loop.
  static bool fast_count()
  {
    static const bool f = [] { const char* e = std::getenv("KMINDEX_FASTCOUNT"); return e && e[0] == '1'; }();
    return f;
  }

  static std::size_t count_threads()
  {
    static const std::size_t t = [] { const char* e = std::getenv("KMINDEX_COUNT_THREADS"); std::size_t v = e ? std::strtoul(e, nullptr, 10) : 1UL; return v ? v : 1UL; }();
    return t;
  }

  query_result::query_result(query_response_t&& qr, std::size_t z, const index_infos& infos, bool pos)
    : m_qr(std::move(qr)), m_z(z), m_infos(infos), m_nbk(m_qr->nbk() - z)
  {
    m_ratios.resize(m_infos.nb_samples(), 0);
    m_counts.resize(m_infos.nb_samples(), 0);

    if (pos)
      m_positions.resize(m_infos.nb_samples());

    if (m_infos.bw() > 1)
    {
      if (pos)
        compute_abs_pos();
      else
        compute_abs();
    }
    else
    {
      if (pos)
        compute_ratios_pos();
      else
        compute_ratios();
    }
  }

  void query_result::compute_ratios_fast()
  {
    const std::uint8_t* data = m_qr->get(0);
    const std::size_t bs = m_qr->block_size(), words = (bs + 7) / 8, N = m_counts.size(), nbk = m_nbk, z = m_z;
    std::uint32_t* counts = m_counts.data();
    auto work = [=](std::size_t w0, std::size_t w1) {
      for (std::size_t i = 0; i < nbk; ++i)
      {
        const std::uint8_t* row = data + i * bs;
        for (std::size_t w = w0; w < w1; ++w)
        {
          const std::size_t off = w * 8, n = std::min<std::size_t>(8, bs - off);
          std::uint64_t acc = ~0ULL;
          for (std::size_t r = 0; r <= z; ++r)
          {
            std::uint64_t v = 0; std::memcpy(&v, row + r * bs + off, n); acc &= v;
          }
          while (acc)
          {
            const std::size_t s = off * 8 + static_cast<std::size_t>(__builtin_ctzll(acc));
            if (s < N) ++counts[s];
            acc &= acc - 1;
          }
        }
      }
    };
    const std::size_t T = std::min(count_threads(), words ? words : 1);
    if (T <= 1) work(0, words);
    else
    {
      std::vector<std::thread> th;
      for (std::size_t k = 0; k < T; ++k) th.emplace_back(work, words * k / T, words * (k + 1) / T);
      for (auto& x : th) x.join();
    }
    for (std::size_t i = 0; i < N; ++i)
      m_ratios[i] = m_counts[i] / static_cast<double>(m_nbk);
    m_qr->free();
  }

  void query_result::compute_ratios()
  {
    if (fast_count()) { compute_ratios_fast(); return; }
    const uint8_t* data = m_qr->get(0);

    std::vector<uint8_t> kres(m_qr->block_size(), 255);
    std::vector<std::uint32_t> count(m_ratios.size(), 0);

    std::size_t block_size = m_qr->block_size();
    std::size_t block_size_z = m_qr->block_size() * m_z;

    for (std::size_t i = 0; i < m_nbk * block_size; i += block_size)
    {
      for (std::size_t j = i; j <= i + block_size_z; j += block_size)
      {
        for(std::size_t k = 0; k < block_size; ++k)
        {
          kres[k] &= data[j+k];
        }
      }

      for (std::size_t s = 0; s < m_ratios.size(); ++s)
      {
        m_counts[s] += static_cast<bool>(BITCHECK(kres, s));
      }

      std::fill(kres.begin(), kres.end(), 255);
    }

    for (std::size_t i = 0; i < m_ratios.size(); ++i)
    {
      m_ratios[i] = m_counts[i] / static_cast<double>(m_nbk);
    }

    m_qr->free();
  }

  void query_result::compute_ratios_pos()
  {
    const uint8_t* data = m_qr->get(0);

    std::vector<uint8_t> kres(m_qr->block_size(), 255);
    std::vector<std::uint32_t> count(m_ratios.size(), 0);

    std::size_t block_size = m_qr->block_size();
    std::size_t block_size_z = m_qr->block_size() * m_z;

    for (auto& v : m_positions)
      v.reserve(m_nbk);

    for (std::size_t i = 0; i < m_nbk * block_size; i += block_size)
    {
      for (std::size_t j = i; j <= i + block_size_z; j += block_size)
      {
        for(std::size_t k = 0; k < block_size; ++k)
        {
          kres[k] &= data[j+k];
        }
      }

      for (std::size_t s = 0; s < m_ratios.size(); ++s)
      {
        auto b = static_cast<bool>(BITCHECK(kres, s));
        m_counts[s] += static_cast<std::uint8_t>(b);
        m_positions[s].push_back(static_cast<std::uint8_t>(b));
      }

      std::fill(kres.begin(), kres.end(), 255);
    }

    for (std::size_t i = 0; i < m_ratios.size(); ++i)
    {
      m_ratios[i] = m_counts[i] / static_cast<double>(m_nbk);
    }

    m_qr->free();
  }

  void query_result::compute_abs_pos()
  {
    const uint8_t* data = m_qr->get(0);

    std::vector<std::uint32_t> kres_abs(m_counts.size(), std::numeric_limits<std::uint32_t>::max());
    std::fill(m_counts.begin(), m_counts.end(), std::numeric_limits<std::uint32_t>::min());

    std::size_t block_size = m_qr->block_size();
    std::size_t block_size_z = m_qr->block_size() * m_z;
    
    for (auto& v : m_positions)
      v.reserve(m_nbk);

    // for each k-mers
    for (std::size_t i = 0; i < m_nbk * block_size; i += block_size)
    {
      // for each s-mers in the k-mer
      for (std::size_t j = i; j <= i + block_size_z; j += block_size)
      {
        auto s = nonstd::span<const std::uint8_t>(&data[j], block_size);
        for(std::size_t k = 0, l = 0; k < m_counts.size(); ++k, l+=m_infos.bw())
        {
          kres_abs[k] = std::min(bitpacker::extract<std::uint32_t>(s, l, m_infos.bw()), kres_abs[k]);
        }
      }

      // at this point, kres is a vector of abundances of a k-mer (min of s-mers) across samples
      for (std::size_t s = 0; s < m_ratios.size(); ++s)
      {
        m_counts[s] += kres_abs[s];
        m_ratios[s] += static_cast<bool>(kres_abs[s]);
        m_positions[s].push_back(kres_abs[s]);
      }

      std::fill(kres_abs.begin(), kres_abs.end(), std::numeric_limits<std::uint32_t>::max());
    }

    for (std::size_t i = 0; i < m_ratios.size(); ++i)
    {
      m_counts[i] = m_counts[i] / m_nbk;
      m_ratios[i] = m_ratios[i] / static_cast<double>(m_nbk);
    }

  }

  void query_result::compute_abs()
  {
    const uint8_t* data = m_qr->get(0);

    std::vector<std::uint32_t> kres_abs(m_counts.size(), std::numeric_limits<std::uint32_t>::max());
    std::fill(m_counts.begin(), m_counts.end(), std::numeric_limits<std::uint32_t>::min());

    std::size_t block_size = m_qr->block_size();
    std::size_t block_size_z = m_qr->block_size() * m_z;

    // for each k-mers
    for (std::size_t i = 0; i < m_nbk * block_size; i += block_size)
    {
      // for each s-mers in the k-mer
      for (std::size_t j = i; j <= i + block_size_z; j += block_size)
      {
        auto s = nonstd::span<const std::uint8_t>(&data[j], block_size);
        for(std::size_t k = 0, l = 0; k < m_counts.size(); ++k, l+=m_infos.bw())
        {
          kres_abs[k] = std::min(bitpacker::extract<std::uint32_t>(s, l, m_infos.bw()), kres_abs[k]);
        }
      }

      // at this point, kres is a vector of abundances of a k-mer (min of s-mers) across samples
      for (std::size_t s = 0; s < m_ratios.size(); ++s)
      {
        m_counts[s] += kres_abs[s];
        m_ratios[s] += static_cast<bool>(kres_abs[s]);
      }

      std::fill(kres_abs.begin(), kres_abs.end(), std::numeric_limits<std::uint32_t>::max());
    }

    for (std::size_t i = 0; i < m_ratios.size(); ++i)
    {
      m_counts[i] = m_counts[i] / m_nbk;
      m_ratios[i] = m_ratios[i] / static_cast<double>(m_nbk);
    }
  }

  std::size_t query_result::nbk() const
  {
    return m_nbk;
  }

  const std::vector<std::uint32_t>& query_result::counts() const
  {
    return m_counts;
  }

  const std::vector<double>& query_result::ratios() const
  {
    return m_ratios;
  }

  const std::vector<std::vector<std::uint8_t>>& query_result::positions() const
  {
    return m_positions;
  }

  const std::string& query_result::name() const
  {
    return m_qr->name();
  }

  void query_result_agg::add(query_result&& r)
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    m_results.push_back(std::move(r));
  }

  std::size_t query_result_agg::size() const
  {
    return m_results.size();
  }

  const query_result_agg::vec_t& query_result_agg::results() const
  {
    return m_results;
  }

  query_result_agg::const_iterator query_result_agg::begin() const
  {
    return m_results.begin();
  }

  query_result_agg::const_iterator query_result_agg::end() const
  {
    return m_results.end();
  }

  void query_result_agg::output(const index_infos& infos,
                                const std::string& output_dir,
                                enum format f,
                                const std::string& qname,
                                double threshold,
                                const std::string& filename)
  {
    fs::create_directory(output_dir);

    const std::string path =
      filename.empty()
        ? fmt::format("{}/{}.{}", output_dir, infos.name(), format_to_fext(f))
        : fmt::format("{}/{}", output_dir, filename);
    std::ofstream out(path);
    if (!out)
      throw kmq_io_error(fmt::format("Cannot open '{}' for writing", path));

    {
      auto formatter = make_formatter(f, threshold, infos.bw());

      if (qname.size() > 0)
      {
        formatter->merge_format(infos, qname, m_results, out);
      }
      else
      {
        formatter->write_headers(out, infos);
        for (auto& r : m_results)
        {
          formatter->format(infos, r, out);
        }
      }
      // The json formatters emit their accumulated document from their
      // destructor; they must go out of scope before the close check below
      // so a failed final write is still reported.
    }
    out.close();
    if (!out)
      throw kmq_io_error(fmt::format("Failed writing '{}'", path));
  }

}
