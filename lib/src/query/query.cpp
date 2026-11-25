#include <kmindex/query/query.hpp>
#include <kmindex/utils.hpp>

namespace kmq {

  query_response::query_response(const std::string& name, std::size_t n, std::size_t nbits, std::size_t width)
    : m_name(std::move(name)), m_block_size(((nbits * width) + 7) / 8), m_n(n)
  {
    m_streamed = stream_chunk() > 0 && width == 1;  // no matrix, see kindex::solve_stream
    if (!m_streamed)
      m_responses.resize(n * m_block_size, 0);
  }

  std::size_t query_response::block_size() const
  {
    return m_block_size;
  }

  std::size_t query_response::nbk() const
  {
    return m_n;  // was m_responses.size() / m_block_size, equal to n while allocated
  }

  const std::string& query_response::name() const
  {
    return m_name;
  }

  std::uint8_t* query_response::get(std::size_t mer_pos)
  {
    return m_responses.data() + (mer_pos * m_block_size);
  }

  void query_response::free()
  {
    free_container(m_responses);
  }

}
