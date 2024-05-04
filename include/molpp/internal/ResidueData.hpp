#ifndef MOLPP_INTERNAL_RESIDUEDATA_HPP
#define MOLPP_INTERNAL_RESIDUEDATA_HPP

#include <molpp/tools/iterators.hpp>

#include <vector>
#include <string>
#include <ranges>
#include <unordered_set>

namespace mol::internal
{

class ResidueData
{
private:
    using indices_type = std::unordered_set<index_t>;

public:
    ResidueData()
    {}

    void set(index_t const index, int const res_id, std::string const& res_name, std::string const& seg_id, std::string const& chain_id)
    {
        residue_id(index) = res_id;
        residue_name(index) = res_name;
        segid(index) = seg_id;
        chain(index) = chain_id;
    }

    size_t size() const
    {
        return m_indices.size();
    }

    size_t size(size_t const index) const
    {
        return m_indices[index].size();
    }

    int& residue_id(size_t const index)
    {
        return m_id[index];
    }

    int const& residue_id(size_t const index) const
    {
        return m_id[index];
    }

    std::string& residue_name(size_t const index)
    {
        return m_name[index];
    }

    std::string const& residue_name(size_t const index) const
    {
        return m_name[index];
    }

    std::string& segid(size_t const index)
    {
        return m_segid[index];
    }

    std::string const& segid(size_t const index) const
    {
        return m_segid[index];
    }

    std::string& chain(size_t const index)
    {
        return m_chain[index];
    }

    std::string const& chain(size_t const index) const
    {
        return m_chain[index];
    }

    void resize(size_t const size)
    {
        m_indices.resize(size);
        m_id.resize(size, -1);
        m_name.resize(size);
        m_segid.resize(size);
        m_chain.resize(size);
    }

    auto const indices(index_t const index) const
    {
        return std::ranges::views::all(m_indices[index]);
    }

    void reset(index_t const index, size_t const new_size = 0)
    {
        indices_type& residue = m_indices[index];
        residue.clear();
        residue.reserve(new_size);
    }

    void add_atom(index_t const residue, index_t const atom)
    {
        m_indices[residue].insert(atom);
    }

    void remove_atom(index_t const residue, index_t const atom)
    {
        m_indices[residue].erase(atom);
    }

private:
    std::vector<int> m_id;
    std::vector<std::string> m_name;
    std::vector<std::string> m_segid;
    std::vector<std::string> m_chain;
    std::vector<indices_type> m_indices;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_RESIDUEDATA_HPP
