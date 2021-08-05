#ifndef RESIDUEDATA_HPP
#define RESIDUEDATA_HPP

#include "tools/iterators.hpp"
#include <vector>
#include <string>
#include <unordered_set>

namespace mol::internal {

class ResidueData
{
private:
    using indices_type = std::unordered_set<size_t>;

public:
    using indices_iterator = indices_type::const_iterator;

    ResidueData()
    {}

    void set(size_t const index, int const res_id, std::string const res_name, std::string const seg_id, std::string const chain_id)
    {
        resid(index) = res_id;
        resname(index) = res_name;
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

    int &resid(size_t const index)
    {
        return m_resid[index];
    }

    std::string &resname(size_t const index)
    {
        return m_resname[index];
    }

    std::string &segid(size_t const index)
    {
        return m_segid[index];
    }

    std::string &chain(size_t const index)
    {
        return m_chain[index];
    }

    void resize(size_t const size)
    {
        m_indices.resize(size);
        m_resid.resize(size, -1);
        m_resname.resize(size);
        m_segid.resize(size);
        m_chain.resize(size);
    }

    Range<indices_iterator> indices(size_t const index) const
    {
        indices_type const &residue = m_indices[index];
        return {residue.cbegin(), residue.cend()};
    }

    void reset(size_t const index, size_t const new_size=0)
    {
        indices_type &residue = m_indices[index];
        residue.clear();
        residue.reserve(new_size);
    }

    void add_atom(size_t const residue, size_t const atom)
    {
        m_indices[residue].insert(atom);
    }

    void remove_atom(size_t const residue, size_t const atom)
    {
        m_indices[residue].erase(atom);
    }

private:
    std::vector<int> m_resid;
    std::vector<std::string> m_resname;
    std::vector<std::string> m_segid;
    std::vector<std::string> m_chain;
    std::vector<indices_type> m_indices;
};

} // namespace mol::internal

#endif // RESIDUEDATA_HPP
