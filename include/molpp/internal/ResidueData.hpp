#ifndef MOLPP_INTERNAL_RESIDUEDATA_HPP
#define MOLPP_INTERNAL_RESIDUEDATA_HPP

#include <vector>
#include <string>

namespace mol::internal
{

class ResidueData
{
public:
    size_t size() const
    {
        return m_id.size();
    }

    int& id(size_t const index)
    {
        return m_id[index];
    }

    int const& id(size_t const index) const
    {
        return m_id[index];
    }

    std::string& name(size_t const index)
    {
        return m_name[index];
    }

    std::string const& name(size_t const index) const
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
        m_id.resize(size, -1);
        m_name.resize(size);
        m_segid.resize(size);
        m_chain.resize(size);
    }

private:
    std::vector<int> m_id;
    std::vector<std::string> m_name;
    std::vector<std::string> m_segid;
    std::vector<std::string> m_chain;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_RESIDUEDATA_HPP
