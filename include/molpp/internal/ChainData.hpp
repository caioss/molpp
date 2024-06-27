#ifndef MOLPP_INTERNAL_CHAINDATA_HPP
#define MOLPP_INTERNAL_CHAINDATA_HPP

#include <molpp/MolppCore.hpp>

#include <vector>
#include <string>

namespace mol::internal
{

class ChainData
{
public:
    size_t size() const
    {
        return m_name.size();
    }

    std::string& name(size_t const index)
    {
        return m_name[index];
    }

    std::string const& name(size_t const index) const
    {
        return m_name[index];
    }

    void resize(size_t const size)
    {
        m_name.resize(size);
    }

private:
    std::vector<std::string> m_name;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_CHAINDATA_HPP
