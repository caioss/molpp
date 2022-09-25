#ifndef ELEMENTSTABLE_HPP
#define ELEMENTSTABLE_HPP

#include <vector>
#include <string>
#include <initializer_list>

namespace mol {

class ElementsTable
{
public:
    struct Element
    {
        int atomic_number;
        float mass;
        float covalent_radius;
        float VDW_radius;
        std::string symbol;
        std::string name;
    };

    ElementsTable(std::initializer_list<Element> data);

    int const &atomic_number(int const atomic) const
    {
        return m_atomic[atomic];
    }

    float const &mass(int const atomic) const
    {
        return m_mass[atomic];
    }

    float const &covalent_radius(int const atomic) const
    {
        return m_covalent[atomic];
    }

    float const &VDW_radius(int const atomic) const
    {
        return m_vdw[atomic];
    }

    std::string const &symbol(int const atomic) const
    {
        return m_symbol[atomic];
    }

    std::string const &name(int const atomic) const
    {
        return m_name[atomic];
    }

    float max_VDW_radius() const
    {
        return m_max_vdw;
    }

    float max_covalent_radius() const
    {
        return m_max_covalent;
    }

    Element operator()(int const atomic) const
    {
        return {m_atomic[atomic], m_mass[atomic], m_covalent[atomic], m_vdw[atomic], m_symbol[atomic], m_name[atomic]};
    }

private:
    float m_max_vdw;
    float m_max_covalent;
    std::vector<int> m_atomic;
    std::vector<float> m_mass;
    std::vector<float> m_covalent;
    std::vector<float> m_vdw;
    std::vector<std::string> m_symbol;
    std::vector<std::string> m_name;
};

// Elements data derived from the "Blue Obelisk Element Repository"
// which is distributed under the MIT license.
ElementsTable const& ELEMENTS_TABLE();

} // namespace mol

#endif // ELEMENTSTABLE_HPP
