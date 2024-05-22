#ifndef MOLPP_INTERNAL_ATOMDATA_HPP
#define MOLPP_INTERNAL_ATOMDATA_HPP

#include <molpp/MolppCore.hpp>

#include <vector>
#include <string>
#include <optional>
#include <unordered_map>

namespace mol::internal
{

class AtomData
{
public:
    AtomData() = delete;

    AtomData(size_t const num_atoms)
    : m_num_atoms{num_atoms}
    , m_atomic(num_atoms, 0)
    , m_occupancy(num_atoms, 0)
    , m_tempfactor(num_atoms, 0)
    , m_mass(num_atoms, 0)
    , m_charge(num_atoms, 0)
    , m_radius(num_atoms, 0)
    , m_name(num_atoms)
    , m_type(num_atoms)
    , m_altloc(num_atoms)
    , m_insertion_code(num_atoms)
    {
        m_residue.reserve(num_atoms);
    }

    size_t size() const
    {
        return m_num_atoms;
    }

    std::optional<index_t> residue(index_t const atom_index) const
    {
        auto iter = m_residue.find(atom_index);
        if (iter == m_residue.end())
        {
            return std::nullopt;
        }

        return iter->second;
    }

    void set_residue(index_t const atom_index, index_t const residue_index)
    {
        m_residue[atom_index] = residue_index;
    }

    int& atomic_number(size_t const index)
    {
        return m_atomic[index];
    }

    int const& atomic_number(size_t const index) const
    {
        return m_atomic[index];
    }

    float& occupancy(size_t const index)
    {
        return m_occupancy[index];
    }

    float const& occupancy(size_t const index) const
    {
        return m_occupancy[index];
    }

    float& temperature_factor(size_t const index)
    {
        return m_tempfactor[index];
    }

    float const& temperature_factor(size_t const index) const
    {
        return m_tempfactor[index];
    }

    float& mass(size_t const index)
    {
        return m_mass[index];
    }

    float const& mass(size_t const index) const
    {
        return m_mass[index];
    }

    float& charge(size_t const index)
    {
        return m_charge[index];
    }

    float const& charge(size_t const index) const
    {
        return m_charge[index];
    }

    float& radius(size_t const index)
    {
        return m_radius[index];
    }

    float const& radius(size_t const index) const
    {
        return m_radius[index];
    }

    std::string& name(size_t const index)
    {
        return m_name[index];
    }

    std::string const& name(size_t const index) const
    {
        return m_name[index];
    }

    std::string& type(size_t const index)
    {
        return m_type[index];
    }

    std::string const& type(size_t const index) const
    {
        return m_type[index];
    }

    std::string& alternate_location(size_t const index)
    {
        return m_altloc[index];
    }

    std::string const& alternate_location(size_t const index) const
    {
        return m_altloc[index];
    }

    std::string& insertion_code(size_t const index)
    {
        return m_insertion_code[index];
    }

    std::string const& insertion_code(size_t const index) const
    {
        return m_insertion_code[index];
    }

private:
    size_t m_num_atoms;
    std::unordered_map<index_t, index_t> m_residue;
    std::vector<int> m_atomic;
    std::vector<float> m_occupancy;
    std::vector<float> m_tempfactor;
    std::vector<float> m_mass;
    std::vector<float> m_charge;
    std::vector<float> m_radius;
    std::vector<std::string> m_name;
    std::vector<std::string> m_type;
    std::vector<std::string> m_altloc;
    std::vector<std::string> m_insertion_code;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_ATOMDATA_HPP
