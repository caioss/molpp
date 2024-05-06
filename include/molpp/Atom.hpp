#ifndef MOLPP_ATOM_HPP
#define MOLPP_ATOM_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/MolecularEntity.hpp>
#include <memory>
#include <vector>
#include <ranges>

namespace mol
{

class Bond;
class Residue;

class Atom : public internal::MolecularEntity
{
public:
    using internal::MolecularEntity::MolecularEntity;

    int resid() const;
    Residue residue();
    index_t residue_id() const;

    int atomic_number() const;
    void set_atomic_number(int const atomic);

    float occupancy() const;
    void set_occupancy(float const occupancy);

    float temperature_factor() const;
    void set_temperature_factor(float const tempfactor);

    float mass() const;
    void set_mass(float const mass);

    float charge() const;
    void set_charge(float const charge);

    float radius() const;
    void set_radius(float const radius);

    std::string const& name() const;
    void set_name(std::string const& name);

    std::string const& type() const;
    void set_type(std::string const& type);

    std::string const& residue_name() const;
    std::string const& segid() const;
    std::string const& chain() const;

    std::string const& alternate_location() const;
    void set_alternate_location(std::string const& altloc);

    std::string const& insertion_code() const;
    void set_insertion_code(std::string const& insertion_code);

    std::shared_ptr<Bond> add_bond(index_t const bonded_to);
    std::shared_ptr<Bond> add_bond(Atom const& bonded_to);
    std::shared_ptr<Bond> bond(index_t const other);
    std::shared_ptr<Bond> bond(Atom const& other);
    std::vector<std::shared_ptr<Bond>> bonds();

    Coord3::ColXpr position();
    Coord3::ConstColXpr position() const;

    auto as_atom_indices() const
    {
        return std::ranges::single_view(index());
    }
};

} // namespace mol

#endif // MOLPP_ATOM_HPP
