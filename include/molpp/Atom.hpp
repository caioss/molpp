#ifndef ATOM_HPP
#define ATOM_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/AtomAggregate.hpp>
#include <memory>
#include <vector>
#include <optional>

namespace mol {

class Bond;
class Residue;

class Atom : public internal::AtomAggregate<Atom>
{
public:
    Atom() = delete;
    using internal::AtomAggregate<Atom>::AtomAggregate;

    int resid() const;
    Residue residue();
    index_t residue_id() const;

    int atomic() const;
    void set_atomic(int const &atomic);

    float occupancy() const;
    void set_occupancy(float const &occupancy);

    float tempfactor() const;
    void set_tempfactor(float const &tempfactor);

    float mass() const;
    void set_mass(float const &mass);

    float charge() const;
    void set_charge(float const &charge);

    float radius() const;
    void set_radius(float const &radius);

    std::string name() const;
    void set_name(std::string const &name);

    std::string type() const;
    void set_type(std::string const &type);

    std::string resname() const;
    std::string segid() const;
    std::string chain() const;

    std::string altloc() const;
    void set_altloc(std::string const &altloc);

    std::shared_ptr<Bond> add_bond(index_t const bonded_to);
    std::shared_ptr<Bond> add_bond(Atom const &bonded_to);
    std::shared_ptr<Bond> bond(index_t const other);
    std::shared_ptr<Bond> bond(Atom const &other);

    std::vector<index_t> atom_indices() const
    {
        return {index()};
    }
};

} // namespace mol

#endif // ATOM_HPP
