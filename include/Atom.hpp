#ifndef ATOM_HPP
#define ATOM_HPP

#include "MolppCore.hpp"
#include <memory>

namespace mol {

class AtomSel;
class Bond;

namespace internal {
    class AtomData;
}

class Atom
{
public:
    Atom(size_t const index, size_t const frame,std::shared_ptr<internal::AtomData> data);

    bool operator==(Atom const &other) const;

    size_t index() const;
    size_t frame() const;

    int resid() const;
    void set_resid(int const &resid);

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
    void set_resname(std::string const &resname);

    std::string segid() const;
    void set_segid(std::string const &segid);

    std::string chain() const;
    void set_chain(std::string const &chain);

    std::string altloc() const;
    void set_altloc(std::string const &altloc);

    std::shared_ptr<Bond> add_bond(size_t const bonded_to);
    std::shared_ptr<Bond> bond(size_t const other);
    std::vector<std::shared_ptr<Bond>> bonds();
    std::shared_ptr<AtomSel> bonded();

    Eigen::Ref<Pos3> coords();

private:
    size_t m_index;
    size_t m_frame;
    std::shared_ptr<internal::AtomData> m_data;
};

} // namespace mol

#endif // ATOM_HPP
