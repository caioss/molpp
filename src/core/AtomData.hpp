#ifndef ATOMDATA_HPP
#define ATOMDATA_HPP

#include "AtomProperties.hpp"
#include "Timestep.hpp"
#include "BondGraph.hpp"
#include <memory>
#include <vector>

namespace mol::internal {

class AtomData : public std::enable_shared_from_this<AtomData>
{
public:
    AtomData() = delete;
    static std::shared_ptr<AtomData> create(size_t const num_atoms);

    size_t size() const { return m_num_atoms; };
    AtomProperties &properties() { return m_properties; }
    BondGraph &bonds() { return m_bonds; }

    size_t num_frames() { return m_timestep.size(); }
    Timestep &timestep(size_t const index) { return m_timestep[index]; }
    void add_timestep(Timestep &&ts);

private:
    AtomData(size_t const num_atoms);

    size_t m_num_atoms;
    AtomProperties m_properties;
    std::vector<Timestep> m_timestep;
    BondGraph m_bonds;
};

} // namespace mol::internal

#endif // ATOMDATA_HPP
