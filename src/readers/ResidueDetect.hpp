#ifndef RESIDUEDETECT_HPP
#define RESIDUEDETECT_HPP

#include "core/MolData.hpp"
#include <map>
#include <string>

namespace mol::internal {

class ResidueDetect
{
private:
    struct Residue
    {
        size_t index;
        size_t count;
        int resid;
        std::string resname;
        std::string segid;
        std::string chain;
    };

    using residues_type = std::map<std::tuple<int, std::string, std::string, std::string>, Residue>;

    residues_type m_residues;
    residues_type::iterator m_iterator;
    Residue m_current;

public:
    ResidueDetect()
    : m_iterator{m_residues.end()}
    {}

    size_t register_atom(int const resid, std::string const &resname, std::string const &segid, std::string const chain)
    {
        if (m_iterator == m_residues.end()
            || m_current.resid != resid
            || m_current.resname != resname
            || m_current.segid != segid
            || m_current.chain != chain)
        {
            size_t const index = m_residues.size();
            m_current = {index, 0, resid, resname, segid, chain};

            std::tuple const key{resid, resname, segid, chain};
            m_iterator = m_residues.insert(std::pair(key, m_current)).first;
        }

        Residue &residue = (*m_iterator).second;
        residue.count++;
        return residue.index;
    }

    void update_residue_data(MolData &mol_data) const
    {
        ResidueData &residues_data = mol_data.residues();
        residues_data.resize(m_residues.size());
        for (auto const &item : m_residues)
        {
            Residue const &residue = item.second;
            residues_data.reset(residue.index, residue.count);
            residues_data.set(residue.index, residue.resid, residue.resname, residue.segid, residue.chain);
        }

        for (size_t index = 0; index < mol_data.size(); ++index)
        {
            size_t const residue_idx = mol_data.properties().residue(index);
            residues_data.add_atom(residue_idx, index);
        }
    }
};

} // namespace mol::internal

#endif // RESIDUEDETECT_HPP
