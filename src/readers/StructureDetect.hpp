#ifndef MOLPP_READERS_STRUCTUREDETECT_HPP
#define MOLPP_READERS_STRUCTUREDETECT_HPP

#include <molpp/MolppCore.hpp>

#include <map>
#include <string>

namespace mol
{

namespace internal
{

class MolData;

class StructureDetect
{
public:
    StructureDetect(MolData& data);
    void register_atom(index_t const atom_index, int const resid, std::string const& resname, std::string const& segment, std::string const& chain);
    void update_residue_data(MolData& mol_data) const;

private:
    struct ResidueKey
    {
        int resid;
        index_t chain;
        index_t segment;
        std::string resname;
    };

    index_t register_chain(std::string const& chain);
    index_t register_segment(std::string const& segment);
    friend bool operator<(ResidueKey const& lhs, ResidueKey const& rhs);

    MolData& m_data;
    std::map<ResidueKey, index_t> m_residues;
    std::map<std::string, index_t> m_chains;
    std::map<std::string, index_t> m_segments;
};

} // namespace internal
} // namespace mol

#endif // MOLPP_READERS_STRUCTUREDETECT_HPP
