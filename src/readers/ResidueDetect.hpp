#ifndef MOLPP_READERS_RESIDUEDETECT_HPP
#define MOLPP_READERS_RESIDUEDETECT_HPP

#include <molpp/MolppCore.hpp>

#include <map>
#include <tuple>
#include <string>

namespace mol
{

namespace internal
{

class MolData;

class ResidueDetect
{
public:
    index_t register_atom(int const resid, std::string const& resname, std::string const& segid, std::string const& chain);
    void update_residue_data(MolData& mol_data) const;

private:
    struct ResidueKey
    {
        int resid;
        index_t chain;
        index_t segment;
        std::string resname;
    };

    struct EntityInfo
    {
        index_t index;
        size_t size;
    };

    friend bool operator<(ResidueKey const& lhs, ResidueKey const& rhs);

    std::map<ResidueKey, EntityInfo> m_residues;
    std::map<std::string, EntityInfo> m_chains;
    std::map<std::string, EntityInfo> m_segments;
    std::vector<std::string> m_chain_name;
    std::vector<std::string> m_segment_name;
};

} // namespace internal
} // namespace mol

#endif // MOLPP_READERS_RESIDUEDETECT_HPP
