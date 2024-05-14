#ifndef MOLPP_READERS_RESIDUEDETECT_HPP
#define MOLPP_READERS_RESIDUEDETECT_HPP

#include <molpp/MolppCore.hpp>

#include <map>
#include <string>

namespace mol
{

namespace internal
{

class MolData;

class ResidueDetect
{
public:
    ResidueDetect();
    index_t register_atom(int const resid, std::string const& resname, std::string const& segid, std::string const chain);
    void update_residue_data(MolData& mol_data) const;

private:
    struct ResidueInfo
    {
        index_t index;
        size_t count;
        int resid;
        std::string resname;
        std::string segid;
        std::string chain;
    };

    using residues_map = std::map<std::tuple<int, std::string, std::string, std::string>, ResidueInfo>;

    residues_map m_residues;
    residues_map::iterator m_iterator;
    ResidueInfo m_current;
};

} // namespace internal
} // namespace mol

#endif // MOLPP_READERS_RESIDUEDETECT_HPP
