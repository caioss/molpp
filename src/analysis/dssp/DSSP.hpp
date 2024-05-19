#ifndef MOLPP_ANALYSIS_DSSP_DSSP_HPP
#define MOLPP_ANALYSIS_DSSP_DSSP_HPP

#include <molpp/Atom.hpp>
#include <molpp/MolSystem.hpp>

#include <vector>
#include <optional>

namespace dssp
{
class Residue;
} // namespace dssp

namespace mol
{

enum SecondaryStructure
{
    Unknown,
    Loop,   //' '
    Helix,  // H
    Bridge, // B
    Strand, // E
    Helix3, // G
    Helix5, // I
    Turn,   // T
    Bend    // S
};

class SSResidue
{
public:
    SSResidue();
    SSResidue(Residue& residue);
    void set_frame(Frame const frame);
    bool is_amino_acid() const;
    SecondaryStructure secondary_structure() const;
    void set_secondary_structure(SecondaryStructure const structure);
    bool is_proline() const;
    bool is_chain_break() const;
    void set_chain_break(bool const is_break);
    std::string const& chain() const;
    std::optional<Atom> const& N() const;
    std::optional<Atom> const& CA() const;
    std::optional<Atom> const& C() const;
    std::optional<Atom> const& O() const;

private:
    bool m_is_proline;
    bool m_is_chain_break;
    SecondaryStructure m_structure;
    std::optional<Atom> m_N;
    std::optional<Atom> m_CA;
    std::optional<Atom> m_C;
    std::optional<Atom> m_O;
    std::string m_chain;
};

class DSSP
{
public:
    DSSP(MolSystem const& molecule);
    void run(Frame frame);
    std::vector<SSResidue> const& residues() const;
    std::vector<SecondaryStructure> secondary_structures() const;

private:
    size_t m_num_prot_aa;
    std::vector<SSResidue> m_residues;
};

} // namespace mol

#endif // MOLPP_ANALYSIS_DSSP_DSSP_HPP
