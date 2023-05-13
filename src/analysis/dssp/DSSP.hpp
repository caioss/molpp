#ifndef DSSP_HPP
#define DSSP_HPP

#include <molpp/Atom.hpp>
#include <molpp/MolSystem.hpp>
#include <vector>

namespace dssp
{
class MResidue;
} // namespace dssp

namespace mol
{

enum SecondaryStructure {
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
    SSResidue() = default;
    SSResidue(Residue& residue);
    void set_frame(Frame const frame);
    bool is_amino_acid() const;
    SecondaryStructure secondary_structure() const;
    bool is_proline() const;
    std::string const& chain() const;


// private: // TODO
    bool m_is_proline;
    Atom N;
    Atom CA;
    Atom C;
    Atom O;
    std::string m_chain;
    dssp::MResidue* m_residue;
};

class DSSP
{
public:
    DSSP(MolSystem& molecule);
    std::vector<SecondaryStructure> run(Frame frame);

private:
    std::vector<SSResidue> m_residues;
};

} // namespace mol

#endif // DSSP_HPP
