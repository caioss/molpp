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
    ~SSResidue();
    SSResidue(SSResidue&&) = default;
    SSResidue(Residue& residue);
    void update(dssp::MResidue* previous);
    void set_frame(Frame const frame);
    void set_structure(SecondaryStructure const structure);
    bool is_amino_acid() const;
    dssp::MResidue& residue();
    SecondaryStructure secondary_structure() const;

private:
    bool m_is_proline;
    Atom N;
    Atom CA;
    Atom C;
    Atom O;
    std::string m_chain;
    std::unique_ptr<dssp::MResidue> m_residue; // TODO try to not use pointers
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
