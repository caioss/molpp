/*
 * Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * This is a port of version 2.2.1 of DSSP to Molpp. DSSP is
 * distributed under the Boost Software License, Version 1.0.
 * Modifications here present were made by Caio S. Souza and are
 * distributed under the BSD 3-Clause License.
 */

#ifndef MOLPP_STRUCTURE_HPP
#define MOLPP_STRUCTURE_HPP

#include "analysis/dssp/DSSP.hpp"
#include <string>
#include <vector>

namespace dssp
{

enum class BridgeType
{
    None,
    Parallel,
    AntiParallel
};

enum class HelixType
{
    None,
    Start,
    End,
    StartAndEnd,
    Middle
};

struct HBond
{
    HBond()
    : residue{0}
    , energy{0.0}
    {}

    size_t residue;
    double energy;
};

class MResidue
{
public:
    MResidue(size_t const index, std::string chain_id, bool const isProline, mol::Point3 N, mol::Point3 CA, mol::Point3 C, mol::Point3 O);

    size_t chain() const;
    bool is_proline() const;
    bool is_chain_break() const;
    mol::Point3 const& N() const;
    mol::Point3 const& CA() const;
    mol::Point3 const& C() const;
    mol::Point3 const& O() const;
    mol::Point3 const& H() const;
    void set_previous(MResidue const& previous);
    HelixType helix_flag(uint32_t const stride) const;
    void set_helix_flag(uint32_t const stride, HelixType const type);
    bool is_helix_start(uint32_t const stride) const;
    bool is_valid_distance(MResidue const& inNext) const;

    bool is_bend;
    size_t index;
    uint32_t sheet;
    mol::SecondaryStructure structure;
    HBond h_bond_donor[2];
    HBond h_bond_acceptor[2];

private:
    bool m_is_proline;
    bool m_is_chain_break;
    size_t m_chain;
    mol::Point3 m_N;
    mol::Point3 m_CA;
    mol::Point3 m_C;
    mol::Point3 m_O;
    mol::Point3 m_H;
    HelixType m_helix_flags[3];
};

class MProtein
{
public:
    MProtein(size_t const num_residues);
    void compute_secondary_structure(bool prefer_pi_helices = true);

    template<class... Args>
    MResidue& emplace_residue(Args&&... args)
    {
        size_t const index = m_residues.size();
        MResidue& residue = m_residues.emplace_back(index, std::forward<Args>(args)...);
        if (index > 0)
        {
            residue.set_previous(m_residues[index - 1]);
        }
        return residue;
    }

private:
    void compute_h_bond_energies();
    void compute_helices(bool prefer_pi_helices);
    void compute_sheets();
    bool test_bond(size_t const first, size_t const second) const;
    BridgeType test_bridge(size_t const first, size_t const second) const;
    double compute_h_bond(MResidue& inDonor, MResidue& inAcceptor);
    bool no_chain_break(size_t const fisrt, size_t const last) const;
    double compute_kappa(size_t const index) const;
    bool is_last(size_t const index) const;

    std::vector<MResidue> m_residues;
};

} // namespace dssp

#endif // MOLPP_STRUCTURE_HPP
