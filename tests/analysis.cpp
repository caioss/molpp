#include "analysis/dssp/DSSP.hpp"
#include "files.hpp"
#include <molpp/MolError.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

std::vector<mol::SecondaryStructure> const& structure_4lad();

TEST(Analysis, DSSP)
{
    MolSystem mol("4lad.pdb"); // TODO we should have a MolSystem already loaded somewhere
    mol.add_trajectory("4lad.pdb");
    DSSP dssp(mol);
    auto exp = structure_4lad();
    auto res = dssp.run(0);
    res = dssp.run(0);
    res = dssp.run(0);
    EXPECT_EQ(dssp.run(0), structure_4lad());

    // TODO run twice
}

std::vector<mol::SecondaryStructure> const& structure_4lad()
{
    static std::vector<mol::SecondaryStructure> structure{
        mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Bend, mol::Bend, mol::Turn, mol::Turn, mol::Loop, mol::Turn, mol::Turn, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Bend, mol::Turn, mol::Turn, mol::Turn, mol::Turn, mol::Loop, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Turn, mol::Turn, mol::Bend, mol::Bend, mol::Loop, mol::Loop, mol::Strand, mol::Strand, mol::Strand, mol::Strand, mol::Loop, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Bridge, mol::Loop, mol::Turn, mol::Turn, mol::Loop, mol::Bridge, mol::Bridge, mol::Loop, mol::Loop, mol::Helix3, mol::Helix3, mol::Helix3, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Turn, mol::Turn, mol::Loop, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Bend, mol::Loop, mol::Loop, mol::Bend, mol::Bend, mol::Loop, mol::Loop, mol::Bend, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Turn, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Turn, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Bend, mol::Bend, mol::Bend, mol::Loop, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Strand, mol::Strand, mol::Strand, mol::Loop, mol::Turn, mol::Turn, mol::Bend, mol::Loop, mol::Strand, mol::Strand, mol::Strand, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Turn, mol::Loop, mol::Bend, mol::Bend, mol::Loop, mol::Bend, mol::Bend, mol::Bend, mol::Loop, mol::Loop, mol::Loop, mol::Loop, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Helix, mol::Turn, mol::Loop, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown, mol::Unknown,
    };

    return structure;
}
