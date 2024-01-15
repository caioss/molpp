#include "files.hpp"
#include "matchers.hpp"
#include "auxiliary.hpp"
#include <molpp/internal/AtomData.hpp>
#include <molpp/internal/MolData.hpp>
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/MolError.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

TEST(Atoms, MolData) {
    size_t const num_atoms { 3 };
    MolData data(num_atoms);
    EXPECT_EQ(data.size(), num_atoms);
    EXPECT_EQ(data.atoms().size(), num_atoms);
    EXPECT_EQ(data.bonds().size(), 0);
    EXPECT_EQ(data.residues().size(), 0);
}

TEST(Atoms, AtomData) {
    AtomData props(1);
    EXPECT_EQ(props.size(), 1);
    EXPECT_EQ(props.residue(0), -1);
    EXPECT_EQ(props.atomic(0), 0);
    EXPECT_EQ(props.occupancy(0), 0);
    EXPECT_EQ(props.tempfactor(0), 0);
    EXPECT_EQ(props.mass(0), 0);
    EXPECT_EQ(props.charge(0), 0);
    EXPECT_EQ(props.radius(0), 0);
    EXPECT_EQ(props.name(0), "");
    EXPECT_EQ(props.type(0), "");
    EXPECT_EQ(props.altloc(0), "");
}
