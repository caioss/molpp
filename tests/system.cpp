#include "MolSystem.hpp"
#include "Atom.hpp"
#include "MolError.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace testing;

#include <iostream>
TEST(molsystem, BasicAssertions) {
    EXPECT_THROW(MolSystem(""), MolError);
    EXPECT_THROW(MolSystem("file.unk"), MolError);
    EXPECT_THROW(MolSystem("no_file.pdb"), MolError);
    MolSystem mol("traj.pdb");

    EXPECT_THROW(mol.add_trajectory("traj.unk"), MolError);
    EXPECT_THROW(mol.add_trajectory("tiny.pdb"), MolError);
    EXPECT_THROW(mol.add_trajectory("no_file.pdb"), MolError);
    EXPECT_NO_THROW(mol.add_trajectory("traj.pdb"));

    auto all_sel = mol.all();
    ASSERT_THAT(all_sel, NotNull());
    EXPECT_EQ(all_sel->size(), 2);

    auto index_sel = mol.select({0});
    ASSERT_THAT(index_sel, NotNull());
    EXPECT_EQ(index_sel->size(), 1);
}
