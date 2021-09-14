#ifndef FILES_HPP
#define FILES_HPP

#include "core/MolData.hpp"
#include "readers/MolReader.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace testing;
using namespace mol::internal;

class PDBFiles
{
public:
    void check()
    {
        ASSERT_THAT(reader, NotNull());
        ASSERT_THAT(tiny, NotNull());
        ASSERT_THAT(big, NotNull());
        ASSERT_THAT(traj, NotNull());
    }

    static std::shared_ptr<MolReader> reader;
    static const std::shared_ptr<MolData> tiny;
    static const std::shared_ptr<MolData> big;
    static const std::shared_ptr<MolData> traj;
};

class Mol2Files
{
public:
    void check()
    {
        ASSERT_THAT(reader, NotNull());
        ASSERT_THAT(flben, NotNull());
    }

    static std::shared_ptr<MolReader> reader;
    static const std::shared_ptr<MolData> flben;
};

#endif // FILES_HPP
