#ifndef FILES_HPP
#define FILES_HPP

#include "core/AtomData.hpp"
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
        ASSERT_THAT(traj, NotNull());
    }

    static std::shared_ptr<MolReader> reader;
    static const std::shared_ptr<AtomData> tiny;
    static const std::shared_ptr<AtomData> traj;
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
    static const std::shared_ptr<AtomData> flben;
};

#endif // FILES_HPP
