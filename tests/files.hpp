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
    PDBFiles(PDBFiles const &) = delete;
    void operator=(PDBFiles const &) = delete;

    static std::shared_ptr<MolReader> const reader();
    static MolData* tiny();
    static MolData* big();
    static MolData* traj();

private:
    PDBFiles();
    static PDBFiles const& instance();

    std::shared_ptr<MolReader> m_reader;
    std::unique_ptr<MolData> m_tiny;
    std::unique_ptr<MolData> m_big;
    std::unique_ptr<MolData> m_traj;
};

class Mol2Files
{
public:
    Mol2Files(Mol2Files const &) = delete;
    void operator=(Mol2Files const &) = delete;

    static std::shared_ptr<MolReader> const reader();
    static std::shared_ptr<MolData> const flben();

private:
    Mol2Files();
    static Mol2Files const& instance();

    std::shared_ptr<MolReader> m_reader;
    std::shared_ptr<MolData> m_flben;
};

#endif // FILES_HPP
