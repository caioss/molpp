#include "files.hpp"
#include <string>

std::unique_ptr<MolData> load_topology(std::shared_ptr<MolReader> reader, std::string const file_name)
{
    return reader->read_topology(file_name);
}

std::unique_ptr<MolData> load_trajectory(std::shared_ptr<MolReader> reader, std::string const file_name, std::unique_ptr<MolData> atom_data)
{
    reader->read_trajectory(file_name, *atom_data);
    return atom_data;
}

/*
 * PDB
 */
std::shared_ptr<MolReader> const PDBFiles::reader()
{
    return instance().m_reader;
}

MolData* PDBFiles::tiny()
{
    return instance().m_tiny.get();
}

MolData* PDBFiles::big()
{
    return instance().m_big.get();
}

MolData* PDBFiles::traj()
{
    return instance().m_traj.get();
}

PDBFiles::PDBFiles()
: m_reader{MolReader::from_file_ext(".pdb")},
  m_tiny{m_reader->read_topology("tiny.pdb")},
  m_big{m_reader->read_topology("4lad.pdb")},
  m_traj{m_reader->read_topology("traj.pdb")}
{
    m_reader->read_trajectory("4lad.pdb", *m_big);
    m_reader->read_trajectory("traj.pdb", *m_traj);
}

PDBFiles const& PDBFiles::instance()
{
    static PDBFiles instance;
    return instance;
}

/*
 * Mol2
 */
std::shared_ptr<MolReader> const Mol2Files::reader()
{
    return instance().m_reader;
}

std::shared_ptr<MolData> const Mol2Files::flben()
{
    return instance().m_flben;
}

Mol2Files::Mol2Files()
: m_reader{MolReader::from_file_ext(".mol2")},
  m_flben{m_reader->read_topology("fluorobenzene.mol2")}
{
    m_reader->read_trajectory("fluorobenzene.mol2", *m_flben);
}

Mol2Files const& Mol2Files::instance()
{
    static Mol2Files instance;
    return instance;
}
