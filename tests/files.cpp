#include "files.hpp"
#include <string>

std::shared_ptr<MolData> load_topology(std::shared_ptr<MolReader> reader, std::string const file_name)
{
    return reader->read_topology(file_name);
}

std::shared_ptr<MolData> load_trajectory(std::shared_ptr<MolReader> reader, std::string const file_name, std::shared_ptr<MolData> atom_data)
{
    reader->read_trajectory(file_name, atom_data);
    return atom_data;
}

/*
 * PDB
 */
std::shared_ptr<MolReader> PDBFiles::reader = MolReader::from_file_ext(".pdb");
const std::shared_ptr<MolData> PDBFiles::tiny = load_topology(PDBFiles::reader, "tiny.pdb");
const std::shared_ptr<MolData> PDBFiles::big = load_topology(PDBFiles::reader, "4lad.pdb");
const std::shared_ptr<MolData> PDBFiles::traj = load_trajectory(PDBFiles::reader, "traj.pdb", load_topology(PDBFiles::reader, "traj.pdb"));

/*
 * Mol2
 */
std::shared_ptr<MolReader> Mol2Files::reader = MolReader::from_file_ext(".mol2");
const std::shared_ptr<MolData> Mol2Files::flben = load_trajectory(Mol2Files::reader, "fluorobenzene.mol2", load_topology(Mol2Files::reader, "fluorobenzene.mol2"));
