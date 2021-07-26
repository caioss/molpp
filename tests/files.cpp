#include "files.hpp"
#include <string>

std::shared_ptr<MolData> load_mol_file(std::shared_ptr<MolReader> reader, std::string const file_name)
{
    auto atom_data = reader->read_topology(file_name);
    reader->read_trajectory(file_name, atom_data);
    return atom_data;
}

/*
 * PDB
 */
std::shared_ptr<MolReader> PDBFiles::reader = MolReader::from_file_ext(".pdb");
const std::shared_ptr<MolData> PDBFiles::tiny = load_mol_file(PDBFiles::reader, "tiny.pdb");
const std::shared_ptr<MolData> PDBFiles::traj = load_mol_file(PDBFiles::reader, "traj.pdb");

/*
 * Mol2
 */
std::shared_ptr<MolReader> Mol2Files::reader = MolReader::from_file_ext(".mol2");
const std::shared_ptr<MolData> Mol2Files::flben = load_mol_file(Mol2Files::reader, "fluorobenzene.mol2");
