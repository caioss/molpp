#include "files.hpp"
#include <string>

std::shared_ptr<AtomData> load_mol_file(std::shared_ptr<MolReader> reader, std::string const file_name)
{
    auto atom_data = reader->read_topology(file_name);
    reader->read_trajectory(file_name, atom_data);
    return atom_data;
}

std::shared_ptr<MolReader> PDBFiles::reader = MolReader::from_file_ext(".pdb");
const std::shared_ptr<AtomData> PDBFiles::tiny = load_mol_file(PDBFiles::reader, "tiny.pdb");
const std::shared_ptr<AtomData> PDBFiles::traj = load_mol_file(PDBFiles::reader, "traj.pdb");
