#include "MolReader.hpp"
#include "MolfileReader.hpp"

using namespace mol::internal;

std::shared_ptr<MolReader> MolReader::from_file_ext(const std::string &file_ext)
{
    if (MolfileReader::can_read(file_ext))
    {
        return std::make_shared<MolfileReader>(file_ext);
    }

    return nullptr;
}

std::shared_ptr<AtomData> MolReader::read_topology(std::string const &file_name)
{
    if (!(has_topology() && open(file_name)))
    {
        return nullptr;
    }

    std::shared_ptr<AtomData> atom_data = read_atoms();
    close();
    return atom_data;
}

bool MolReader::read_trajectory(std::string const &file_name, std::shared_ptr<AtomData> atom_data, int begin, int end, int step)
{
    if (!(has_trajectory() && open(file_name)))
    {
        return false;
    }

    begin = (begin < 0) ? 0 : begin;
    step = (step < 1) ? 1 : step;
    int current;
    for (current = 0; current < begin; ++current)
    {
        if (!skip_timestep(atom_data))
        {
            close();
            return true;
        }
    }

    while (end < 0 || current++ < end)
    {
        if (!read_timestep(atom_data))
        {
            break;
        }

        for (int i = 1; i < step; ++i, ++current)
        {
            if (!skip_timestep(atom_data))
            {
                break;
            }
        }
    }

    close();
    return true;
}
