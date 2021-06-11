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
    if (!has_topology() || open(file_name) != SUCCESS)
    {
        return nullptr;
    }

    std::shared_ptr<AtomData> atom_data = read_atoms();
    close();
    return atom_data;
}

MolReader::Status MolReader::read_trajectory(std::string const &file_name, std::shared_ptr<AtomData> atom_data, int begin, int end, int step)
{
    // Sanity checks
    if (!has_trajectory())
    {
        return INVALID;
    }

    Status status = open(file_name);
    if (status != SUCCESS)
    {
        return status;
    }

    status = check_timestep_read(atom_data);
    if (status != SUCCESS)
    {
        close();
        return status;
    }

    // Initial skip
    begin = (begin < 0) ? 0 : begin;
    step = (step < 1) ? 1 : step;
    int current;
    for (current = 0; current < begin; ++current)
    {
        status = skip_timestep(atom_data);
        if (status != SUCCESS)
        {
            close();
            return (status == END) ? SUCCESS : status;
        }
    }

    // Read frames
    status = SUCCESS;
    while (end < 0 || current++ < end)
    {
        status = read_timestep(atom_data);
        if (status != SUCCESS)
        {
            break;
        }

        for (int i = 1; i < step; ++i, ++current)
        {
            status = skip_timestep(atom_data);
            if (status != SUCCESS)
            {
                break;
            }
        }
    }

    close();
    return (status == END) ? SUCCESS : status;
}
