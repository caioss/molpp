#include "molfile.h"
#include "MolfileReader.hpp"
#include "core/MolData.hpp"
#include "ResidueDetect.hpp"
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/MolError.hpp>
#include <map>
#include <tuple>
#include <regex>
#include <string>
#include <unordered_map>

using namespace mol::internal;

/*
 * Plugins initialization
 */
class MolfilePlugins
{
public:
    MolfilePlugins(MolfilePlugins const &) = delete;
    void operator=(MolfilePlugins const &) = delete;
    static MolfilePlugins &getInstance()
    {
        static MolfilePlugins instance;
        return instance;
    }
    molfile_plugin_t *find_plugin(std::string file_ext)
    {
        auto it = m_extensions.find(file_ext);
        return (it == m_extensions.end()) ? nullptr : it->second;
    }
    bool has_extension(std::string file_ext)
    {
        return m_extensions.count(file_ext);
    }

private:
    MolfilePlugins()
    {
        // Register plugins
        if (pdbplugin_init() == VMDPLUGIN_SUCCESS) pdbplugin_register(this, register_cb);
        if (mol2plugin_init() == VMDPLUGIN_SUCCESS) mol2plugin_register(this, register_cb);
        if (psfplugin_init() == VMDPLUGIN_SUCCESS) psfplugin_register(this, register_cb);
        if (gromacsplugin_init() == VMDPLUGIN_SUCCESS) gromacsplugin_register(this, register_cb);

        // Register extensions
        std::regex regexz(",");
        std::sregex_token_iterator end;

        for (molfile_plugin_t *plugin : m_plugins)
        {
            std::string const plugin_ext(plugin->filename_extension);
            std::sregex_token_iterator ext_iter(plugin_ext.begin(), plugin_ext.end(), regexz, -1);
            while (ext_iter != end)
            {
                m_extensions.insert({"." + std::string(*ext_iter++), plugin});
            }
        }
    }

    static int register_cb(void *c_pointer, vmdplugin_t *p)
    {
        if (p == nullptr)
        {
            return VMDPLUGIN_ERROR;
        }
        MolfilePlugins *self = static_cast<MolfilePlugins *>(c_pointer);
        molfile_plugin_t *plugin = (molfile_plugin_t *)p;
        self->m_plugins.push_back(plugin);
        return VMDPLUGIN_SUCCESS;
    }

    std::vector<molfile_plugin_t *> m_plugins;
    std::unordered_multimap<std::string, molfile_plugin_t *> m_extensions;

};

MolfileReader::MolfileReader(std::string const &file_ext)
: m_num_atoms { 0 },
  m_handle { nullptr }
{
    // Find correspondent plugin
    m_plugin = MolfilePlugins::getInstance().find_plugin(file_ext);
    if (!m_plugin)
    {
        throw mol::MolError("Unknown plugin for extension " + file_ext);
    }
    m_name = m_plugin->name;
}

MolfileReader::~MolfileReader()
{
    close();
}

bool MolfileReader::can_read(const std::string &file_ext)
{
    return MolfilePlugins::getInstance().has_extension(file_ext);
}

bool MolfileReader::has_topology() const
{
    return m_plugin->read_structure != nullptr;
}

bool MolfileReader::has_trajectory_metadata() const
{
    return m_plugin->read_timestep_metadata != nullptr;
}

bool MolfileReader::has_trajectory() const
{
    return m_plugin->read_next_timestep != nullptr;
}

bool MolfileReader::has_bonds() const
{
    return m_plugin->read_bonds != nullptr;
}

MolReader::Status MolfileReader::open(const std::string &file_name)
{
    if (m_handle)
    {
        return MolReader::INVALID;
    }

    m_handle = m_plugin->open_file_read(file_name.c_str(), m_name.c_str(), &m_num_atoms);

    if (!m_handle || m_num_atoms <= 0)
    {
        close();
        return MolReader::FAILED;
    }

    return MolReader::SUCCESS;
}

void MolfileReader::close()
{
    if (m_handle)
    {
        m_plugin->close_file_read(m_handle);
    }
    m_handle = nullptr;
    m_num_atoms = 0;
}

std::unique_ptr<MolData> MolfileReader::read_atoms()
{
    if (!m_handle)
    {
        throw mol::MolError("No opened file");
    }

    // Read data and allocate atoms
    std::vector<molfile_atom_t> molfile_atoms(m_num_atoms);
    int flags = MOLFILE_BADOPTIONS;
    if (m_plugin->read_structure(m_handle, &flags, molfile_atoms.data()) != MOLFILE_SUCCESS || flags == (int)MOLFILE_BADOPTIONS)
    {
        return nullptr;
    }

    std::unique_ptr<MolData> mol_data = std::make_unique<MolData>(m_num_atoms);
    ResidueDetect residues;

    // Loop over all atoms
    for (index_t i = 0; i < (size_t)m_num_atoms; ++i)
    {
        /*
         * Atoms properties
         */
        molfile_atom_t const &mol_atom = molfile_atoms[i];

        Atom atom(i, {}, mol_data.get());
        atom.set_name(mol_atom.name);
        atom.set_type(mol_atom.type);

        if (flags & MOLFILE_INSERTION)
        {
            atom.set_altloc(mol_atom.altloc);
        }

        if (flags & MOLFILE_OCCUPANCY)
        {
            atom.set_occupancy(mol_atom.occupancy);
        }

        if (flags & MOLFILE_BFACTOR)
        {
            atom.set_tempfactor(mol_atom.bfactor);
        }

        if (flags & MOLFILE_MASS)
        {
            atom.set_mass(mol_atom.mass);
        }

        if (flags & MOLFILE_CHARGE)
        {
            atom.set_charge(mol_atom.charge);
        }

        if (flags & MOLFILE_RADIUS)
        {
            atom.set_radius(mol_atom.radius);
        }

        if (flags & MOLFILE_ATOMICNUMBER)
        {
            atom.set_atomic(mol_atom.atomicnumber);
        }

        // Detect residues
        mol_data->atoms().residue(i) = residues.register_atom(mol_atom.resid, mol_atom.resname, mol_atom.segid, mol_atom.chain);
    }

    // Update residues data
    residues.update_residue_data(*mol_data);

    /*
     * Bonds
     */
    if (has_bonds())
    {
        int num_bonds = 0, num_types = 0;
        int *from = nullptr, *to = nullptr, *type = nullptr;
        float *order = nullptr;
        char **type_name = nullptr;

        int rc = m_plugin->read_bonds(m_handle, &num_bonds, &from, &to, &order,
                                      &type, &num_types, &type_name);
        if (rc == MOLFILE_SUCCESS && num_bonds > 0)
        {
            BondData &bond_graph = mol_data->bonds();
            bond_graph.set_incomplete(flags & MOLFILE_BONDSSPECIAL);

            for (index_t i = 0; i < (size_t)num_bonds; ++i)
            {
                auto bond = bond_graph.add_bond(from[i] - 1, to[i] - 1);
                if (!bond)
                {
                    continue;
                }

                bond->set_guessed(false);
                if (order)
                {
                    float const atom_order = order[i];
                    if (atom_order > 1 && atom_order < 2)
                    {
                        bond->set_order(0);
                        bond->set_aromatic(true);
                    }
                    else
                    {
                        bond->set_order(atom_order);
                        bond->set_guessed_order(false);
                    }

                }
            }
        }
    }

    return mol_data;
}

MolReader::Status MolfileReader::check_timestep_read(MolData& mol_data)
{
    if (!m_handle)
    {
        return INVALID;
    }

    if (mol_data.size() != (size_t)m_num_atoms)
    {
        return WRONG_ATOMS;
    }

    return SUCCESS;
}

MolReader::Status MolfileReader::skip_timestep(MolData& mol_data)
{
    switch (m_plugin->read_next_timestep(m_handle, mol_data.size(), nullptr))
    {
        case MOLFILE_SUCCESS:
            return SUCCESS;

        case MOLFILE_EOF:
            return END;

        default:
            return FAILED;
    }
}

MolReader::Status MolfileReader::read_timestep(MolData& mol_data)
{
    Timestep ts(mol_data.size());
    molfile_timestep_t mol_ts;
    mol_ts.coords = ts.coords().data();
    mol_ts.physical_time = 0.0;

    switch (m_plugin->read_next_timestep(m_handle, mol_data.size(), &mol_ts))
    {
        case MOLFILE_SUCCESS:
            break;

        case MOLFILE_EOF:
            return END;

        default:
            return FAILED;
    }

    mol_data.trajectory().add_timestep(std::move(ts));

    return SUCCESS;
}
