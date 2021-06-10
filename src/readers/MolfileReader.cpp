#include "MolfileReader.hpp"
#include "core/AtomData.hpp"
#include "Atom.hpp"
#include "MolError.hpp"
#include "molfile.h"
#include <regex>
#include <string>
#include <filesystem>
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
    std::unordered_map<std::string, molfile_plugin_t *> plugins;

private:
    MolfilePlugins()
    {
        if (pdbplugin_init() == VMDPLUGIN_SUCCESS) pdbplugin_register(this, register_cb);
    }

    static int register_cb(void *c_pointer, vmdplugin_t *p)
    {
        if (p == nullptr)
        {
            return VMDPLUGIN_ERROR;
        }
        MolfilePlugins *self = static_cast<MolfilePlugins *>(c_pointer);
        molfile_plugin_t *plugin = (molfile_plugin_t *)p;
        self->plugins[plugin->name] = plugin;
        return VMDPLUGIN_SUCCESS;
    }
};

MolfileReader::MolfileReader(std::string const &name)
: m_num_atoms { 0 },
  m_handle { nullptr },
  m_name(name)
{
    // Find correspondent plugin
    MolfilePlugins &plugins = MolfilePlugins::getInstance();
    if (plugins.plugins.count(name) == 0)
    {
        throw mol::MolError("Unknown plugin: " + name);
    }
    m_plugin = plugins.plugins[name];

    // Populate (real) extensions
    std::regex regexz(",");
    std::string const extensions(m_plugin->filename_extension);
    std::sregex_token_iterator ext_iter(extensions.begin(), extensions.end(),
                                        regexz, -1);
    std::sregex_token_iterator end;
    while (ext_iter != end)
    {
        m_extensions.insert("." + std::string(*ext_iter++));
    }
}

MolfileReader::~MolfileReader()
{
    close();
}

bool MolfileReader::can_read(const std::string &file_name) const
{
    std::string file_ext(std::filesystem::path(file_name).extension());
    return m_extensions.count(file_ext);
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

bool MolfileReader::open(const std::string &file_name)
{
    m_handle = m_plugin->open_file_read(file_name.c_str(),
                                        m_name.c_str(),
                                        &m_num_atoms);
    if (!m_handle || m_num_atoms <= 0)
    {
        m_num_atoms = 0;
        return false;
    }

    return true;
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

std::shared_ptr<AtomData> MolfileReader::read_atoms()
{
    std::vector<molfile_atom_t> molfile_atoms(m_num_atoms);
    int flags = MOLFILE_BADOPTIONS;
    if (m_plugin->read_structure(m_handle, &flags, molfile_atoms.data())
        || flags == (int)MOLFILE_BADOPTIONS)
    {
        return nullptr;
    }
    std::shared_ptr<AtomData> atom_data = AtomData::create(m_num_atoms);

    for (size_t i = 0; i < (size_t)m_num_atoms; ++i)
    {
        molfile_atom_t const &mol_atom = molfile_atoms[i];

        Atom atom = atom_data->index(i, 0);
        atom.set_name(mol_atom.name);
        atom.set_type(mol_atom.type);
        atom.set_resname(mol_atom.resname);
        atom.set_resid(mol_atom.resid);
        atom.set_segid(mol_atom.segid);
        atom.set_chain(mol_atom.chain);

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

    }

    if (has_bonds())
    {
        // TODO
    }

    return atom_data;
}

bool MolfileReader::skip_timestep(std::shared_ptr<AtomData> atom_data)
{
    return m_plugin->read_next_timestep(m_handle, atom_data->size(), nullptr) == MOLFILE_SUCCESS;
}

bool MolfileReader::read_timestep(std::shared_ptr<AtomData> atom_data)
{
    Timestep ts(atom_data->size());
    molfile_timestep_t mol_ts;
    mol_ts.coords = ts.coords().data();
    mol_ts.physical_time = 0.0;
    int rc = m_plugin->read_next_timestep(m_handle, atom_data->size(), &mol_ts);
    if (rc != MOLFILE_SUCCESS)
    {
        return false;
    }

    atom_data->add_timestep(std::move(ts));

    return true;
}
