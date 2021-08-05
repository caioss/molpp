#include "Residue.hpp"
#include "AtomSel.hpp"
#include "core/MolData.hpp"

using namespace mol;

int Residue::resid() const
{
    return m_data->residues().resid(m_index);
}

void Residue::set_resid(int const &resid)
{
    m_data->residues().resid(m_index) = resid;
}

std::string Residue::resname() const
{
    return m_data->residues().resname(m_index);
}

void Residue::set_resname(std::string const &resname)
{
    m_data->residues().resname(m_index) = resname;
}

std::string Residue::segid() const
{
    return m_data->residues().segid(m_index);
}

void Residue::set_segid(std::string const &segid)
{
    m_data->residues().segid(m_index) = segid;
}

std::string Residue::chain() const
{
    return m_data->residues().chain(m_index);
}

void Residue::set_chain(std::string const &chain)
{
    m_data->residues().chain(m_index) = chain;
}

void Residue::add_atom(size_t index)
{
    size_t const old_res = m_data->properties().residue(index);
    m_data->residues().remove_atom(old_res, index);
    m_data->residues().add_atom(m_index, index);
    m_data->properties().residue(index) = m_index;
}

void Residue::add_atom(Atom &atom)
{
    add_atom(atom.index());
}

std::shared_ptr<AtomSel> Residue::atoms()
{
    auto indices = m_data->residues().indices(m_index);
    return std::make_shared<AtomSel>(m_data->size(), std::vector<size_t>(indices.begin(), indices.end()), m_data);
}

size_t Residue::size() const
{
    return m_data->residues().size(m_index);
}
