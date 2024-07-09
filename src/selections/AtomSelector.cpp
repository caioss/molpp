#include <molpp/AtomSelector.hpp>
#include <molpp/internal/MolData.hpp>
#include "selections/SelectionNode.hpp"
#include "selections/SelectionStack.hpp"
#include "selections/SelectionIndices.hpp"
#include "selections/SelectionParser.hpp"

namespace mol
{

AtomSelector::AtomSelector(std::string const& selection, internal::MolData& data)
: m_data{data}
{
    parse(selection);
}

AtomSel AtomSelector::apply(Frame frame)
{
    internal::SelectionIndices flags;
    size_t const num_atoms = m_data.size<Atom>();
    flags.available->reserve(num_atoms);
    for (index_t atom_idx = 0; atom_idx < num_atoms; atom_idx++)
    {
        flags.available->insert(atom_idx);
    }

    internal::SelectionStack stack;
    stack.push_node(m_tree);
    stack.push_indices(flags);
    while (!stack.empty_nodes())
    {
        std::shared_ptr<internal::SelectionNode> node = stack.pop_node();
        node->evaluate(stack, m_data, frame);
    }

    AtomSel sel(*(flags.selected), m_data);
    sel.set_frame(frame);

    return sel;
}

void AtomSelector::parse(std::string const& selection)
{
    m_tree = internal::default_parser().parse(selection);
}

} // namespace mol
