#include <molpp/AtomSelector.hpp>
#include <molpp/internal/MolData.hpp>
#include "selections/SelectionStack.hpp"
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
    internal::SelectionFlags flags;
    for (index_t atom_idx = 0; atom_idx < m_data.size<Atom>(); atom_idx++)
    {
        flags.mask->insert(atom_idx);
    }

    internal::SelectionStack sel_stack(m_tree);
    sel_stack.evaluate(m_data, flags, frame);

    AtomSel sel(*(flags.selected), m_data);
    sel.set_frame(frame);

    return sel;
}

void AtomSelector::parse(std::string const& selection)
{
    m_tree = internal::SEL_PARSER.parse(selection);
}

} // namespace mol
