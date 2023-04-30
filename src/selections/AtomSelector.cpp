#include <molpp/MolError.hpp>
#include <molpp/AtomSelector.hpp>
#include "core/MolData.hpp"
#include "selections/SelectionStack.hpp"
#include "selections/SelectionParser.hpp"
#include <set>
#include <stack>

using namespace mol;
using namespace mol::internal;

AtomSelector::AtomSelector(std::string const& selection, MolData* data)
: m_data(data)
{
    parse(selection);
}

AtomSel AtomSelector::apply(Frame frame)
{
    SelectionFlags flags;
    for (index_t atom_idx = 0; atom_idx < m_data->size(); atom_idx++)
    {
        flags.mask->insert(atom_idx);
    }

    SelectionStack sel_stack(m_tree);
    sel_stack.evaluate(*m_data, flags, frame);

    AtomSel sel(*(flags.selected), m_data);
    if (frame)
    {
        sel.set_frame(frame);
    }

    return sel;
}

void AtomSelector::parse(std::string const& selection)
{
    m_tree = SEL_PARSER.parse(selection);
}
