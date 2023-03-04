#include "selections/selections.hpp"
#include "core/MolData.hpp"
#include <molpp/MolError.hpp>

using namespace mol;
using namespace mol::internal;

void SelectionStack::evaluate(MolData const& data, SelectionFlags& flags, Frame frame)
{
    // We use an explicit stack to allow arbitrarily sized selections
    m_callers.clear();
    m_flags.clear();
    push(m_root, flags);

    while (!m_callers.empty())
    {
        std::shared_ptr<SelectionNode> node = m_callers.front();
        m_callers.pop_front();
        node->evaluate(*this, data, frame);
    }
}
