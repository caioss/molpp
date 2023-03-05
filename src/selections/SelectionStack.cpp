#include "selections/SelectionStack.hpp"
#include "selections/SelectionNode.hpp"

using namespace mol;
using namespace mol::internal;

SelectionStack::SelectionStack(std::shared_ptr<SelectionNode> root)
: m_root(root)
{
}

void SelectionStack::push_flags(SelectionFlags const& flags)
{
    m_flags.push_front(flags);
}

void SelectionStack::push(std::shared_ptr<SelectionNode> node, SelectionFlags const& flags)
{
    m_callers.push_front(node);
    push_flags(flags);
}

SelectionFlags SelectionStack::pop_flags()
{
    SelectionFlags flags = m_flags.front();
    m_flags.pop_front();
    return flags;
}

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
