#include "SelectionStack.hpp"

namespace mol::internal
{

void mol::internal::SelectionStack::clear()
{
    m_nodes.clear();
    m_indices.clear();
}

bool mol::internal::SelectionStack::empty_nodes() const
{
    return m_nodes.empty();
}

bool mol::internal::SelectionStack::empty_indices() const
{
    return m_indices.empty();
}

void SelectionStack::push_indices(SelectionIndices const& indices)
{
    m_indices.push_back(indices);
}

void SelectionStack::push_node(std::shared_ptr<SelectionNode> node)
{
    m_nodes.push_back(node);
}

SelectionIndices SelectionStack::pop_indices()
{
    SelectionIndices indices = m_indices.back();
    m_indices.pop_back();
    return indices;
}

std::shared_ptr<SelectionNode> mol::internal::SelectionStack::pop_node()
{
    std::shared_ptr<SelectionNode> node = m_nodes.back();
    m_nodes.pop_back();
    return node;
}

} // namespace mol::internal
