#ifndef MOLPP_SELECTIONS_SELECTIONSTACK_HPP
#define MOLPP_SELECTIONS_SELECTIONSTACK_HPP

#include "selections/SelectionIndices.hpp"

#include <memory>
#include <deque>

namespace mol::internal
{

class SelectionNode;

class SelectionStack
{
public:
    void clear();
    bool empty_nodes() const;
    bool empty_indices() const;
    void push_indices(SelectionIndices const& indices);
    void push_node(std::shared_ptr<SelectionNode> node);
    SelectionIndices pop_indices();
    std::shared_ptr<SelectionNode> pop_node();

private:
    std::deque<std::shared_ptr<SelectionNode>> m_nodes;
    std::deque<SelectionIndices> m_indices;
};

} // namespace mol::internal

#endif // MOLPP_SELECTIONS_SELECTIONSTACK_HPP
