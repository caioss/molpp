#include "selections/boolean.hpp"
#include "selections/SelectionStack.hpp"
#include "selections/SelectionIndices.hpp"

using namespace mol;
using namespace mol::internal;

void OrSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const
{
    SelectionIndices indices = stack.pop_indices();
    stack.push_node(right);
    stack.push_indices(indices);

    stack.push_node(left);
    stack.push_indices(indices);
}

void AndSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const
{
    SelectionIndices indices = stack.pop_indices();
    std::shared_ptr<std::unordered_set<index_t>> parcial = std::make_shared<std::unordered_set<index_t>>();

    // Short-circuit
    stack.push_node(right);
    stack.push_indices({parcial, indices.selected});

    stack.push_node(left);
    stack.push_indices({indices.available, parcial});
}

class NotImpl : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const override
    {
        SelectionIndices inverted = stack.pop_indices();
        SelectionIndices result = stack.pop_indices();

        std::copy_if(inverted.available->begin(), inverted.available->end(), std::inserter(*result.selected, result.selected->end()), [&](index_t const index) {
            return !inverted.selected->contains(index);
        });
    }
};

void NotSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const
{
    SelectionIndices indices = stack.pop_indices();
    SelectionIndices inverted(indices.available, std::make_shared<std::unordered_set<index_t>>());

    // Process selection separetely and then combine inside NotImpl
    stack.push_node(std::make_shared<NotImpl>());
    stack.push_indices(indices);
    stack.push_indices(inverted);

    stack.push_node(left);
    stack.push_indices(inverted);
}

void AllSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const
{
    SelectionIndices indices = stack.pop_indices();
    *indices.selected = *indices.available;
}
