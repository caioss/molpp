#include "selections/boolean.hpp"
#include "selections/SelectionStack.hpp"

using namespace mol;
using namespace mol::internal;

void OrSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const
{
    SelectionFlags flags = stack.pop_flags();
    stack.push(right, flags);
    stack.push(left, flags);
}

void AndSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const
{
    SelectionFlags flags = stack.pop_flags();
    std::shared_ptr<std::set<index_t>> parcial = std::make_shared<std::set<index_t>>();
    // Short-circuit
    stack.push(right, {parcial, flags.selected});
    stack.push(left, {flags.mask, parcial});
}

class NotImpl : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const override
    {
        SelectionFlags inverted = stack.pop_flags();
        SelectionFlags result = stack.pop_flags();
        std::set_difference(inverted.mask->begin(), inverted.mask->end(),
                            inverted.selected->begin(), inverted.selected->end(),
                            std::inserter(*(result.selected), result.selected->begin()));
    }
};

void NotSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const
{
    SelectionFlags flags = stack.pop_flags();
    SelectionFlags inverted(flags.mask, std::make_shared<std::set<index_t>>());
    // Process selection separetely and then combine inside NotImpl
    stack.push(std::make_shared<NotImpl>(), flags);
    stack.push_flags(inverted);
    stack.push(left, inverted);
}
