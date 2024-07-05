#include "selections/boolean.hpp"
#include "selections/SelectionStack.hpp"
#include "selections/SelectionIndices.hpp"

namespace mol::internal
{

void OrSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const
{
    // Just merge the results of the two operands
    SelectionIndices indices = stack.pop_indices();
    stack.push_node(right);
    stack.push_indices(indices);

    stack.push_node(left);
    stack.push_indices(indices);
}

void AndSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const
{
    SelectionIndices indices = stack.pop_indices();

    // Short-circuit is implemented by using the results of one operand as the available set for the other
    SelectionIndices::IndexSet merged = SelectionIndices::make_index_set();

    stack.push_node(right);
    stack.push_indices({merged, indices.selected});

    stack.push_node(left);
    stack.push_indices({indices.available, merged});
}

class InvertSelection : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const override
    {
        SelectionIndices source = stack.pop_indices();
        SelectionIndices result = stack.pop_indices();

        for (index_t const index : *result.available)
        {
            if (!source.selected->contains(index))
            {
                result.selected->insert(index);
            }
        }
    }
};

void NotSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const
{
    SelectionIndices indices = stack.pop_indices();

    // These will contain the inverse of the desired selection
    SelectionIndices inverse(indices.available, SelectionIndices::make_index_set());

    // Process the operand and then invert the results using InvertSelection
    stack.push_node(std::make_shared<InvertSelection>());
    stack.push_indices(indices);
    stack.push_indices(inverse);

    stack.push_node(left);
    stack.push_indices(inverse);
}

void AllSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, Frame /*frame*/) const
{
    SelectionIndices indices = stack.pop_indices();
    *indices.selected = *indices.available;
}

} // namespace mol::internal
