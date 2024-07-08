#ifndef MOLPP_SELECTIONS_PROPERTYSELECTION_HPP
#define MOLPP_SELECTIONS_PROPERTYSELECTION_HPP

#include "selections/SelectionNode.hpp"
#include "selections/SelectionStack.hpp"
#include "selections/SelectionIndices.hpp"

namespace mol::internal
{

template<class Derived>
class PropertySelection : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& data, Frame /*frame*/) const override
    {
        Derived const& derived = *static_cast<Derived const*>(this);
        SelectionIndices indices = stack.pop_indices();

        for (index_t atom_index : *(indices.available))
        {
            if (derived.evaluate_atom(atom_index, data))
            {
                indices.selected->insert(atom_index);
            }
        }
    }
};

} // namespace mol::internal

#endif // MOLPP_SELECTIONS_PROPERTYSELECTION_HPP
