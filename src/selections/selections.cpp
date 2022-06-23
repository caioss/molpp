#include "selections/selections.hpp"
#include "core/MolData.hpp"
#include <molpp/MolError.hpp>

using namespace mol;
using namespace mol::internal;

void SelectionStack::evaluate(MolData const& data, SelectionFlags& flags, std::optional<size_t> frame)
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

void OrSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, std::optional<size_t> /*frame*/) const
{

    SelectionFlags flags = stack.pop_flags();
    stack.push(right, flags);
    stack.push(left, flags);
}

void AndSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, std::optional<size_t> /*frame*/) const
{
    SelectionFlags flags = stack.pop_flags();
    std::shared_ptr<std::set<size_t>> parcial = std::make_shared<std::set<size_t>>();
     // Short-circuit
    stack.push(right, {parcial, flags.selected});
    stack.push(left, {flags.mask, parcial});
}

class NotImpl : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& /*data*/, std::optional<size_t> /*frame*/) const override
    {
        SelectionFlags inverted = stack.pop_flags();
        SelectionFlags result = stack.pop_flags();
        std::set_difference(inverted.mask->begin(), inverted.mask->end(),
                            inverted.selected->begin(), inverted.selected->end(),
                            std::inserter(*(result.selected), result.selected->begin()));
    }
};

void NotSelection::evaluate(SelectionStack& stack, MolData const& /*data*/, std::optional<size_t> /*frame*/) const
{
    SelectionFlags flags = stack.pop_flags();
    SelectionFlags inverted(flags.mask, std::make_shared<std::set<size_t>>());
    // Process selection separetely and then combine inside NotImpl
    stack.push(std::make_shared<NotImpl>(), flags);
    stack.push_flags(inverted);
    stack.push(left, inverted);
}

void NumPropSelection::add_number(SelNumber const& number)
{
    m_numbers.push_back(number);
}

void NumPropSelection::add_range(SelNumberRange const& range)
{
    m_ranges.push_back(range);
}

void ResidSelection::evaluate(SelectionStack& stack, MolData const& data, std::optional<size_t> /*frame*/) const
{
    SelectionFlags flags = stack.pop_flags();
    AtomData const& atom_data = data.properties();
    ResidueData const& res_data = data.residues();

    for (size_t atom_idx : *(flags.mask))
    {
        size_t const res_idx = atom_data.residue(atom_idx);
        int const resid = res_data.resid(res_idx);
        if (has(resid))
        {
            flags.selected->insert(atom_idx);
        }
    }
}
