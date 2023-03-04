#ifndef BOOLEAN_HPP
#define BOOLEAN_HPP

#include "selections/selections.hpp"

namespace mol::internal
{

class OrSelection : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& data, Frame frame) const override;
};

class AndSelection : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& data, Frame frame) const override;
};

class NotSelection : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& /*data*/, Frame frame) const override;
};

} // namespace mol::internal

#endif // BOOLEAN_HPP
