#ifndef MOLPP_SELECTIONS_SELECTIONNODE_HPP
#define MOLPP_SELECTIONS_SELECTIONNODE_HPP

#include <molpp/Common.hpp>

#include <memory>

namespace mol::internal
{

class MolData;
class SelectionStack;

class SelectionNode
{
public:
    virtual ~SelectionNode()
    {}

    virtual void evaluate(SelectionStack& stack, MolData const& data, Frame frame) const = 0;

    std::shared_ptr<SelectionNode> left;
    std::shared_ptr<SelectionNode> right;
};

} // namespace mol::internal

#endif // MOLPP_SELECTIONS_SELECTIONNODE_HPP
