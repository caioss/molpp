#ifndef SELECTIONSTACK_HPP
#define SELECTIONSTACK_HPP

#include <molpp/MolppCore.hpp>
#include <set>
#include <memory>
#include <forward_list>

namespace mol::internal
{

class MolData;
class SelectionNode;

struct SelectionFlags
{
    SelectionFlags(std::shared_ptr<std::set<index_t>> prev_mask, std::shared_ptr<std::set<index_t>> prev_selected)
    : mask{prev_mask}
    , selected{prev_selected}
    {
    }

    SelectionFlags()
    : mask{std::make_shared<std::set<index_t>>()}
    , selected{std::make_shared<std::set<index_t>>()}
    {
    }

    // Flags must be ordered
    std::shared_ptr<std::set<index_t>> mask;
    std::shared_ptr<std::set<index_t>> selected;
};

class SelectionStack
{
public:
    SelectionStack(std::shared_ptr<SelectionNode> root);
    void push_flags(SelectionFlags const& flags);
    void push(std::shared_ptr<SelectionNode> node, SelectionFlags const& flags);
    SelectionFlags pop_flags();
    void evaluate(MolData const& data, SelectionFlags& flags, Frame frame);

private:
    std::shared_ptr<SelectionNode> m_root;
    std::forward_list<std::shared_ptr<SelectionNode>> m_callers;
    std::forward_list<SelectionFlags> m_flags;
};

} // namespace mol::internal

#endif // SELECTIONSTACK_HPP
