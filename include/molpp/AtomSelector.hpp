#ifndef MOLPP_ATOMSELECTOR_HPP
#define MOLPP_ATOMSELECTOR_HPP

#include <molpp/Common.hpp>
#include <molpp/AtomSel.hpp>

#include <memory>

namespace mol
{

namespace internal
{
class MolData;
class SelectionNode;
} // namespace internal

class AtomSelector
{
public:
    AtomSelector() = delete;
    AtomSelector(std::string const& selection, internal::MolData& data);
    AtomSel apply(Frame frame);

private:
    void parse(std::string const& selection);

    internal::MolData& m_data;
    std::shared_ptr<internal::SelectionNode> m_tree;
};

} // namespace mol

#endif // MOLPP_ATOMSELECTOR_HPP
