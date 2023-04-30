#ifndef ATOMSELECTOR_HPP
#define ATOMSELECTOR_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/AtomSel.hpp>
#include <memory>

namespace mol {

namespace internal {
class MolData;
class SelectionNode;
}

class AtomSelector
{
public:
    AtomSelector() = delete;
    AtomSelector(std::string const& selection, internal::MolData* data);
    AtomSel apply(Frame frame);

private:
    void parse(std::string const& selection);

    internal::MolData* m_data;
    std::shared_ptr<internal::SelectionNode> m_tree;
};

} // namespace mol

#endif // ATOMSELECTOR_HPP
