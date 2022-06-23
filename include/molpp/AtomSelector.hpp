#ifndef ATOMSELECTOR_HPP
#define ATOMSELECTOR_HPP

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
    AtomSelector(std::string const& selection, std::shared_ptr<internal::MolData> data);
    AtomSel apply(std::optional<size_t> frame);

private:
    void parse(std::string const& selection);

    std::shared_ptr<internal::MolData> m_data;
    std::shared_ptr<internal::SelectionNode> m_tree;
};

} // namespace mol

#endif // ATOMSELECTOR_HPP
