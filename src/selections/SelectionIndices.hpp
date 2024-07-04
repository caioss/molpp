#ifndef MOLPP_SELECTIONS_SELECTIONINDICES_HPP
#define MOLPP_SELECTIONS_SELECTIONINDICES_HPP

#include <molpp/Common.hpp>

#include <memory>
#include <unordered_set>

namespace mol::internal
{

struct SelectionIndices
{
    SelectionIndices(std::shared_ptr<std::unordered_set<index_t>> in_available, std::shared_ptr<std::unordered_set<index_t>> in_selected)
    : available{in_available}
    , selected{in_selected}
    {}

    SelectionIndices()
    : available{std::make_shared<std::unordered_set<index_t>>()}
    , selected{std::make_shared<std::unordered_set<index_t>>()}
    {}

    std::shared_ptr<std::unordered_set<index_t>> available;
    std::shared_ptr<std::unordered_set<index_t>> selected;

    friend bool operator==(SelectionIndices const& lhs, SelectionIndices const& rhs)
    {
        return lhs.available == rhs.available && lhs.selected == rhs.selected;
    }
};

} // namespace mol::internal

#endif // MOLPP_SELECTIONS_SELECTIONINDICES_HPP
