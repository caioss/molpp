#ifndef MOLPP_SELECTIONS_SELECTIONINDICES_HPP
#define MOLPP_SELECTIONS_SELECTIONINDICES_HPP

#include <molpp/Common.hpp>

#include <memory>
#include <unordered_set>

namespace mol::internal
{

struct SelectionIndices
{
    using IndexSet = std::shared_ptr<std::unordered_set<index_t>>;

    SelectionIndices(IndexSet in_available, IndexSet in_selected)
    : available{in_available}
    , selected{in_selected}
    {}

    SelectionIndices()
    : available{make_index_set()}
    , selected{make_index_set()}
    {}

    static IndexSet make_index_set()
    {
        return std::make_shared<std::unordered_set<index_t>>();
    }

    IndexSet available;
    IndexSet selected;

    friend bool operator==(SelectionIndices const& lhs, SelectionIndices const& rhs)
    {
        return lhs.available == rhs.available && lhs.selected == rhs.selected;
    }
};

} // namespace mol::internal

#endif // MOLPP_SELECTIONS_SELECTIONINDICES_HPP
