#ifndef MOLPP_SELECTIONS_PROPERTIES_HPP
#define MOLPP_SELECTIONS_PROPERTIES_HPP

#include "selections/PropertySelection.hpp"
#include "selections/NumberSet.hpp"

namespace mol::internal
{

class ResidSelection : public PropertySelection<ResidSelection>
{
public:
    ResidSelection(NumberSet&& resids);
    bool evaluate_atom(index_t const atom_index, MolData const& data) const;

private:
    NumberSet m_resids;
};

} // namespace mol::internal

#endif // MOLPP_SELECTIONS_PROPERTIES_HPP
