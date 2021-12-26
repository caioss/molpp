#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include "core/MolData.hpp"
#include <vector>

using namespace mol;
using namespace mol::internal;

std::shared_ptr<MolData> aux_moldata();

template <class Type>
std::vector<Type> view2vector(auto const &view)
{
    return {view.begin(), view.end()};
}

#endif // AUXILIARY_HPP
