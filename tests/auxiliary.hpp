#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include "core/MolData.hpp"
#include <vector>

using namespace mol;
using namespace mol::internal;

MolData create_moldata(size_t const num_res, size_t const num_res_atoms, size_t const num_chains, size_t const num_segments, size_t const num_frames);

template <class Type>
std::vector<Type> view2vector(auto const &view)
{
    return {view.begin(), view.end()};
}

#endif // AUXILIARY_HPP
