#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <molpp/internal/MolData.hpp>

#include <vector>

using namespace mol;
using namespace mol::internal;

MolData create_moldata(size_t const num_res, size_t const num_res_atoms, size_t const num_chains, size_t const num_segments, size_t const num_frames);

template<std::ranges::range Range>
constexpr auto view2vector(Range&& r)
{
    using elem_t = std::decay_t<std::ranges::range_value_t<Range>>;
    return std::vector<elem_t>{r.begin(), r.end()};
}

#endif // AUXILIARY_HPP
