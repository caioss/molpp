#ifndef MOLPP_SEGMENTSEL_HPP
#define MOLPP_SEGMENTSEL_HPP

#include <molpp/internal/Sel.hpp>
#include <molpp/Segment.hpp>

namespace mol
{

class SegmentSel;

namespace internal
{

template<>
struct SelTraits<SegmentSel>
{
    using entity_type = Segment;
};

} // namespace internal

class SegmentSel : public internal::Sel<SegmentSel>
{
public:
    SegmentSel() = delete;
    using internal::Sel<SegmentSel>::Sel;

    template<class>
    friend class internal::Sel;
};

} // namespace mol

#endif // MOLPP_SEGMENTSEL_HPP
