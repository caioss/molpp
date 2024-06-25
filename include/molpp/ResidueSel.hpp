#ifndef MOLPP_RESIDUESEL_HPP
#define MOLPP_RESIDUESEL_HPP

#include <molpp/internal/Sel.hpp>
#include <molpp/Residue.hpp>

namespace mol
{

class ResidueSel;

namespace internal
{

template<>
struct SelTraits<ResidueSel>
{
    using entity_type = Residue;
};

} // namespace internal

class ResidueSel : public internal::Sel<ResidueSel>
{
public:
    ResidueSel() = delete;
    using internal::Sel<ResidueSel>::Sel;

    template<class>
    friend class internal::Sel;
};

} // namespace mol

#endif // MOLPP_RESIDUESEL_HPP
