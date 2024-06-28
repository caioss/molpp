#ifndef MOLPP_CHAINSEL_HPP
#define MOLPP_CHAINSEL_HPP

#include <molpp/internal/Sel.hpp>
#include <molpp/Chain.hpp>

namespace mol
{

class ChainSel;

namespace internal
{

template<>
struct SelTraits<ChainSel>
{
    using entity_type = Chain;
};

} // namespace internal

class ChainSel : public internal::Sel<ChainSel>
{
public:
    ChainSel() = delete;
    using internal::Sel<ChainSel>::Sel;

    template<class>
    friend class internal::Sel;
};

} // namespace mol

#endif // MOLPP_CHAINSEL_HPP
