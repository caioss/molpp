#ifndef ATOMAGGREGATE_HPP
#define ATOMAGGREGATE_HPP

#include <molpp/internal/BaseAtomAggregate.hpp>
#include <memory>
#include <vector>
#include <concepts>

namespace mol::internal {

template <class Derived>
class AtomAggregate;

template <class T>
concept AtomAggregateDerived = requires(T t)
{
    std::derived_from<T, AtomAggregate<T>>;
    {t.atom_indices()} -> std::same_as<std::vector<index_t>>;
};

template <class Derived>
class AtomAggregate : public BaseAtomAggregate
{
public:
    using BaseAtomAggregate::coords_type;

    AtomAggregate() = delete;

    AtomAggregate(index_t const index, Frame const frame, internal::MolData* data)
    requires AtomAggregateDerived<Derived>
    : BaseAtomAggregate(index, frame, data)
    {}

    coords_type coords()
    {
        Derived &derived = static_cast<Derived &>(*this);
        return coords(derived.atom_indices());
    }

    std::vector<std::shared_ptr<Bond>> bonds()
    {
        Derived &derived = static_cast<Derived &>(*this);
        return bonds(derived.atom_indices());
    }

protected:
    using BaseAtomAggregate::coords;
    using BaseAtomAggregate::bonds;
    using BaseAtomAggregate::data;

    template <class, class>
    friend class Sel;
};

} // namespace mol::internal

#endif // ATOMAGGREGATE_HPP
