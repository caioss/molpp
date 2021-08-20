#ifndef ATOMAGGREGATE_HPP
#define ATOMAGGREGATE_HPP

#include <molpp/internal/BaseAtomAggregate.hpp>
#include <memory>
#include <vector>
#include <concepts>

namespace mol {

class AtomSel;

namespace internal {

template <class Derived>
class AtomAggregate;

template <class T>
concept is_AtomAggregate = requires(T t)
{
    std::derived_from<T, AtomAggregate<T>>;
    {t.atom_indices()} -> std::same_as<std::vector<size_t>>;
};

template <class Derived>
class AtomAggregate : public BaseAtomAggregate
{
public:
    using BaseAtomAggregate::coords_type;

    AtomAggregate() = delete;
    using BaseAtomAggregate::BaseAtomAggregate;

    coords_type coords()
    requires is_AtomAggregate<Derived>
    {
        Derived &derived = static_cast<Derived &>(*this);
        return coords(derived.atom_indices());
    }

    std::vector<std::shared_ptr<Bond>> bonds()
    requires is_AtomAggregate<Derived>
    {
        Derived &derived = static_cast<Derived &>(*this);
        return bonds(derived.atom_indices());
    }

    std::shared_ptr<AtomSel> bonded()
    requires is_AtomAggregate<Derived>
    {
        Derived &derived = static_cast<Derived &>(*this);
        return bonded(derived.atom_indices());
    }

    std::shared_ptr<AtomSel> atoms()
    requires is_AtomAggregate<Derived>
    {
        Derived &derived = static_cast<Derived &>(*this);
        return atoms(derived.atom_indices());
    }

protected:
    using BaseAtomAggregate::coords;
    using BaseAtomAggregate::bonds;
    using BaseAtomAggregate::bonded;
    using BaseAtomAggregate::atoms;
    using BaseAtomAggregate::data;
    using BaseAtomAggregate::cdata;
};

} // namespace internal
} // namespace mol

#endif // ATOMAGGREGATE_HPP
