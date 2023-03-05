#ifndef PROPERTIES_HPP
#define PROPERTIES_HPP

#include "selections/SelectionNode.hpp"
#include "selections/numbers.hpp"
#include <vector>

namespace mol::internal
{

class PropSelection : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& data, Frame /*frame*/) const override;

protected:
    virtual bool selected(index_t atom_idx, MolData const& data) const = 0;
};

class NumPropSelection : public PropSelection
{
public:
    void add_number(SelNumber const& number)
    {
        m_numbers.push_back(number);
    }

    void add_range(SelNumberRange const& range)
    {
        m_ranges.push_back(range);
    }

    template <class T>
    requires std::integral<T> || std::floating_point<T>
    bool has(T const& value) const
    {
        for (auto const& range : m_ranges)
        {
            if (range.has(value))
            {
                return true;
            }
        }

        for (auto const &number : m_numbers)
        {
            if (number == value)
            {
                return true;
            }
        }
        return false;
    }

private:
    std::vector<SelNumber> m_numbers;
    std::vector<SelNumberRange> m_ranges;
};

class ResidSelection : public NumPropSelection
{
protected:
    bool selected(index_t atom_idx, MolData const& data) const override;
};

} // namespace mol::internal

#endif // PROPERTIES_HPP
