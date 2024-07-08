#ifndef MOLPP_SELECTIONS_NUMBERSET_HPP
#define MOLPP_SELECTIONS_NUMBERSET_HPP

#include "selections/numbers.hpp"

#include <vector>

namespace mol::internal
{

class NumberSet
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

    template<class T>
    requires std::integral<T> || std::floating_point<T>
    bool has(T const value) const
    {
        for (auto const& range : m_ranges)
        {
            if (range.has(value))
            {
                return true;
            }
        }

        for (auto const& number : m_numbers)
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

} // namespace mol::internal

#endif // MOLPP_SELECTIONS_NUMBERSET_HPP
