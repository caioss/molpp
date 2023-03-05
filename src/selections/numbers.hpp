#ifndef NUMBERS_HPP
#define NUMBERS_HPP

#include <string>
#include <concepts>
#include "tools/math.hpp"

namespace mol::internal
{

class SelNumber
{
public:
    SelNumber(std::string const& text)
    : m_value(std::stod(text))
    {
    }

    bool operator==(std::floating_point auto const& rhs) const
    {
        // High epsilon allowed by molecular properties
        return essentially_equal(m_value, rhs, 1e-4);
    }

    bool operator==(std::integral auto const& rhs) const
    {
        return m_value == rhs;
    }

private:
    double m_value;
};

class SelNumberRange
{
public:
    SelNumberRange(std::string const& first, std::string const& last)
    : m_first(std::stoi(first))
    , m_last(std::stoi(last))
    {
    }

    bool has(std::integral auto const& query) const
    {
        return (query >= m_first) && (query <= m_last);
    }

    template<std::floating_point T>
    bool has(T const& query) const
    {
        T const first = m_first;
        T const last = m_last;
        // High epsilon allowed by molecular properties
        return essentially_equal(first, query, 1e-4) ||
               essentially_equal(last, query, 1e-4) ||
               (query > first && query < last);
    }

private:
    int m_first;
    int m_last;
};

} // namespace mol::internal

#endif // NUMBERS_HPP
