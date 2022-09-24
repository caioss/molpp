#ifndef SELECTIONS_HPP
#define SELECTIONS_HPP

#include <molpp/MolppCore.hpp>
#include "tools/math.hpp"
#include <set>
#include <memory>
#include <vector>
#include <string>
#include <concepts>
#include <optional>
#include <forward_list>

namespace mol::internal {

class MolData;
class SelectionNode;

class SelNumber
{
public:
    SelNumber(std::string const& text)
    : m_value(std::stod(text))
    {}

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
    : m_first(std::stoi(first)),
      m_last(std::stoi(last))
    {}

    bool has(std::integral auto const& query) const
    {
        return (query >= m_first) && (query <= m_last);
    }

    template <std::floating_point T>
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

struct SelectionFlags
{
    SelectionFlags(std::shared_ptr<std::set<index_t>> prev_mask, std::shared_ptr<std::set<index_t>> prev_selected)
    : mask{prev_mask},
      selected{prev_selected}
    {}

    SelectionFlags()
    : mask{std::make_shared<std::set<index_t>>()},
      selected{std::make_shared<std::set<index_t>>()}
    {}

    // Flags must be ordered
    std::shared_ptr<std::set<index_t>> mask;
    std::shared_ptr<std::set<index_t>> selected;
};

class SelectionStack
{
public:
    SelectionStack(std::shared_ptr<SelectionNode> root)
    : m_root(root)
    {}

    void push_flags(SelectionFlags const& flags)
    {
        m_flags.push_front(flags);
    }

    void push(std::shared_ptr<SelectionNode> node, SelectionFlags const& flags)
    {
        m_callers.push_front(node);
        push_flags(flags);
    }

    SelectionFlags pop_flags()
    {
        SelectionFlags flags = m_flags.front();
        m_flags.pop_front();
        return flags;
    }

    void evaluate(MolData const& data, SelectionFlags& flags, Frame frame);

private:
    std::shared_ptr<SelectionNode> m_root;
    std::forward_list<std::shared_ptr<SelectionNode>> m_callers;
    std::forward_list<SelectionFlags> m_flags;
};

class SelectionNode
{
public:
    virtual ~SelectionNode() {}
    virtual void evaluate(SelectionStack& stack, MolData const& data, Frame frame) const = 0;

    std::shared_ptr<SelectionNode> left;
    std::shared_ptr<SelectionNode> right;
};

class OrSelection : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& data, Frame frame) const override;
};

class AndSelection : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& data, Frame frame) const override;
};

class NotSelection : public SelectionNode
{
public:
    void evaluate(SelectionStack& stack, MolData const& /*data*/, Frame frame) const override;
};

class NumPropSelection : public SelectionNode
{
public:
    void add_number(SelNumber const& number);
    void add_range(SelNumberRange const& range);

    template <class T>
    requires std::integral<T> || std::floating_point<T>
    bool has(T const& value) const
    {
        for (auto const &range : m_ranges)
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
public:
    void evaluate(SelectionStack& stack, MolData const& data, Frame frame) const override;
};

} // namespace mol::internal

#endif // SELECTIONS_HPP
