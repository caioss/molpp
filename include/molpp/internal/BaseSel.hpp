#ifndef BASESEL_HPP
#define BASESEL_HPP

#include <molpp/MolError.hpp>
#include <memory>
#include <numeric>
#include <algorithm>

namespace mol::internal {

class MolData;

bool check_valid_frame(size_t const frame, std::shared_ptr<MolData> data);

template <class Type>
class BaseSel
{
private:
    template <class ItType>
    class Iterator;

public:
    using value_type = Type;
    using iterator = Iterator<Type>;
    using const_iterator = Iterator<const Type>;

    BaseSel() = delete;
    BaseSel(size_t const total_size, std::shared_ptr<MolData> data)
    : m_data(data),
      m_frame { 0 },
      m_selected(total_size, true),
      m_indices(total_size)
    {
        std::iota(m_indices.begin(), m_indices.end(), 0);
    }

    BaseSel(size_t const total_size, std::vector<size_t> const &indices, std::shared_ptr<MolData> data)
    : m_data(data),
      m_frame { 0 },
      m_selected(total_size, false),
      m_indices(indices)
    {
        update_indices(total_size);
    }

    BaseSel(size_t const total_size, std::vector<size_t> &&indices, std::shared_ptr<MolData> data)
    : m_data(data),
      m_frame { 0 },
      m_selected(total_size, false),
      m_indices(0)
    {
        std::swap(m_indices, indices);
        update_indices(total_size);
    }

    iterator begin()
    {
        return iterator(m_data, m_indices.begin(), m_frame);
    }

    iterator end()
    {
        return iterator(m_data, m_indices.end(), m_frame);
    }

    const_iterator cbegin()
    {
        return const_iterator(m_data, m_indices.begin(), m_frame);
    }

    const_iterator cend()
    {
        return const_iterator(m_data, m_indices.end(), m_frame);
    }

    Type operator[](size_t const index)
    {
        return Type(m_indices[index], m_frame, m_data);
    }

    Type at(size_t const index)
    {
        if (index >= m_indices.size())
        {
            throw mol::MolError("Out of bounds index: " + std::to_string(index));
        }

        return Type(m_indices[index], m_frame, m_data);
    }

    size_t size() const
    {
        return m_indices.size();
    }

    size_t frame() const
    {
        return m_frame;
    }

    void set_frame(size_t const frame)
    {
        if (!check_valid_frame(frame, m_data))
        {
            throw mol::MolError("Out of bounds frame: " + std::to_string(frame));
        }
        m_frame = frame;
    }

    bool contains(size_t const index) const
    {
        if (index >= m_selected.size())
        {
            return false;
        }
        else
        {
            return m_selected[index];
        }
    }

    std::vector<size_t> const &indices() const
    {
        return m_indices;
    }

    std::vector<bool> const &selected() const
    {
        return m_selected;
    }

protected:
    std::shared_ptr<MolData> m_data;

private:
    void update_indices(size_t const total_size)
    {
        std::sort(m_indices.begin(), m_indices.end());
        for (size_t index : m_indices)
        {
            if (index >= total_size)
            {
                throw mol::MolError("Out of bounds selection index: " + std::to_string(index));
            }

            m_selected[index] = true;
        }
    }

    size_t m_frame;
    std::vector<bool> m_selected;
    std::vector<size_t> m_indices;

    template <class ItType>
    class Iterator
    {
    private:
        using indices_iterator = std::vector<size_t>::iterator;

    public:
        using iterator_category = indices_iterator::iterator_category;
        using difference_type = indices_iterator::difference_type;
        using value_type = ItType;
        using pointer = ItType *;
        using reference = ItType &;

        Iterator(std::shared_ptr<MolData> data, indices_iterator begin, size_t frame)
        : m_frame { frame },
          m_data(data),
          m_current(begin)
        {}

        value_type operator*() const
        {
            return ItType(*m_current, m_frame, m_data);
        }

        Iterator& operator++()
        {
            m_current++;
            return *this;
        }

        Iterator operator++(int)
        {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        difference_type operator-(const Iterator& other)
        {
            return m_current - other.m_current;
        }

        bool operator==(const Iterator& other)
        {
            return m_current == other.m_current;
        };

        bool operator!=(const Iterator& other)
        {
            return m_current != other.m_current;
        };

    private:
        size_t m_frame;
        std::shared_ptr<MolData> m_data;
        indices_iterator m_current;
    };
};

} // namespace mol::internal

#endif // BASESEL_HPP
