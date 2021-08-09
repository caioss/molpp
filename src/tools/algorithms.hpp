#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace mol::internal {

template <class Container>
class BFS
{
public:
    using node_type = typename Container::node_type;

    bool run(Container const &container, node_type const &start, node_type const &end)
    {
        m_visited.clear();
        m_parent.clear();
        std::queue<node_type> queue;
        queue.push(start);
        m_visited.insert(start);

        while (!queue.empty())
        {
            node_type current = queue.front();
            queue.pop();

            if (current == end)
            {
                return true;
            }

            for (node_type const &node : container.adjacency(current))
            {
                if (m_mask.size() && !m_mask.count(node))
                {
                    continue;
                }

                if (!m_visited.count(node))
                {
                    m_visited.insert(node);
                    m_parent[node] = current;
                    queue.push(node);
                }
            }
        }

        return false;
    }

    std::unordered_set<node_type> const &mask() const
    {
        return m_mask;
    }

    void set_mask(std::unordered_set<node_type> const &mask)
    {
        m_mask = mask;
    }

    std::unordered_set<node_type> const &visited() const
    {
        return m_visited;
    }

    std::unordered_map<node_type, node_type> const &parent_map() const
    {
        return m_parent;
    }

private:
    std::unordered_set<node_type> m_mask;
    std::unordered_set<node_type> m_visited;
    std::unordered_map<node_type, node_type> m_parent;
};

} // namespace mol::internal

#endif // ALGORITHMS_HPP
