#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include <set>
#include <queue>
#include <concepts>
#include <unordered_map>

namespace mol::internal {

template <class Container>
class BreadthFirstTraversal
{
public:
    using node_type = typename Container::node_type;

    template <std::invocable<node_type const&> StopPredicate, std::invocable<node_type const&> FilterPredicate>
    bool run(Container const &container, node_type const &start, StopPredicate stop, FilterPredicate filter)
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

            if (stop(current))
            {
                return true;
            }

            for (node_type const &node : container.adjacency(current))
            {
                if (!filter(node))
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

    std::set<node_type> const &visited() const
    {
        return m_visited;
    }

    std::unordered_map<node_type, node_type> const &parent_map() const
    {
        return m_parent;
    }

private:
    std::set<node_type> m_visited;
    std::unordered_map<node_type, node_type> m_parent;
};

} // namespace mol::internal

#endif // ALGORITHMS_HPP
