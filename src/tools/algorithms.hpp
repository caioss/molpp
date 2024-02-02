#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include <list>
#include <queue>
#include <ranges>
#include <concepts>
#include <unordered_set>
#include <unordered_map>

namespace mol::internal
{

template<class Container>
requires requires(Container obj) {
    typename Container::node_type;
    {
        obj.adjacency(std::declval<typename Container::node_type>())
    } -> std::ranges::range;
}
class BreadthFirstTraversal
{
public:
    using node_type = typename Container::node_type;

    BreadthFirstTraversal(Container const& container)
    : m_container(container)
    {}

    template<std::predicate<node_type> Stop, std::predicate<node_type> Filter>
    bool run(node_type const& start, Stop stop, Filter filter)
    {
        m_visited.clear();
        m_parent.clear();
        if (!filter(start))
        {
            return false;
        }
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

            for (node_type const& node : m_container.adjacency(current))
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

    std::unordered_set<node_type> const& visited() const
    {
        return m_visited;
    }

    std::unordered_map<node_type, node_type> const& parent_map() const
    {
        return m_parent;
    }

private:
    Container const& m_container;
    std::unordered_set<node_type> m_visited;
    std::unordered_map<node_type, node_type> m_parent;
};

template<class Container>
class ConnectedComponents
{
public:
    using node_type = typename Container::node_type;

    ConnectedComponents(Container const& container)
    : m_container(container)
    {}

    template<std::invocable<node_type const&> Predicate>
    size_t run(Predicate filter)
    {
        auto const container_nodes = m_container.nodes();
        std::unordered_set<node_type> not_visited(container_nodes.begin(), container_nodes.end());
        m_components.clear();
        BreadthFirstTraversal bfs(m_container);

        while (!not_visited.empty())
        {
            node_type const start = *(not_visited.begin());
            if (!filter(start))
            {
                not_visited.erase(start);
                continue;
            }

            bfs.run(start, [](auto){return false;}, filter);

            for (auto const& node : bfs.visited())
            {
                not_visited.erase(node);
            }

            m_components.push_back(bfs.visited());
        }

        return m_components.size();
    }

    std::list<std::unordered_set<node_type>> const& components()
    {
        return m_components;
    }

private:
    Container const& m_container;
    std::list<std::unordered_set<node_type>> m_components;
};

} // namespace mol::internal

#endif // ALGORITHMS_HPP
