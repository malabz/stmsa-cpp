#include "GraphTester.hpp"
#include "../Utils/Graph.hpp"

#include <iostream>

void utils::test_graph()
{
    constexpr size_t node_num = 4;

    utils::AdjacencyList graph(node_num);
    graph.add_edge(0, 1, 2);
    graph.add_edge(0, 3, 5);
    graph.add_edge(1, 2, 3);
    graph.add_edge(3, 2, 1);

    auto result = graph.get_longest_path();
    std::copy(result.cbegin(), result.cend(), std::ostream_iterator<size_t>(std::cout, " ")); std::cout << '\n';
}
