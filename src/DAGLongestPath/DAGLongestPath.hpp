#include <vector>

namespace dag_longest_path
{

    template <typename MatrixType>
    std::vector<size_t> longest_path_of(const MatrixType& graph, size_t node_num)
    {
        auto distance    = new size_t[node_num]();
        auto path_record = new size_t[node_num];

        for (size_t to = 1; to != node_num; ++to)
            for (size_t station = 0; station != to; ++station)
                if (graph[station][to])
                {
                    unsigned curr = distance[station] + graph[station][to];
                    if (curr > distance[to])
                    {
                        distance[to] = curr;
                        path_record[to] = station;
                    }
                }

        size_t end_node = 0;
        for (size_t i = 1; i != node_num; ++i)
            if (distance[end_node] < distance[i]) end_node = i;

        std::vector<size_t> path;
        for (size_t i = end_node; i != 0; i = path_record[i])
            path.push_back(i);

        for (size_t i = 0, mid = path.size() >> 1; i != mid; ++i)
            std::swap(path[i], path[path.size() - i - 1]);

        // std::copy(distance, distance + node_num, std::ostream_iterator<size_t>(std::cout, ", ")); std::cout << '\n';
        // std::copy(path_record, path_record + node_num, std::ostream_iterator<size_t>(std::cout, ", ")); std::cout << '\n';
        // std::copy(path.cbegin(), path.cend(), std::ostream_iterator<size_t>(std::cout, ", ")); std::cout << '\n';
            
        delete[] distance;
        delete[] path_record;

        return path;
    }
    
}
