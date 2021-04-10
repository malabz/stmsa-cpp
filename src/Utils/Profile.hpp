#pragma once

#include <iterator>
#include <vector>

namespace utils
{

    template <unsigned char max_element>
    class Profile
    {
    public:
        using self_type = Profile<max_element>;

        template<typename InputIterator>
        Profile(InputIterator first, InputIterator last, size_t index) noexcept
        {
            constexpr size_t len = max_element + 1;
            
            unsigned frequencies[len];
            memset(frequencies, 0, sizeof(frequencies));

            size_t n = std::distance(first, last);
            for (size_t i = 0; i != n; ++i, ++first) ++frequencies[(*first)[index]];

            for (size_t i = 0; i != len; ++i)
                _relative_frequencies[i] = static_cast<double>(frequencies[i]) / n;

            _master = 0;
            for (unsigned char i = 1; i != len; ++i)
                if (_relative_frequencies[i] > _relative_frequencies[_master])
                    _master = i;
        }

        Profile(const self_type& rhs) noexcept = default;

        Profile(self_type&& rhs) noexcept = default;

        self_type& operator=(const self_type& rhs) noexcept = default;

        self_type& operator=(self_type&& rhs) noexcept = default;

        Profile()
        {
            std::fill(_relative_frequencies, _relative_frequencies + max_element + 1, 2.);
            _master = -1;
        }

        ~Profile() = default;

        double differs_from(const self_type& rhs) const noexcept
        {
            constexpr size_t len = max_element + 1;

            double difference = 0;
            for (size_t i = 0; i != len; ++i)
                difference += std::abs(_relative_frequencies[i] - rhs._relative_frequencies[i]);

            return difference;
        }

        operator size_t() const noexcept
        {
            return _master;
        }

    private:
        double _relative_frequencies[max_element + 1];
        unsigned char _master;
    };

    template <typename CharMatrix>
    std::vector<unsigned char> abstract(const CharMatrix& matrix, size_t row, size_t clm, size_t max_element);

    template <typename SuffixTreeType, typename RandomAccessIterator1, typename RandomAccessIterator2>
    std::vector<std::array<size_t, 3>> get_similar_substrings(const SuffixTreeType& suffix_tree, RandomAccessIterator1 lhs_first,
            RandomAccessIterator2 rhs_first, RandomAccessIterator2 rhs_last, double tolerance, size_t threshold);

    template <typename SuffixTreeType, typename RandomAccessIterator1, typename RandomAccessIterator2>
    std::vector<size_t> _search_for_prefix(const SuffixTreeType& suffix_tree, RandomAccessIterator1 lhs_first,
            RandomAccessIterator2 rhs_first, RandomAccessIterator2 rhs_last, double tolerance, size_t threshold);

}

template <typename CharMatrix>
std::vector<unsigned char> utils::abstract(const CharMatrix& matrix, size_t row, size_t clm, size_t max_element)
{
    std::vector<unsigned char> result(clm);

    size_t count_size = sizeof(size_t) * (max_element + 1);
    size_t* count = new size_t[max_element + 1];
    memset(count, 0, count_size);

    for (size_t j = 0; j != clm; ++j)
    {
        for (size_t i = 0; i != row; ++i) ++count[matrix[i][j]];

        unsigned char master = 0;
        for (size_t k = 1; k != max_element + 1; ++k)
            if (count[k] > count[master]) master = k;

        result[j] = master;
        memset(count, 0, count_size);
    }

    delete[] count;
    return result;
}

template <typename SuffixTreeType, typename RandomAccessIterator1, typename RandomAccessIterator2>
std::vector<std::array<size_t, 3>> utils::get_similar_substrings(const SuffixTreeType& suffix_tree, RandomAccessIterator1 lhs_first,
        RandomAccessIterator2 rhs_first, RandomAccessIterator2 rhs_last, double tolerance, size_t threshold)
{
        std::vector<std::array<size_t, 3>> similar_substrings;
        const size_t rhs_len = rhs_last - rhs_first;

        for (size_t rhs_index = 0; rhs_index < rhs_len; )
        {
            auto found = _search_for_prefix(suffix_tree, lhs_first, rhs_first + rhs_index, rhs_last, tolerance, threshold);

            if (found.size() == 0)
            {
                ++rhs_index;
            }
            else
            {
                for (size_t i = 1; i != found.size(); ++i)
                    similar_substrings.push_back(std::array<size_t, 3>({ found[i], rhs_index, found[0] }));

                rhs_index += found[0] - threshold + 1;
            }
        }

        return similar_substrings;
}

template <typename SuffixTreeType, typename RandomAccessIterator1, typename RandomAccessIterator2>
std::vector<size_t> utils::_search_for_prefix(const SuffixTreeType& suffix_tree, RandomAccessIterator1 lhs_first,
        RandomAccessIterator2 rhs_first, RandomAccessIterator2 rhs_last, double tolerance, size_t threshold)
{
    using Node = typename SuffixTreeType::Node;
    using value_type = typename SuffixTreeType::value_type;
    if (tolerance <= 0) return std::vector<size_t>();

    size_t common_prefix_length = 0;
    for (Node* last_node = suffix_tree.root, * curr_node = last_node->children[*rhs_first]; ;
            last_node = curr_node, curr_node = curr_node->children[*rhs_first]) // unsigned operator
    {
        if (curr_node == nullptr)
            return common_prefix_length < threshold ? std::vector<size_t>() : suffix_tree.get_all_beginning_with(last_node, common_prefix_length);

        const double curr_tolerance = tolerance / (common_prefix_length + curr_node->length);
        for (RandomAccessIterator1 lhs_begin = lhs_first + curr_node->first, lhs_end = lhs_begin + curr_node->length;
                lhs_begin != lhs_end && rhs_first != rhs_last; ++lhs_begin, ++rhs_first, ++common_prefix_length)
            if (lhs_begin->differs_from(*rhs_first) > curr_tolerance)
                return common_prefix_length < threshold ? std::vector<size_t>() : suffix_tree.get_all_beginning_with(curr_node, common_prefix_length);

        if (curr_node->children == nullptr || rhs_first == rhs_last)
            return common_prefix_length < threshold ? std::vector<size_t>() : suffix_tree.get_all_beginning_with(curr_node, common_prefix_length);
    }

}
