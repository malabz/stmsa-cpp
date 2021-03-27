#include "StarAligner.hpp"
#include "../Utils/Pseudo.hpp"
#include "../Utils/Utils.hpp"
#include "../PairwiseAlignment/PairwiseAlignment.hpp"
#include "../DAGLongestPath/DAGLongestPath.hpp"

#include <numeric>

std::vector<std::vector<unsigned char>> star_alignment::StarAligner::align(const std::vector<sequence_type>& sequences)
{
    return StarAligner(sequences)._align();
}

star_alignment::StarAligner::StarAligner(const std::vector<sequence_type>& sequences):
        _sequences(sequences),
        _row(_sequences.size()),
        _lengths(_set_lengths()),
        _centre(_set_centre()),
        _centre_len(_centre.size())
{}

std::vector<size_t> star_alignment::StarAligner::_set_lengths() const
{
    std::vector<size_t> lengths(_row);

    for (size_t i = 0; i != _row; ++i) lengths[i] = _sequences[i].size();
    return lengths;
}

// to do
std::vector<unsigned char> star_alignment::StarAligner::_set_centre() const
{
    size_t centre_index = 0;

    // for (size_t i = 1; i != _row; ++i)
    //     if (_lengths[i] > _lengths[centre_index])
    //         centre_index = i;
    return _sequences[centre_index];
}

std::vector<std::vector<unsigned char>> star_alignment::StarAligner::_align() const
{
    return _insert_gaps(_merge_results(_pairwise_align()));
}

auto star_alignment::StarAligner::_pairwise_align() const -> std::vector<std::array<std::vector<indel>, 2>>
{
    constexpr size_t threshold = 15;

    suffixtree::SuffixTree<unsigned char> st(_centre.cbegin(), _centre.cend(), nucleic_acid_pseudo::MAX_ELE);
    std::vector<std::array<std::vector<indel>, 2>> all_pairwise_gaps;

    for (size_t i = 0; i != _row; ++i)
    {
        auto identical_substrings = _optimal_path(st.get_identical_substrings(_sequences[i].cbegin(), _sequences[i].cend(), threshold));

        std::vector<quadra> intervals;
        intervals.reserve(identical_substrings.size() + 1);

        if (identical_substrings[0][0] || identical_substrings[0][1])
            intervals.push_back(quadra
                    ({ 0, identical_substrings[0][0], 0, identical_substrings[0][1] }));

        for (size_t j = 0, end_index = identical_substrings.size() - 1; j != end_index; ++j)
            if (identical_substrings[j][0] + identical_substrings[j][2] != identical_substrings[j + 1][0] ||
                identical_substrings[j][1] + identical_substrings[j][2] != identical_substrings[j + 1][1])
                intervals.push_back(quadra
                ({
                    identical_substrings[j][0] + identical_substrings[j][2], identical_substrings[j + 1][0],
                    identical_substrings[j][1] + identical_substrings[j][2], identical_substrings[j + 1][1]
                }));

        if (identical_substrings.back()[0] + identical_substrings.back()[2] != _centre_len || 
            identical_substrings.back()[1] + identical_substrings.back()[2] != _lengths[i])
            intervals.push_back(quadra
            ({
                identical_substrings.back()[0] + identical_substrings.back()[2], _centre_len,
                identical_substrings.back()[1] + identical_substrings.back()[2], _lengths[i]
            }));
        
        std::array<std::vector<indel>, 2> pairwise_gaps;
        for (size_t j = 0; j != intervals.size(); ++j)
        {
            const size_t centre_begin = intervals[j][0];
            const size_t centre_end = intervals[j][1];
            const size_t sequence_begin = intervals[j][2];
            const size_t sequence_end = intervals[j][3];

            auto [lhs_gaps, rhs_gaps] = pairwise_alignment::needleman_wunsh( _centre.cbegin() + centre_begin, _centre.cbegin() + centre_end,
                    _sequences[i].cbegin() + sequence_begin, _sequences[i].cbegin() + sequence_end, pairwise_alignment::default_scoring_matrix);

            _converse(lhs_gaps, pairwise_gaps[0], centre_begin);
            _converse(rhs_gaps, pairwise_gaps[1], sequence_begin);
        }
        all_pairwise_gaps.push_back(pairwise_gaps);
    }

    return all_pairwise_gaps;
}

auto star_alignment::StarAligner::_optimal_path(const std::vector<triple>& identical_substrings)
        -> std::vector<triple> 
{
    if (identical_substrings.empty()) return std::vector<triple>();

    const size_t pair_num = identical_substrings.size();
    size_t** graph = new size_t*[pair_num + 1];
    for (size_t i = 0; i != pair_num + 1; ++i) graph[i] = new size_t[pair_num + 1]();

    for (size_t i = 0; i != pair_num; ++i)
    for (size_t j = 0; j != pair_num; ++j)
        if (i != j && identical_substrings[i][0] + identical_substrings[i][2] < identical_substrings[j][0] + identical_substrings[j][2]
                   && identical_substrings[i][1] + identical_substrings[i][2] < identical_substrings[j][1] + identical_substrings[j][2])
        {
            const int possible_overlap = std::max(
                    static_cast<int>(identical_substrings[i][0] + identical_substrings[i][2]) - static_cast<int>(identical_substrings[j][0]), 
                    static_cast<int>(identical_substrings[i][1] + identical_substrings[i][2]) - static_cast<int>(identical_substrings[j][1]));

            graph[i + 1][j + 1] = identical_substrings[j][2];
            if (possible_overlap > 0) graph[i + 1][j + 1] -= possible_overlap;
        }

    for (size_t i = 0; i != pair_num; ++i)
        graph[0][i + 1] = identical_substrings[i][2];

    const auto optimal_path = dag_longest_path::longest_path_of(graph, pair_num + 1);

    std::vector<triple> optimal_identical_substrings;
    optimal_identical_substrings.reserve(optimal_path.size());
    optimal_identical_substrings.push_back(triple( { identical_substrings[optimal_path[0] - 1][0],
                                                     identical_substrings[optimal_path[0] - 1][1],
                                                     identical_substrings[optimal_path[0] - 1][2] }));

    for (size_t i = 0; i < optimal_path.size() - 1; ++i)
    {
        size_t new_len = graph[optimal_path[i]][optimal_path[i + 1]];
        size_t old_len = identical_substrings[optimal_path[i + 1] - 1][2];
        int difference = static_cast<int>(old_len) - static_cast<int>(new_len);

        size_t lhs_first = identical_substrings[optimal_path[i + 1] - 1][0];
        size_t rhs_first = identical_substrings[optimal_path[i + 1] - 1][1];
        if (difference > 0)
        { lhs_first += difference; rhs_first += difference; }

        optimal_identical_substrings.push_back(triple({ lhs_first, rhs_first, new_len }));
    }

    for (size_t i = 0; i != pair_num + 1; ++i) delete[] graph[i];
    delete[] graph;

    return optimal_identical_substrings;
}

void star_alignment::StarAligner::_converse(const std::vector<size_t>& src_gaps, std::vector<indel>& des_gaps, size_t start)
{
    for (size_t i = 0; i != src_gaps.size(); ++i)
        if (src_gaps[i])
        {
            if (des_gaps.size() && des_gaps.back().index == start + i)
                des_gaps.back().number += src_gaps[i];
            else
                des_gaps.push_back(indel({ start + i, src_gaps[i] }));
        }
}

auto star_alignment::StarAligner::_merge_results(const std::vector<std::array<std::vector<indel>, 2>>& pairwise_gaps) const
        -> std::vector<std::vector<indel>>
{
    std::vector<indel> final_centre_gaps;
    for (size_t i = 0; i != _row; ++i)
    {
        const auto& curr_centre_gaps = pairwise_gaps[i][0];
        for (size_t lhs_pointer = 0, rhs_pointer = 0; rhs_pointer != curr_centre_gaps.size(); )
        {
            if (lhs_pointer == final_centre_gaps.size())
            {
                final_centre_gaps.insert(final_centre_gaps.cend(), curr_centre_gaps.cbegin() + rhs_pointer, curr_centre_gaps.cend());
                break;
            }

            if (final_centre_gaps[lhs_pointer].index == curr_centre_gaps[rhs_pointer].index)
            {
                if (final_centre_gaps[lhs_pointer].number < curr_centre_gaps[rhs_pointer].number)
                    final_centre_gaps[lhs_pointer].number = curr_centre_gaps[rhs_pointer].number;
                ++lhs_pointer;
                ++rhs_pointer;
            }
            else if (final_centre_gaps[lhs_pointer].index < curr_centre_gaps[rhs_pointer].index)
            {
                ++lhs_pointer;
            }
            else
            {
                final_centre_gaps.insert(final_centre_gaps.cbegin() + lhs_pointer, curr_centre_gaps[rhs_pointer]);
                ++lhs_pointer; // because of the insert above
                ++rhs_pointer;
            }
        }
    }

    std::vector<std::vector<indel>> final_sequence_gaps;
    final_sequence_gaps.reserve(_row);
    for (size_t i = 0; i != _row; ++i)
    {
        const auto& curr_centre_gaps = pairwise_gaps[i][0];
        const auto& curr_sequence_gaps = pairwise_gaps[i][1];

        std::vector<indel> centre_addition = _minus(final_centre_gaps, curr_centre_gaps);

        std::vector<indel> sequence_addition;
        for (size_t centre_index = 0, sequence_index = 0, centre_gaps_index = 0, sequence_gaps_index = 0, centre_addition_index = 0;
                centre_addition_index != centre_addition.size(); ++centre_addition_index)
        {
            const auto curr_addition = centre_addition[centre_addition_index]; // current addition pending process

            while (centre_index < curr_addition.index)
            {
                size_t centre_distance = centre_gaps_index < curr_centre_gaps.size() ?
                        curr_centre_gaps[centre_gaps_index].index - centre_index : std::numeric_limits<size_t>::max();
                size_t sequence_distance = sequence_gaps_index < curr_sequence_gaps.size() ?
                        curr_sequence_gaps[sequence_gaps_index].index - sequence_index : std::numeric_limits<size_t>::max();

                size_t step = std::min({ sequence_distance, centre_distance, curr_addition.index - centre_index }); // assure centre_index <= curr_addtion.index
                centre_index += step;
                sequence_index += step;

                if (centre_gaps_index < curr_centre_gaps.size() && curr_centre_gaps[centre_gaps_index].index == centre_index)
                    sequence_index += curr_centre_gaps[centre_gaps_index++].number;

                else if (sequence_gaps_index < curr_sequence_gaps.size() && curr_sequence_gaps[sequence_gaps_index].index == sequence_index)
                    centre_index += curr_sequence_gaps[sequence_gaps_index++].number;
            }

            if (sequence_addition.size() && sequence_index == sequence_addition.back().index)
                sequence_addition.back().number += curr_addition.number;
            else
                sequence_addition.push_back(indel({ sequence_index, curr_addition.number }));
        }

        final_sequence_gaps.push_back(_add(curr_sequence_gaps, sequence_addition));
    }

    return final_sequence_gaps;
}

auto star_alignment::StarAligner::_add(const std::vector<indel>& lhs, const std::vector<indel>& rhs)
        -> std::vector<indel>
{
    std::vector<indel> sum;

    size_t lhs_index = 0, rhs_index = 0;
    for (; lhs_index != lhs.size() && rhs_index != rhs.size(); )
        if (lhs[lhs_index].index < rhs[rhs_index].index)
        {
            sum.push_back(lhs[lhs_index]);
            ++lhs_index;
        }
        else if (lhs[lhs_index].index > rhs[rhs_index].index)
        {
            sum.push_back(rhs[rhs_index]);
            ++rhs_index;
        }
        else
        {
            sum.push_back(indel({ lhs[lhs_index].index, lhs[lhs_index].number + rhs[rhs_index].number }));
            ++lhs_index;
            ++rhs_index;
        }

    sum.insert(sum.cend(), lhs.cbegin() + lhs_index, lhs.cend());
    sum.insert(sum.cend(), rhs.cbegin() + rhs_index, rhs.cend());
    return sum;
}

// assume lhs >= rhs
auto star_alignment::StarAligner::_minus(const std::vector<indel>& lhs, const std::vector<indel>& rhs)
        -> std::vector<indel>
{
    std::vector<indel> difference;
    difference.reserve(lhs.size());

    size_t lhs_index = 0;
    for (size_t rhs_index = 0; rhs_index != rhs.size(); ++lhs_index, ++rhs_index)
    {
        while (lhs[lhs_index].index != rhs[rhs_index].index)
            difference.push_back(lhs[lhs_index++]);

        size_t distance = lhs[lhs_index].number - rhs[rhs_index].number;
        if (distance) difference.push_back(indel({ lhs[lhs_index].index, distance }));
    }

    difference.insert(difference.cend(), lhs.cbegin() + lhs_index, lhs.cend());
    return difference;
}

auto star_alignment::StarAligner::_insert_gaps(const std::vector<std::vector<indel>>& gaps) const
        -> std::vector<sequence_type>
{
    size_t length = _lengths[0];
    for (const auto gap : gaps[0]) length += gap.number;

    std::vector<sequence_type> aligned;
    aligned.reserve(_row);

    for (size_t i = 0; i != _row; ++i)
        aligned.push_back(_insert_gaps(_sequences[i], gaps[i], length));
    return aligned;
}

auto star_alignment::StarAligner::_insert_gaps(const sequence_type& sequence, const std::vector<indel>& gaps, size_t length)
        -> sequence_type
{
    sequence_type gapped;
    gapped.reserve(length);

    size_t sequence_index = 0;
    for (const auto gap : gaps)
    {
        const size_t new_sequence_index = sequence_index + (gap.index - sequence_index);
        gapped.insert(gapped.cend(), sequence.cbegin() + sequence_index, sequence.cbegin() + new_sequence_index);
        gapped.insert(gapped.cend(), gap.number, nucleic_acid_pseudo::GAP);

        sequence_index = new_sequence_index;
    }
    gapped.insert(gapped.cend(), sequence.cbegin() + sequence_index, sequence.cend());
    return gapped;
}
