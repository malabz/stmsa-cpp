#include "StarAligner.hpp"
#include "../Utils/Pseudo.hpp"
#include "../Utils/Utils.hpp"
#include "../Utils/Graph.hpp"
#include "../PairwiseAlignment/NeedlemanWunshReusable.hpp"

std::vector<std::vector<unsigned char>> star_alignment::StarAligner::align(const std::vector<sequence_type> &sequences)
{
    return StarAligner(sequences)._align();
}

auto star_alignment::StarAligner::get_gaps(const std::vector<sequence_type> &sequences) -> std::vector<std::vector<utils::Insertion>>
{
    return StarAligner(sequences)._get_gaps();
}

star_alignment::StarAligner::StarAligner(const std::vector<sequence_type> &sequences)
        : _sequences(sequences)
        , _row(_sequences.size())
        , _lengths(_set_lengths())
        , _centre(_set_centre())
        , _centre_len(_centre.size())
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

auto star_alignment::StarAligner::_get_gaps() const
        -> std::vector<std::vector<utils::Insertion>>
{
    return _merge_results(_pairwise_align());
}

auto star_alignment::StarAligner::_pairwise_align() const -> std::vector<std::array<std::vector<utils::Insertion>, 2>>
{
    static constexpr size_t threshold = 15;

    suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER> st(_centre.cbegin(), _centre.cend(), nucleic_acid_pseudo::GAP);
    std::vector<std::array<std::vector<utils::Insertion>, 2>> all_pairwise_gaps;

    using iter = std::vector<unsigned char>::const_iterator;
    pairwise_alignment::NeedlemanWunshReusable<iter, iter, decltype(pairwise_alignment::default_scoring_matrix)>
    nw(pairwise_alignment::default_scoring_matrix, pairwise_alignment::default_gap_open, pairwise_alignment::default_gap_extention);

    for (size_t i = 0; i != _row; ++i)
    {
        auto identical_substrings = _optimal_path(st.get_identical_substrings(_sequences[i].cbegin(), _sequences[i].cend(), threshold));

        std::vector<quadra> intervals;
        intervals.reserve(identical_substrings.size() + 1);

        if (identical_substrings.empty())
        {
            intervals.push_back(quadra({ 0, _centre.size(), 0, _sequences[i].size() }));
        }
        else
        {
            if (identical_substrings[0][0] || identical_substrings[0][1])
                intervals.push_back(quadra({ 0, identical_substrings[0][0], 0, identical_substrings[0][1] }));

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
        }

        std::array<std::vector<utils::Insertion>, 2> pairwise_gaps;
        for (size_t j = 0; j != intervals.size(); ++j)
        {
            const size_t centre_begin = intervals[j][0];
            const size_t centre_end = intervals[j][1];
            const size_t sequence_begin = intervals[j][2];
            const size_t sequence_end = intervals[j][3];

            unsigned flag = 0;
            if (centre_begin == 0) flag |= pairwise_alignment::LEFT_ENDING;
            if (centre_end == _centre_len) flag |= pairwise_alignment::RIGHT_ENDING;

            auto [lhs_gaps, rhs_gaps] = nw(_centre.cbegin() + centre_begin, _centre.cbegin() + centre_end,
                    _sequences[i].cbegin() + sequence_begin, _sequences[i].cbegin() + sequence_end, flag);

            _append(lhs_gaps, pairwise_gaps[0], centre_begin);
            _append(rhs_gaps, pairwise_gaps[1], sequence_begin);
        }

        size_t sequence_gap_num = 0;
        for (auto gap : pairwise_gaps[0]) sequence_gap_num += gap.number;

        all_pairwise_gaps.push_back(pairwise_gaps);
    }

    return all_pairwise_gaps;
}

auto star_alignment::StarAligner::_optimal_path(const std::vector<triple> &identical_substrings)
        -> std::vector<triple> 
{
    std::vector<triple> optimal_identical_substrings;
    if (identical_substrings.empty()) return optimal_identical_substrings;

    const size_t pair_num = identical_substrings.size();
    utils::AdjacencyList graph(pair_num + 1);

    for (size_t i = 0; i != pair_num; ++i)
    for (size_t j = 0; j != pair_num; ++j)
        if (i != j && identical_substrings[i][0] + identical_substrings[i][2] < identical_substrings[j][0] + identical_substrings[j][2]
                   && identical_substrings[i][1] + identical_substrings[i][2] < identical_substrings[j][1] + identical_substrings[j][2])
        {
            const int possible_overlap = std::max(
                    static_cast<int>(identical_substrings[i][0] + identical_substrings[i][2]) - static_cast<int>(identical_substrings[j][0]), 
                    static_cast<int>(identical_substrings[i][1] + identical_substrings[i][2]) - static_cast<int>(identical_substrings[j][1]));

            unsigned weight = identical_substrings[j][2];
            if (possible_overlap > 0) weight -= possible_overlap;
            graph.add_edge(i + 1, j + 1, weight);
        }

    for (size_t i = 0; i != pair_num; ++i)
        graph.add_edge(0, i + 1, identical_substrings[i][2]);

    const auto optimal_path = graph.get_longest_path();

    optimal_identical_substrings.reserve(optimal_path.size());
    optimal_identical_substrings.push_back(triple( { identical_substrings[optimal_path[0] - 1][0],
                                                     identical_substrings[optimal_path[0] - 1][1],
                                                     identical_substrings[optimal_path[0] - 1][2] }));

    for (size_t i = 0; i < optimal_path.size() - 1; ++i)
    {
        size_t new_len = graph.get_weight(optimal_path[i], optimal_path[i + 1]);
        size_t old_len = identical_substrings[optimal_path[i + 1] - 1][2];
        int difference = static_cast<int>(old_len) - static_cast<int>(new_len);

        size_t lhs_first = identical_substrings[optimal_path[i + 1] - 1][0];
        size_t rhs_first = identical_substrings[optimal_path[i + 1] - 1][1];
        if (difference > 0)
        { lhs_first += difference; rhs_first += difference; }

        optimal_identical_substrings.push_back(triple({ lhs_first, rhs_first, new_len }));
    }

    return optimal_identical_substrings;
}

void star_alignment::StarAligner::_append(const std::vector<size_t> &src_gaps, std::vector<utils::Insertion> &des_gaps, size_t start)
{
    for (size_t i = 0; i != src_gaps.size(); ++i)
        if (src_gaps[i])
        {
            if (des_gaps.size() && des_gaps.back().index == start + i)
                des_gaps.back().number += src_gaps[i];
            else
                des_gaps.push_back(utils::Insertion({ start + i, src_gaps[i] }));
        }
}

auto star_alignment::StarAligner::_merge_results(const std::vector<std::array<std::vector<utils::Insertion>, 2>> &pairwise_gaps) const
        -> std::vector<std::vector<utils::Insertion>>
{
    std::vector<utils::Insertion> final_centre_gaps;
    for (size_t i = 0; i != _row; ++i)
    {
        const auto &curr_centre_gaps = pairwise_gaps[i][0];
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

    std::vector<std::vector<utils::Insertion>> final_sequence_gaps;
    final_sequence_gaps.reserve(_row);
    for (size_t i = 0; i != _row; ++i)
    {
        const auto &curr_centre_gaps = pairwise_gaps[i][0];
        const auto &curr_sequence_gaps = pairwise_gaps[i][1];

        std::vector<utils::Insertion> centre_addition;
        centre_addition.reserve(final_centre_gaps.size());
        utils::Insertion::minus(final_centre_gaps.cbegin(), final_centre_gaps.cend(),
                            curr_centre_gaps.cbegin(),  curr_centre_gaps.cend(),
                            std::back_inserter(centre_addition));

        std::vector<utils::Insertion> sequence_addition;
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
                sequence_addition.push_back(utils::Insertion({ sequence_index, curr_addition.number }));
        }

        std::vector<utils::Insertion> indels_of_current_sequence;
        indels_of_current_sequence.reserve(curr_sequence_gaps.size() + sequence_addition.size());
        utils::Insertion::plus(curr_sequence_gaps.cbegin(), curr_sequence_gaps.cend(),
                           sequence_addition.cbegin(), sequence_addition.cend(),
                           std::back_inserter(indels_of_current_sequence));
        final_sequence_gaps.push_back(indels_of_current_sequence);
    }

    return final_sequence_gaps;
}

auto star_alignment::StarAligner::_insert_gaps(const std::vector<std::vector<utils::Insertion>> &gaps) const
        -> std::vector<sequence_type>
{
    size_t length = _lengths[0];
    for (const auto gap : gaps[0]) length += gap.number;

    std::vector<sequence_type> aligned;
    aligned.reserve(_row);

    for (size_t i = 0; i != _row; ++i)
    {
        aligned.emplace_back(length);
        utils::Insertion::insert_gaps(_sequences[i].cbegin(), _sequences[i].cend(),
                gaps[i].cbegin(), gaps[i].cend(), aligned.back().begin(), nucleic_acid_pseudo::GAP);
    }
    return aligned;
}
