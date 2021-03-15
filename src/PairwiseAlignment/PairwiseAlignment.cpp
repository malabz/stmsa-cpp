#include "PairwiseAlignment.hpp"
#include "../Utils/Pseudo.hpp"

#include <numeric>

std::vector<unsigned char>
pairwise_alignment::insert_gaps(const std::vector<unsigned char>& sequence, const std::vector<size_t>& gaps)
{
    size_t total = std::accumulate(gaps.cbegin(), gaps.cend(), static_cast<size_t>(0));
    std::vector<unsigned char> result(total + sequence.size());
    std::cout << result.size() << ' ' << total << ' ' << sequence.size() << ' ' << gaps.size() << std::endl;

    size_t result_index = 0;
    for (size_t i = 0; i != sequence.size(); ++i)
    {
        for (size_t j = 0; j != gaps[i]; ++j) result[result_index++] = nucleic_acid_pseudo::GAP;
        result[result_index++] = sequence[i];
    }
    for (size_t j = 0; j != gaps.back(); ++j) result[result_index++] = nucleic_acid_pseudo::GAP;

    return result;
}
