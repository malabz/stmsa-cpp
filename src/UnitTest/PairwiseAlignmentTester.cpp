#include "../PairwiseAlignment/NeedlemanWunsh.hpp"
#include "../PairwiseAlignment/NeedlemanWunshReusable.hpp"
#include "PairwiseAlignmentTester.hpp"
#include "../Utils/Pseudo.hpp"

#include <iostream>

void pairwise_alignment::needleman_wunsh_test()
{
    needleman_wunsh_test(std::vector<unsigned char>({ 2, 3, 2, 3 }),
                         std::vector<unsigned char>({ 2, 3 }), DEFAULT);

    needleman_wunsh_test(std::vector<unsigned char>({ 2, 3, 2, 3 }),
                         std::vector<unsigned char>({ 2, 3 }), LEFT_ENDING);

    needleman_wunsh_test(std::vector<unsigned char>({ 2, 3, 2, 3 }),
                         std::vector<unsigned char>({ 2, 3 }), RIGHT_ENDING);
}

void pairwise_alignment::needleman_wunsh_test(const std::vector<unsigned char>& lhs,
                                              const std::vector<unsigned char>& rhs, unsigned flag)
{
    auto [lhs_gaps, rhs_gaps] = needleman_wunsh(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend(), flag);

    // using iter_type = std::vector<unsigned char>::const_iterator;
    // static NeedlemanWunshReusable<iter_type, iter_type, decltype(default_scoring_matrix)> nw(default_scoring_matrix, default_gap_open, default_gap_extention);
    // auto [lhs_gaps, rhs_gaps] = nw(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend(), flag);

    auto lhs_aligned = insert_gaps(lhs, lhs_gaps);
    auto rhs_aligned = insert_gaps(rhs, rhs_gaps);

    std::copy(lhs_aligned.cbegin(), lhs_aligned.cend(), std::ostream_iterator<unsigned int>(std::cout)); std::cout << '\n';
    std::copy(rhs_aligned.cbegin(), rhs_aligned.cend(), std::ostream_iterator<unsigned int>(std::cout)); std::cout << '\n';
}
