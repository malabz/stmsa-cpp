#include "../PairwiseAlignment/PairwiseAlignment.hpp"
#include "PairwiseAlignmentTester.hpp"
#include "../Utils/Pseudo.hpp"

#include <vector>


void pairwise_alignment::needleman_wunsh_test()
{
    std::vector<unsigned char> lhs({ 2, 3, 2, 3, 4 });
    std::vector<unsigned char> rhs({ 2, 3, 4 });

    auto [lhs_gaps, rhs_gaps] = needleman_wunsh(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend(), default_scoring_matrix);
    auto lhs_aligned = insert_gaps(lhs, lhs_gaps);
    auto rhs_aligned = insert_gaps(rhs, rhs_gaps);

    print_sequence(std::cout, lhs_aligned.cbegin(), lhs_aligned.cend()); std::cout << '\n';
    print_sequence(std::cout, rhs_aligned.cbegin(), rhs_aligned.cend()); std::cout << '\n';
}
