#include "../PairwiseAlignment/PairwiseAlignment.hpp"
#include "PairwiseAlignmentTester.hpp"
#include "../Utils/Pseudo.hpp"

#include <vector>


void pairwise_alignment::needleman_wunsh_test()
{
    static int scoring_matrix[nucleic_acid_pseudo::NUMBER][nucleic_acid_pseudo::NUMBER]
    {
        { 0, 0, 0, 0, 0, 0 },
        { 0, 7, -3, -3, -3, 7 },
        { 0, -3, 7, -3, -3, 7 },
        { 0, -3, -3, 7, -3, 7 },
        { 0, -3, -3, -3, 7, 7 },
        { 0, 7, 7, 7, 7, 7 },
    };

    std::vector<unsigned char> lhs({ 2, 3, 2, 3, 4 });
    std::vector<unsigned char> rhs({ 2, 3, 4 });

    auto result = needleman_wunsh(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend(), scoring_matrix);
    auto lhs_result = insert_gaps(lhs, result.first);
    auto rhs_result = insert_gaps(rhs, result.second);

    print_sequence(std::cout, lhs_result.cbegin(), lhs_result.cend()); std::cout << std::endl;
    print_sequence(std::cout, rhs_result.cbegin(), rhs_result.cend()); std::cout << std::endl;
}
