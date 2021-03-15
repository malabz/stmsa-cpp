#pragma once

#include "NeedlemanWunsh.hpp"

#include <vector>
#include <iostream>

namespace pairwise_alignment
{

    std::vector<unsigned char>
    insert_gaps(const std::vector<unsigned char>& sequence, const std::vector<size_t>& gaps);

    template <typename RandomAccessIterator1,
              typename RandomAccessIterator2,
              typename ScoringMatrixType>
    auto needleman_wunsh(RandomAccessIterator1 lhs_first, RandomAccessIterator1 lhs_last,
                         RandomAccessIterator2 rhs_first, RandomAccessIterator2 rhs_last,
                         const ScoringMatrixType& scoring_matrix)
    {
        return NeedlemanWunsh<RandomAccessIterator1, RandomAccessIterator2, ScoringMatrixType>
                (lhs_first, lhs_last, rhs_first, rhs_last, scoring_matrix)._align();
    }

    template <typename InputIterator>
    inline void print_sequence(std::ostream& is, InputIterator first, InputIterator last)
    {
        for (; first != last; ++first) is << static_cast<unsigned int>(*first);
    }

}
