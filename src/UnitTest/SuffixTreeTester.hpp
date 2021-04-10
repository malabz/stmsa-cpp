#pragma once

#include "../SuffixTree/SuffixTree.hpp"

#include <iostream>
#include <iomanip>

namespace suffixtree
{

    static void test();

    static void test_profile(const char* lhs_file_path, const char* rhs_file_path);

    template<typename MatrixType>
    void print_matrix(const MatrixType& matrix, size_t row, size_t clm)
    {
        for (size_t i = 0; i != row; ++i)
        {
            for (size_t j = 0; j != clm; ++j)
                std::cout << std::right << std::setw(8) << matrix[i][j];
            std::cout << '\n';
        }
    }

}

