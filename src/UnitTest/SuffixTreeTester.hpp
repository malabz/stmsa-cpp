#pragma once

#include "../SuffixTree/SuffixTree.hpp"

namespace suffixtree
{

    class SuffixTreeTester
    {
    public:

        static void test();

        static void test_profile(const char* lhs_file_path, const char* rhs_file_path);

        template<typename MatrixType>
        static void print_matrix(const MatrixType& matrix, size_t row, size_t clm);
    };

}

