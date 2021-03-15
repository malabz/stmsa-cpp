#pragma once

#include "../SuffixTree/SuffixTree.hpp"

namespace suffixtree
{

    class SuffixTreeTester
    {
    public:

        static void test();

        static void test_profile(const char* lhs_file_path, const char* rhs_file_path);

        template<typename T>
        static void print(const std::vector<std::vector<T>>& vt);
    };

}

