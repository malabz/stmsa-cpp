#pragma once

#include "../SuffixTree/SuffixTreeBuilder.hpp"

namespace suffixtree
{

    class SuffixTreeBuilderTester
    {
    public:
        static void test();

        static void test_profile(const char* lhs_file_path, const char* rhs_file_path);

        static void test_substring(const char* lhs_file_path, const char* rhs_file_path);

        template<typename T>
        static void print(const std::vector<std::vector<T>>& vt);
    };

}

