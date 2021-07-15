#pragma once

#include "Fasta.hpp"
#include "Pseudo.hpp"
#include "Insertion.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iostream>

namespace utils
{
    std::string remove_white_spaces(const std::string &str);

    unsigned char to_pseudo(char c);

    std::vector<unsigned char> to_pseudo(const std::string &str);
    std::string from_pseudo(const std::vector<unsigned char> &pseu);

    void print_duration(std::chrono::system_clock::time_point time_point);
    void print_duration(std::chrono::system_clock::time_point time_point, const std::string &info);

    template<typename InputIterator, typename OutputIterator>
    void transform_to_pseudo(InputIterator src_first, InputIterator src_last, OutputIterator des)
    {
        std::vector<unsigned char> (*op)(const std::string &) = &to_pseudo;
        std::transform(src_first, src_last, des, op);
    }

    template<typename InputIterator, typename OutputIterator>
    void transform_from_pseudo(InputIterator src_first, InputIterator src_last, OutputIterator des)
    {
        std::string (*op)(const std::vector<unsigned char> &) = &from_pseudo;
        std::transform(src_first, src_last, des, op);
    }

    template<typename InputIterator>
    InputIterator iter_of_max(InputIterator first, InputIterator last)
    {
        auto result = first;

        for (; first != last; ++first) if (*result < *first) result = first;
        return result;
    }

    std::vector<std::vector<unsigned char>> read_to_pseudo(std::istream &is);

    void insert_and_write(std::ostream &os, std::istream &is, const std::vector<std::vector<Insertion>> &insertions);

}
