#pragma once

#include "../Utils/Pseudo.hpp"

#include <vector>
#include <iostream>

namespace pairwise_alignment
{
    constexpr unsigned DEFAULT      = 0x0;
    constexpr unsigned LEFT_ENDING  = 0x1;
    constexpr unsigned RIGHT_ENDING = 0x2;

    constexpr int default_scoring_matrix[nucleic_acid_pseudo::NUMBER][nucleic_acid_pseudo::NUMBER]
    { // '-' dismissed
        { 0,  0,  0,  0,  0, 0 },
        { 0,  7, -3, -3, -3, 7 },
        { 0, -3,  7, -3, -3, 7 },
        { 0, -3, -3,  7, -3, 7 },
        { 0, -3, -3, -3,  7, 7 },
        { 0,  7,  7,  7,  7, 7 },
    };

    constexpr int default_gap_open      = -11;
    constexpr int default_gap_extention = -2;

    std::vector<unsigned char>
    insert_gaps(const std::vector<unsigned char>& sequence, const std::vector<size_t>& gaps);

}
