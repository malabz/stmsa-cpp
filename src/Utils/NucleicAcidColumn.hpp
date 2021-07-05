#pragma once

#include "Pseudo.hpp"

#include <iterator>
#include <vector>
#include <algorithm>

namespace utils
{

    class NucleicAcidColumn
    {
    public:
        template <typename MatrixType>
        NucleicAcidColumn(const MatrixType& matrix, size_t num_row, size_t col)
        {
            unsigned int cnt[_number] = {};
            for (size_t i = 0; i != num_row; ++i)
                ++cnt[matrix[i][col]];

            double relative_frequencies[_number] = {};
            for (size_t i = 0; i != _number; ++i)
                relative_frequencies[i] = static_cast<double>(cnt[i]) / num_row;

            _split_unknown(relative_frequencies);

            _feature[nucleic_acid_pseudo::GAP] = nucleic_acid_pseudo::GAP;
            _feature[nucleic_acid_pseudo::C  ] = nucleic_acid_pseudo::C;
            _feature[nucleic_acid_pseudo::G  ] = nucleic_acid_pseudo::G;
            _feature[nucleic_acid_pseudo::A  ] = nucleic_acid_pseudo::A;
            _feature[nucleic_acid_pseudo::T  ] = nucleic_acid_pseudo::T;
            std::sort(_feature, _feature + (_number - 1),
                    [&relative_frequencies](unsigned char lhs, unsigned char rhs)
                    { return relative_frequencies[lhs] > relative_frequencies[rhs]; });

            double sum = 0;
            for (size_t i = 0; ; ++i)
            {
                sum += relative_frequencies[_feature[i]];
                if (sum > _thresholds[i])
                {
                    _available_count = i + 1;
                    break;
                }
            }
        }

        operator unsigned char() const noexcept;

    private:
        static void _split_unknown(double* relative_frequencies) noexcept;

        static constexpr size_t _number = nucleic_acid_pseudo::NUMBER;
        static constexpr double _thresholds[]{ 0.9, 0.99, 0.999, 0 };

        size_t _available_count;
        unsigned char _feature[_number - 1];
    };

}
