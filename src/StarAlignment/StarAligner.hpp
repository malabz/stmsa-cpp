#pragma once

#include "../SuffixTree/SuffixTree.hpp"

#include <vector>
#include <array>

namespace star_alignment
{

    class StarAligner
    {
    private:
        using triple = std::array<size_t, 3>;
        using quadra = std::array<size_t, 4>;
        using sequence_type = std::vector<unsigned char>;

    public:
        static std::vector<sequence_type> align(const std::vector<sequence_type>& sequences);
        static std::vector<triple> _optimal_path(const std::vector<triple>& identical_substrings);

    private:
        struct indel
        {
            size_t index;
            size_t number;

            bool operator==(const indel& rhs)  { return index == rhs.index && number == rhs.number; }
        };

        StarAligner(const std::vector<sequence_type>& sequences);

        std::vector<sequence_type> _align() const;

        std::vector<size_t> _set_lengths() const;
        sequence_type _set_centre() const;

        // main steps
        std::vector<std::array<std::vector<indel>, 2>> _pairwise_align() const;
        std::vector<std::vector<indel>> _merge_results(const std::vector<std::array<std::vector<indel>, 2>>& pairwise_gaps) const;
        std::vector<sequence_type> _insert_gaps(const std::vector<std::vector<indel>>& gaps) const;

        // support
        static void _converse(const std::vector<size_t>& src_gaps, std::vector<indel>& des_gaps, size_t start);
        static std::vector<indel> _add(const std::vector<indel>& lhs, const std::vector<indel>& rhs);
        static std::vector<indel> _minus(const std::vector<indel>& lhs, const std::vector<indel>& rhs);
        static sequence_type _insert_gaps(const sequence_type& sequence, const std::vector<indel>& gaps, size_t length);

        const std::vector<sequence_type>& _sequences;
        const size_t _row;
        const std::vector<size_t> _lengths;

        const sequence_type _centre;
        const size_t _centre_len;

    };

}
