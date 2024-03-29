#pragma once

#include "../SuffixTree/SuffixTree.hpp"
#include "../Utils/Utils.hpp"

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
        static std::vector<sequence_type> align(const std::vector<sequence_type> &sequences);
        static std::vector<std::vector<utils::Insertion>> get_gaps(const std::vector<sequence_type> &sequences);

        static std::vector<triple> _optimal_path(const std::vector<triple> &common_substrings);

    private:
        StarAligner(const std::vector<sequence_type> &sequences);

        std::vector<sequence_type> _align() const;
        std::vector<std::vector<utils::Insertion>> _get_gaps() const;

        std::vector<size_t> _set_lengths() const;
        sequence_type _set_centre() const;

        // main steps
        std::vector<std::array<std::vector<utils::Insertion>, 2>> _pairwise_align() const;
        std::vector<std::vector<utils::Insertion>> _merge_results(const std::vector<std::array<std::vector<utils::Insertion>, 2>> &pairwise_gaps) const;
        std::vector<sequence_type> _insert_gaps(const std::vector<std::vector<utils::Insertion>> &gaps) const;

        // support
        static void _append(const std::vector<size_t> &src_gaps, std::vector<utils::Insertion> &des_gaps, size_t start);

        const std::vector<sequence_type> &_sequences;
        const size_t _row;
        const std::vector<size_t> _lengths;

        const sequence_type _centre;
        const size_t _centre_len;

    };

}
