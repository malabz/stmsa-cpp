#pragma once

#include "PairwiseAlignment.hpp"

#include <vector>
#include <algorithm>

namespace pairwise_alignment
{

    template <typename RandomAccessIterator1, typename RandomAccessIterator2,
              typename ScoringMatrixType>
    class NeedlemanWunsh
    {
    public:
        template <typename RandomAccessIterator1, typename RandomAccessIterator2,
                  typename ScoringMatrixType>
        friend
        std::tuple<std::vector<size_t>, std::vector<size_t>>
        needleman_wunsh(RandomAccessIterator1 lhs_first, RandomAccessIterator1 lhs_last,
                        RandomAccessIterator2 rhs_first, RandomAccessIterator2 rhs_last,
                        unsigned flag, const ScoringMatrixType& scoring_matrix,
                        int gap_open, int gap_extention);

    private:
        using gap_vector_type = std::vector<size_t>;

        // n-dimension int matrix
        using d1im = std::vector<int>;
        using d2im = std::vector<d1im>;
        using d3im = std::vector<d2im>;

        // n-dimension path matrix
        using d1pm = std::vector<unsigned char>;
        using d2pm = std::vector<d1pm>;
        using d3pm = std::vector<d2pm>;

        NeedlemanWunsh(RandomAccessIterator1 lhs_first, RandomAccessIterator1 lhs_last,
                       RandomAccessIterator2 rhs_first, RandomAccessIterator2 rhs_last,
                       unsigned flag, const ScoringMatrixType& scoring_matrix,
                       int gap_open, int gap_extention):

                _lhs_first(lhs_first), _lhs_last(lhs_last),
                _rhs_first(rhs_first), _rhs_last(rhs_last),

                _lhs_len(std::distance(lhs_first, lhs_last)),
                _rhs_len(std::distance(rhs_first, rhs_last)),

                _scoring_matrix(scoring_matrix),
                _l_ending(flag & LEFT_ENDING),
                _r_ending(flag & RIGHT_ENDING),

                _op(gap_open),
                _l_op(_l_ending ? 0 : gap_open),
                _r_op(_r_ending ? 0 : gap_open),

                _ex(gap_extention),
                _l_ex(_l_ending ? 0 : gap_extention),
                _r_ex(_r_ending ? 0 : gap_extention),

                _dp_matrix(NUM, d2im(_lhs_len + 1, d1im(_rhs_len + 1, 0))),
                _pa_matrix(NUM, d2pm(_lhs_len + 1, d1pm(_rhs_len + 1))),

                _lhs_gaps(_lhs_len + 1, 0),
                _rhs_gaps(_rhs_len + 1, 0)
        {}

        std::tuple<gap_vector_type, gap_vector_type> _align()
        {
            _inisialise();
            _do_dp();
            _trace_back();
            return std::make_tuple(std::move(_lhs_gaps), std::move(_rhs_gaps));
        }

        void _inisialise()
        {
            // _dp_matrix

            for (size_t i = 0; i <= _lhs_len; ++i)
            {
                auto score = _l_op + i * _l_ex;

                _dp_matrix[VER][i][0] = score;
                _dp_matrix[HOR][i][0] = NEGATIVE_INFINITY;
                _dp_matrix[DIA][i][0] = score;
            }

            for (size_t j = 0; j <= _rhs_len; ++j)
            {
                auto score = _l_op + j * _l_ex;

                _dp_matrix[VER][0][j] = NEGATIVE_INFINITY;
                _dp_matrix[HOR][0][j] = score;
                _dp_matrix[DIA][0][j] = score;
            }

            _dp_matrix[VER][0][0] = 0;
            _dp_matrix[HOR][0][0] = 0;
            _dp_matrix[DIA][0][0] = 0;

            // _pa_matrix

            for (size_t i = 0; i <= _lhs_len; ++i) { _pa_matrix[VER][i][0] = VER; _pa_matrix[DIA][i][0] = VER; }
            for (size_t j = 0; j <= _rhs_len; ++j) { _pa_matrix[HOR][0][j] = HOR; _pa_matrix[DIA][0][j] = HOR; }
        }

        void _do_dp()
        {
            for (size_t i = 1; i <= _lhs_len; ++i)
            for (size_t j = 1; j <= _rhs_len; ++j)
            {
                const bool is_right = j == _rhs_len;
                const int op = is_right ? _r_op : _op;
                const int ex = is_right ? _r_ex : _ex;

                int scores[NUM];
                unsigned char curr_path;

                scores[VER] = _dp_matrix[VER][i - 1][j];
                // scores[HOR] = _dp_matrix[HOR][i - 1][j] + op; // 
                scores[DIA] = _dp_matrix[DIA][i - 1][j] + op;
                curr_path = _index_of_max({ scores[VER], NEGATIVE_INFINITY, scores[DIA] });
                _dp_matrix[VER][i][j] = scores[curr_path] + ex;
                _pa_matrix[VER][i][j] = curr_path;

                // scores[VER] = _dp_matrix[VER][i][j - 1] + op; // 
                scores[HOR] = _dp_matrix[HOR][i][j - 1];
                scores[DIA] = _dp_matrix[DIA][i][j - 1] + op;
                curr_path = _index_of_max({ NEGATIVE_INFINITY, scores[HOR], scores[DIA] });
                _dp_matrix[HOR][i][j] = scores[curr_path] + ex;
                _pa_matrix[HOR][i][j] = curr_path;

                scores[VER] = _dp_matrix[VER][i - 1][j - 1];
                scores[HOR] = _dp_matrix[HOR][i - 1][j - 1];
                scores[DIA] = _dp_matrix[DIA][i - 1][j - 1];
                curr_path = _index_of_max({ scores[VER], scores[HOR], scores[DIA] });
                _dp_matrix[DIA][i][j] = scores[curr_path] + _scoring_matrix[_lhs_first[i - 1]][_rhs_first[j - 1]];
                _pa_matrix[DIA][i][j] = curr_path;
            }
        }

        void _trace_back()
        {
            int lhs_index = _lhs_len;
            int rhs_index = _rhs_len;
            unsigned char curr_path = _index_of_max({ _dp_matrix[VER][lhs_index][rhs_index],
                                                      _dp_matrix[HOR][lhs_index][rhs_index],
                                                      _dp_matrix[DIA][lhs_index][rhs_index] });

            while (lhs_index > 0 || rhs_index > 0)
                switch (curr_path)
                {
                case VER:
                    curr_path = _pa_matrix[curr_path][lhs_index--][rhs_index];
                    ++_rhs_gaps[rhs_index];
                    break;

                case HOR:
                    curr_path = _pa_matrix[curr_path][lhs_index][rhs_index--];
                    ++_lhs_gaps[lhs_index];
                    break;

                default:
                    curr_path = _pa_matrix[curr_path][lhs_index--][rhs_index--];
                    break;
                }
        }

        template <typename T>
        static unsigned char _index_of_max(std::initializer_list<T> il)
        {
            auto begin = il.begin();
            auto iter = il.begin();

            for (++begin; begin != il.end(); ++begin)
                if (*iter < *begin) iter = begin;
            return std::distance(il.begin(), iter);
        }

        // path
        static constexpr unsigned char VER = 0; // vertical
        static constexpr unsigned char HOR = 1; // horizontal
        static constexpr unsigned char DIA = 2; // diagonal
        static constexpr unsigned char NUM = 3; // number

        // util
        static constexpr int NEGATIVE_INFINITY = std::numeric_limits<int>::min() / 2;

        // scoring
        const ScoringMatrixType& _scoring_matrix;
        const bool               _l_ending, _r_ending;
        const int                _op, _l_op, _r_op;
        const int                _ex, _l_ex, _r_ex;

        // sequences
        const RandomAccessIterator1 _lhs_first, _lhs_last;
        const RandomAccessIterator2 _rhs_first, _rhs_last;
        const size_t                _lhs_len,   _rhs_len;

        // dynamic programing
        d3im _dp_matrix;
        d3pm _pa_matrix; // path

        // result
        gap_vector_type _lhs_gaps;
        gap_vector_type _rhs_gaps;

    };

    template <typename RandomAccessIterator1, typename RandomAccessIterator2,
              typename ScoringMatrixType = decltype(default_scoring_matrix)>
    std::tuple<std::vector<size_t>, std::vector<size_t>>
    needleman_wunsh(RandomAccessIterator1 lhs_first, RandomAccessIterator1 lhs_last,
                    RandomAccessIterator2 rhs_first, RandomAccessIterator2 rhs_last,
                    unsigned flag = DEFAULT,
                    const ScoringMatrixType& scoring_matrix = default_scoring_matrix,
                    int gap_open = default_gap_open,
                    int gap_extention = default_gap_extention)
    {
        return NeedlemanWunsh<RandomAccessIterator1, RandomAccessIterator2, ScoringMatrixType>
                (lhs_first, lhs_last, rhs_first, rhs_last, flag, scoring_matrix, gap_open, gap_extention)._align();
    }

}
