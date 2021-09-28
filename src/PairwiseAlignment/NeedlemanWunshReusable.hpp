#pragma once

#include "PairwiseAlignment.hpp"

#include <vector>
#include <algorithm>

namespace pairwise_alignment
{

    template <typename RandomAccessIterator1, typename RandomAccessIterator2,
              typename ScoringMatrixType>
    class NeedlemanWunshReusable
    {
    public:
        using gap_vector_type = std::vector<size_t>;

        NeedlemanWunshReusable(const ScoringMatrixType &scoring_matrix, int gap_open, int gap_extention)
                : _scoring_matrix(scoring_matrix)
                , _op(gap_open)
                , _ex(gap_extention)
                , _lhs_capacity(0)
                , _rhs_capacity(0)
        {}

        ~NeedlemanWunshReusable()
        {
            _clear();
        }

        std::tuple<gap_vector_type, gap_vector_type>
        operator ()(RandomAccessIterator1 lhs_first, RandomAccessIterator1 lhs_last,
              RandomAccessIterator2 rhs_first, RandomAccessIterator2 rhs_last, unsigned flag)
        {
            _inisialise(lhs_first, lhs_last, rhs_first, rhs_last, flag);
            _do_dp();
            _trace_back();
            return std::make_tuple(std::move(_lhs_gaps), std::move(_rhs_gaps));
        }

    private:
        void _clear()
        {
            if (_lhs_capacity)
            {
                for (size_t i = 0; i != NUM; ++i)
                {
                    for (size_t j = 0; j != _lhs_capacity; ++j)
                    {
                        delete[] _dp_matrix[i][j];
                        delete[] _pa_matrix[i][j];
                    }

                    delete[] _dp_matrix[i];
                    delete[] _pa_matrix[i];
                }

                delete[] _dp_matrix;
                delete[] _pa_matrix;
            }
        }

        void _inisialise(RandomAccessIterator1 lhs_first, RandomAccessIterator1 lhs_last,
                         RandomAccessIterator2 rhs_first, RandomAccessIterator2 rhs_last, unsigned flag)
        {
            _lhs_first = lhs_first; _lhs_last  = lhs_last; _lhs_len = std::distance(lhs_first, lhs_last);
            _rhs_first = rhs_first; _rhs_last  = rhs_last; _rhs_len = std::distance(rhs_first, rhs_last);

            bool _l_ending = flag & LEFT_ENDING;
            bool _r_ending = flag & RIGHT_ENDING;

            _l_op = _l_ending ? 0 : _op;
            _r_op = _r_ending ? 0 : _op;

            _l_ex = _l_ending ? 0 : _ex;
            _r_ex = _r_ending ? 0 : _ex;

            if (_lhs_capacity > _lhs_len && _rhs_capacity > _rhs_len)
            {
                for (size_t i = 0; i != NUM; ++i)
                    for (size_t j = 0; j <= _lhs_len; ++j)
                        memset(_dp_matrix[i][j], 0, _rhs_len + 1);
            }
            else
            {
                _clear();
                _lhs_gaps.clear();
                _rhs_gaps.clear();

                _dp_matrix = new int**[NUM];
                _pa_matrix = new unsigned char**[NUM];
                for (size_t i = 0; i != NUM; ++i)
                {
                    _dp_matrix[i] = new int*[_lhs_len + 1];
                    _pa_matrix[i] = new unsigned char*[_lhs_len + 1];

                    for (size_t j = 0; j <= _lhs_len; ++j)
                    {
                        _dp_matrix[i][j] = new int[_rhs_len + 1]();
                        _pa_matrix[i][j] = new unsigned char[_rhs_len + 1];
                    }
                }

                _lhs_capacity = _lhs_len + 1;
                _rhs_capacity = _rhs_len + 1;
            }

            _lhs_gaps.resize(_lhs_len + 1); std::fill(_lhs_gaps.begin(), _lhs_gaps.end(), 0);
            _rhs_gaps.resize(_rhs_len + 1); std::fill(_rhs_gaps.begin(), _rhs_gaps.end(), 0);

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
            auto lhs_index = _lhs_len;
            auto rhs_index = _rhs_len;
            unsigned char curr_path = _index_of_max({ _dp_matrix[VER][lhs_index][rhs_index],
                                                      _dp_matrix[HOR][lhs_index][rhs_index],
                                                      _dp_matrix[DIA][lhs_index][rhs_index] });

            while (lhs_index > 0 && rhs_index > 0)
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

            for (; lhs_index > 0; --lhs_index) ++_rhs_gaps[rhs_index];
            for (; rhs_index > 0; --rhs_index) ++_lhs_gaps[lhs_index];
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
        const ScoringMatrixType &_scoring_matrix;

        int _op, _l_op, _r_op;
        int _ex, _l_ex, _r_ex;

        // sequences
        RandomAccessIterator1 _lhs_first, _lhs_last;
        RandomAccessIterator2 _rhs_first, _rhs_last;

        size_t _lhs_len, _lhs_capacity;
        size_t _rhs_len, _rhs_capacity;

        // dynamic programing
        int           ***_dp_matrix;
        unsigned char ***_pa_matrix; // path

        // result
        gap_vector_type _lhs_gaps;
        gap_vector_type _rhs_gaps;

    };

}
