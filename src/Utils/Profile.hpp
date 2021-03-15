#pragma once

#include <iterator>
#include <vector>

namespace utils
{

    template <unsigned char max_element>
    class Profile
    {
    public:
        using self_type = Profile<max_element>;

        template<typename InputIterator>
        Profile(InputIterator first, InputIterator last, size_t index) noexcept
        {
            constexpr size_t len = max_element + 1;
            
            unsigned frequencies[len];
            memset(frequencies, 0, sizeof(frequencies));

            size_t n = std::distance(first, last);
            for (size_t i = 0; i != n; ++i, ++first) ++frequencies[(*first)[index]];

            for (size_t i = 0; i != len; ++i)
                _relative_frequencies[i] = static_cast<double>(frequencies[i]) / n;

            _master = 0;
            for (unsigned char i = 1; i != len; ++i)
                if (_relative_frequencies[i] > _relative_frequencies[_master])
                    _master = i;
        }

        Profile(const self_type& rhs) noexcept
        {
            _copy_from(rhs);
        }

        Profile(self_type&& rhs) noexcept
        {
            _copy_from(rhs);
        }

        self_type& operator=(const self_type& rhs) noexcept
        {
            _copy_from(rhs);
            return *this;
        }

        self_type& operator=(self_type&& rhs) noexcept
        {
            _copy_from(rhs);
            return *this;
        }

        Profile() = delete;

        ~Profile() = default;

        operator size_t() const noexcept
        {
            return _master;
        }

    private:
        inline void _copy_from(const self_type& rhs) noexcept
        {
            _master = rhs._master;
            memcpy(_relative_frequencies, rhs._relative_frequencies, sizeof(_relative_frequencies));
        }

        double _relative_frequencies[max_element + 1];
        unsigned char _master;
    };

    template <typename CharMatrix>
    std::vector<unsigned char> abstract(const CharMatrix& matrix, size_t row, size_t clm, size_t max_element);

}

template <typename CharMatrix>
std::vector<unsigned char> utils::abstract(const CharMatrix& matrix, size_t row, size_t clm, size_t max_element)
{
    std::vector<unsigned char> result(clm);

    size_t count_size = sizeof(size_t) * (max_element + 1);
    size_t* count = new size_t[max_element + 1];
    memset(count, 0, count_size);

    for (size_t j = 0; j != clm; ++j)
    {
        for (size_t i = 0; i != row; ++i) ++count[matrix[i][j]];

        unsigned char master = 0;
        for (size_t k = 1; k != max_element + 1; ++k)
            if (count[k] > count[master]) master = k;

        result[j] = master;
        memset(count, 0, count_size);
    }

    delete[] count;
    return result;
}
