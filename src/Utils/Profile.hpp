#pragma once

#include "Pseudo.hpp"

#include <iterator>
#include <vector>

namespace utils
{

    // template <size_t T>
    class Profile
    {
    public:
        friend bool operator==(const Profile& lhs, const Profile& rhs) noexcept;
        friend bool operator!=(const Profile& lhs, const Profile& rhs) noexcept;

        template<typename InputIterator>
        Profile(InputIterator first, InputIterator last, size_t index) noexcept;

        Profile(const Profile& rhs) noexcept;

        Profile(const std::initializer_list<double>& il, unsigned char master) noexcept;

        Profile() = default;

        Profile& operator=(const Profile& rhs) noexcept;

        operator unsigned() const noexcept;

        static const Profile end_mark;

    private:
        static const double _THRESHOLD;

        double _relative_frequencies[pseudo::NUMBER];
        unsigned char _master;
    };

    template <typename CharMatrix>
    std::vector<unsigned char> abstract(const CharMatrix& matrix, size_t row, size_t clm, size_t max_element);

}

template <typename InputIterator>
utils::Profile::Profile(InputIterator first, InputIterator last, size_t index) noexcept
{
    unsigned frequencies[pseudo::NUMBER];
    memset(frequencies, 0, sizeof(frequencies));

    size_t n = std::distance(first, last);
    for (size_t i = 0; i != n; ++i, ++first) ++frequencies[(*first)[index]];

    for (size_t i = 0; i != pseudo::NUMBER; ++i)
        _relative_frequencies[i] = static_cast<double>(frequencies[i]) / n;

    _master = 0;
    for (unsigned char i = 1; i != pseudo::NUMBER; ++i)
        if (_relative_frequencies[i] > _relative_frequencies[_master])
            _master = i;
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
        for (size_t i = 0; i != clm; ++i) ++count[matrix[i][j]];

        unsigned char master = 0;
        for (size_t k = 1; k != max_element + 1; ++k)
            if (count[k] > count[master]) master = k;

        result[i] = master;
        memset(count, 0, count_size);
    }

    delete[] max_element;
    return result;
}
