#pragma once

#include "Pseudo.hpp"

#include <iterator>

namespace utils
{

    class Profile
    {
    public:
        friend bool operator==(const Profile& lhs, const Profile& rhs) noexcept;
        friend bool operator!=(const Profile& lhs, const Profile& rhs) noexcept;

        template<typename InputIterator>
        Profile(InputIterator first, InputIterator last, size_t index) noexcept;

        Profile(const Profile& rhs) noexcept;

        Profile(const std::initializer_list<double>& il, char master) noexcept;

        Profile() = default;

        Profile& operator=(const Profile& rhs) noexcept;

        operator unsigned() const noexcept;

        static const Profile end_mark;

    private:
        static const double _THRESHOLD;

        double _relative_frequencies[pseudo::NUMBER];
        char   _master;
    };

}

template<typename InputIterator>
utils::Profile::Profile(InputIterator first, InputIterator last, size_t index) noexcept
{
    unsigned frequencies[pseudo::NUMBER];
    memset(frequencies, 0, sizeof(frequencies));

    size_t n = std::distance(first, last);
    for (size_t i = 0; i != n; ++i, ++first) ++frequencies[(*first)[index]];

    for (size_t i = 0; i != pseudo::NUMBER; ++i)
        _relative_frequencies[i] = static_cast<double>(frequencies[i]) / n;

    _master = 0;
    for (char i = 1; i != pseudo::NUMBER; ++i)
        if (_relative_frequencies[i] > _relative_frequencies[_master])
            _master = i;
}

