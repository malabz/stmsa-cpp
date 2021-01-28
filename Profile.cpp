#include "Profile.h"

Profile::Profile(const Profile& rhs) noexcept
{
    memcpy(this, &rhs, sizeof(Profile));
}

Profile& Profile::operator=(const Profile& rhs) noexcept
{
    memcpy(this, &rhs, sizeof(Profile));
    return *this;
}

Profile::Profile(const std::initializer_list<double>& il, char prevail) noexcept : _prevail(prevail)
{
    std::copy_n(il.begin(), NUMBER, _relative_frequencies);
}

bool operator==(const Profile& lhs, const Profile& rhs) noexcept
{
    double diff = 0;
    for (size_t i = 0; i != NUMBER - 1; ++i) // degree of freedom
        diff += std::abs(rhs._relative_frequencies[i] - lhs._relative_frequencies[i]);
    return diff < Profile::_THRESHOLD;
}

bool operator!=(const Profile& lhs, const Profile& rhs) noexcept
{
    return !(lhs == rhs);
}

Profile::operator unsigned() const noexcept
{
    return static_cast<unsigned>(_prevail);
}

const double Profile::_THRESHOLD = .2;

const Profile Profile::end_mark(std::initializer_list<double>{ 2, 2, 2, 2, 2, 2 }, NUMBER);

