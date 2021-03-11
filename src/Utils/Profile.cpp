#include "Profile.hpp"

utils::Profile::Profile(const Profile& rhs) noexcept
{
    memcpy(this, &rhs, sizeof(Profile));
}

utils::Profile& utils::Profile::operator=(const Profile& rhs) noexcept
{
    memcpy(this, &rhs, sizeof(Profile));
    return *this;
}

utils::Profile::Profile(const std::initializer_list<double>& il, char master) noexcept : _master(master)
{
    std::copy_n(il.begin(), pseudo::NUMBER, _relative_frequencies);
}

bool utils::operator==(const Profile& lhs, const Profile& rhs) noexcept
{
    double diff = 0;
    for (size_t i = 0; i != pseudo::NUMBER - 1; ++i) // degree of freedom
        diff += std::abs(rhs._relative_frequencies[i] - lhs._relative_frequencies[i]);
    return diff < Profile::_THRESHOLD;
}

bool utils::operator!=(const Profile& lhs, const Profile& rhs) noexcept
{
    return !(lhs == rhs);
}

utils::Profile::operator unsigned() const noexcept
{
    return static_cast<unsigned>(_master);
}

const double utils::Profile::_THRESHOLD = .55;

const utils::Profile utils::Profile::end_mark(std::initializer_list<double>{ 2, 2, 2, 2, 2, 2 }, pseudo::NUMBER);

