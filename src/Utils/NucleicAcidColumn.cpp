#include "Utils.hpp"
#include "NucleicAcidColumn.hpp"

#include <string>
#include <fstream>

void utils::NucleicAcidColumn::_split_unknown(double* relative_frequencies) noexcept
{
    using namespace nucleic_acid_pseudo;

    const double c = relative_frequencies[C];
    const double g = relative_frequencies[G];
    const double a = relative_frequencies[A];
    const double t = relative_frequencies[T];
    const double n = relative_frequencies[N];
    const double one = n / (c + g + a + t);

    relative_frequencies[C] += one * c;
    relative_frequencies[G] += one * g;
    relative_frequencies[A] += one * a;
    relative_frequencies[T] += one * t;
    relative_frequencies[N] = 0;
}

utils::NucleicAcidColumn::operator unsigned char() const noexcept
{
    static constexpr unsigned char base[]{ 0, 5, 25, 85, 205 };

    bool flags[5] = {};
    unsigned char index = base[_available_count - 1];
    for (size_t i = 0; i != _available_count; ++i)
    {
        unsigned cnt = 0;
        for (size_t j = 0; j != _feature[i]; ++j) if (!flags[j]) ++cnt;
        flags[_feature[i]] = true;

        for (unsigned char j = 0, permutation = 4 - i; j != _available_count - i - 1; ++j, --permutation)
            cnt *= permutation;
        index += cnt;
    }

    return index;
}
