#pragma once

#include <vector>

namespace pairwise_alignment
{

    void needleman_wunsh_test();

    void needleman_wunsh_test(const std::vector<unsigned char>& lhs, const std::vector<unsigned char>& rhs, unsigned flag);

}
