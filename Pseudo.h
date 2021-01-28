#pragma once

#include <string>
#include <algorithm>

namespace pseudo
{

    enum PseudoChar : char
    {
        GAP = 0,
        C = 1,
        G = 2,
        A = 3,
        T = 4,
        UNKNOWN = 5,
        NUMBER = 6
    };

    inline char to_pseudo(char c);

    std::string to_pseudo(const std::string& str);

    template<typename InputIterator, typename OutputIterator>
    void tranform_to_pseudo(InputIterator src_first, InputIterator src_last, OutputIterator des)
    {
        std::string (*op)(const std::string&) = &to_pseudo;
        std::transform(src_first, src_last, des, op);
    }

    void inplace_tranform_to_pseudo();

}

