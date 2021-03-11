#pragma once

#include <string>
#include <vector>
#include <algorithm>

namespace utils
{

    std::string remove_white_spaces(const std::string& str);

    void err_exit(const std::initializer_list<std::string>& msg);

    inline unsigned char to_pseudo(char c);

    std::vector<unsigned char> to_pseudo(const std::string& str);

    template<typename InputIterator, typename OutputIterator>
    void transform_to_pseudo(InputIterator src_first, InputIterator src_last, OutputIterator des)
    {
        std::vector<unsigned char> (*op)(const std::string&) = &to_pseudo;
        std::transform(src_first, src_last, des, op);
    }

    template<typename InputIterator>
    InputIterator iter_of_max(InputIterator first, InputIterator last)
    {
        auto result = first;

        for (; first != last; ++first) if (*result < *first) result = first;
        return result;
    }

}
