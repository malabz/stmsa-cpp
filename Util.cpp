#include "Util.h"

#include <string>
#include <regex>
#include <iostream>

inline std::string _remove_white_spaces(const std::string& str)
{
    static std::regex white_spaces("\s+");
    return std::regex_replace(str, white_spaces, "");
}

void err_exit(const std::initializer_list<std::string>& msg)
{
    for (const auto& i : msg) std::cout << i; std::cout << std::endl;
    exit(1);
}

