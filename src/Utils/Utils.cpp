#include "Utils.hpp"
#include "Pseudo.hpp"

#include <regex>
#include <iostream>
#include <limits>

inline std::string utils::remove_white_spaces(const std::string& str)
{
    static const std::regex white_spaces("\\s+");
    return std::regex_replace(str, white_spaces, "");
}

void utils::err_exit(const std::initializer_list<std::string>& msg)
{
    for (const auto& i : msg) std::cout << i; std::cout << std::endl;
    exit(1);
}

std::vector<unsigned char> utils::to_pseudo(const std::string& str)
{
    std::vector<unsigned char> pseu;
    pseu.reserve(str.size());

    for (auto i : str) pseu.push_back(to_pseudo(i));
    return pseu;
}

unsigned char* _get_map()
{
    using namespace nucleic_acid_pseudo;

    static unsigned char map[std::numeric_limits<unsigned char>::max()];
    memset(map, UNKNOWN, sizeof(map));

    map['-'] = GAP;
    map['c'] = map['C'] = C;
    map['g'] = map['G'] = G;
    map['a'] = map['A'] = A;
    map['t'] = map['T'] = map['u'] = map['U'] = T;

    return map;
}

static const unsigned char* _map = _get_map();

inline unsigned char utils::to_pseudo(char ch)
{
    return _map[ch];
}
