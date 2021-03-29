#include "Utils.hpp"
#include "Pseudo.hpp"
#include "Fasta.hpp"

#include <regex>
#include <iostream>
#include <limits>

std::string utils::remove_white_spaces(const std::string& str)
{
    static const std::regex white_spaces("\\s+");
    return std::regex_replace(str, white_spaces, "");
}

void utils::err_exit(const std::string& info)
{
    std::cout << info << '\n';
    exit(1);
}

std::vector<unsigned char> utils::to_pseudo(const std::string& str)
{
    std::vector<unsigned char> pseu;
    pseu.reserve(str.size());

    for (auto i : str) pseu.push_back(to_pseudo(i));
    return pseu;
}

std::string utils::from_pseudo(const std::vector<unsigned char>& pseu)
{
    static constexpr char map[nucleic_acid_pseudo::NUMBER]{ '-', 'c', 'g', 'a', 't', 'n' };

    std::string str;
    str.reserve(pseu.size());

    for (auto i : pseu) str.push_back(map[i]);
    return str;
}

unsigned char* _get_map()
{
    using namespace nucleic_acid_pseudo;

    static unsigned char map[std::numeric_limits<unsigned char>::max()];
    memset(map, UNKNOWN, sizeof(map));

    // map['-'] = GAP; // we could not process sequences with '-'
    map['c'] = map['C'] = C;
    map['g'] = map['G'] = G;
    map['a'] = map['A'] = A;
    map['t'] = map['T'] = map['u'] = map['U'] = T;

    return map;
}

static const unsigned char* _map = _get_map();

unsigned char utils::to_pseudo(char ch)
{
    return _map[ch];
}

void utils::print_duration(std::chrono::system_clock::time_point time_point, const std::string& info)
{
    using namespace std::chrono;
    std::cout << info << ": " << duration_cast<microseconds>(system_clock::now() - time_point).count() << " us";
}

void utils::print_duration(std::chrono::system_clock::time_point time_point)
{
    using namespace std::chrono;
    std::cout << duration_cast<microseconds>(system_clock::now() - time_point).count();
}
