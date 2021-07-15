#include "Utils.hpp"
#include "Pseudo.hpp"
#include "Fasta.hpp"

#include <regex>
#include <iostream>
#include <limits>

std::string utils::remove_white_spaces(const std::string &str)
{
    static const std::regex white_spaces("\\s+");
    return std::regex_replace(str, white_spaces, "");
}

std::vector<unsigned char> utils::to_pseudo(const std::string &str)
{
    std::vector<unsigned char> pseu;
    pseu.reserve(str.size());

    for (auto i : str) pseu.push_back(to_pseudo(i));
    return pseu;
}

std::string utils::from_pseudo(const std::vector<unsigned char> &pseu)
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
    memset(map, N, sizeof(map));

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

void utils::print_duration(std::chrono::system_clock::time_point time_point, const std::string &info)
{
    using namespace std::chrono;
    std::cout << info << ": " << duration_cast<microseconds>(system_clock::now() - time_point).count() << " us";
}

void utils::print_duration(std::chrono::system_clock::time_point time_point)
{
    using namespace std::chrono;
    std::cout << duration_cast<microseconds>(system_clock::now() - time_point).count();
}

std::vector<std::vector<unsigned char>> utils::read_to_pseudo(std::istream &is)
{
    std::vector<std::vector<unsigned char>> sequences;

    std::string each_line, each_sequence;

    while (std::getline(is, each_line))
        if (each_line.size() && each_line[0] == '>')
            break;

    while (std::getline(is, each_line))
    {
        if (each_line.size() == 0) continue;

        if (each_line[0] == '>')
        {
            sequences.push_back(to_pseudo(each_sequence));
            each_sequence.clear();
        }
        else
        {
            each_sequence += each_line;
        }
    }
    sequences.push_back(to_pseudo(each_sequence));

    return sequences;
}

void utils::insert_and_write(std::ostream &os, std::istream &is, const std::vector<std::vector<Insertion>> &insertions)
{
    const size_t sequence_number = insertions.size();
    std::string each_line, each_sequence;

    while (std::getline(is, each_line))
        if (each_line.size() && each_line[0] == '>')
        {
            os << each_line << '\n';
            break;
        }

    std::string each_sequence_aligned;
    for (unsigned count = 0, length; std::getline(is, each_line); )
    {
        if (each_line.size() == 0) continue;

        if (each_line[0] == '>')
        {
            if (count == 0)
            {
                length = each_sequence.size();
                for (auto insertion : insertions[0])
                    length += insertion.number;
                each_sequence_aligned.reserve(length);
                each_sequence.reserve(length);
            }

            utils::Insertion::insert_gaps(each_sequence.cbegin(), each_sequence.cend(),
                    insertions[count].cbegin(), insertions[count].cend(), std::back_inserter(each_sequence_aligned), '-');
            os << each_sequence_aligned;
            if (++count != sequence_number) os << '\n';

            each_sequence.clear();
            each_sequence_aligned.clear();

            os << each_line << '\n';
        }
        else
        {
            each_sequence += each_line;
        }
    }

    utils::Insertion::insert_gaps(each_sequence.cbegin(), each_sequence.cend(),
            insertions.back().cbegin(), insertions.back().cend(), std::back_inserter(each_sequence_aligned), '-');
    os << each_sequence_aligned << '\n';
}
