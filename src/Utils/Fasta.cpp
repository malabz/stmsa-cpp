#include "Fasta.hpp"
#include "Utils.hpp"

#include <fstream>

utils::Fasta::Fasta(const std::string& infile) :
    size(_read(infile)),
    sequences(_sequences),
    identifications(_identifications)
{}

void utils::Fasta::write_to(std::ostream& os, bool with_identification) const
{
    if (with_identification)
        write_to(os, _sequences.cbegin(), _sequences.cend(), _identifications.cbegin());
    else
        write_to(os, _sequences.cbegin(), _sequences.cend());
}

template<typename InputIterator>
static void utils::Fasta::write_to(std::ostream& os, InputIterator sequence_first, InputIterator sequence_last)
{
    if (sequence_first == sequence_last) return;

    using difference_type = decltype(std::distance(sequence_first, sequence_last));
    const difference_type len = std::distance(sequence_first, sequence_last);
    for (difference_type i = 0; i != len - 1; ++sequence_first, ++i)
        os << *sequence_first << std::endl;
    os << *sequence_first;
}

template<typename InputIterator1, typename InputIterator2>
static void utils::Fasta::write_to(std::ostream &os, InputIterator1 sequence_first, InputIterator1 sequence_last,
    InputIterator2 identification_first)
{
    if (sequence_first == sequence_last) return;

    using difference_type = decltype(std::distance(sequence_first, sequence_last));
    const difference_type len = std::distance(sequence_first, sequence_last);
    for (difference_type i = 0; i != len - 1; ++sequence_first, ++identification_first, ++i)
    {
        os << '<' << *identification_first << std::endl;
        os << *sequence_first << std::endl;
    }
    os << '<' << *identification_first << std::endl;
    os << *sequence_first;
}

size_t utils::Fasta::_read(const std::string& infile)
{
    std::ifstream ifs(infile);

    if (!ifs) err_exit({ "cannot access file ", infile });
    std::cout << "reading " << infile << "..." << std::flush;

    std::string line, current;
    while (std::getline(ifs, line))
        if (line.size() && line[0] == '>')
        {
            _identifications.push_back(line.substr(1));
            break;
        }
    while (std::getline(ifs, line))
    {
        if (line.size() == 0) continue;

        if (line[0] == '>')
        {
            _identifications.push_back(line.substr(1));
            _sequences.push_back(std::move(current));
        }
        else
        {
            current += line;
        }
    }
    _sequences.push_back(current);

    std::cout << "finished" << std::endl;
    return _sequences.size();
}

