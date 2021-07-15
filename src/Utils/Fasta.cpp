#include "Fasta.hpp"

utils::Fasta::Fasta(std::istream& is) { _read(is); }

void utils::Fasta::write_to(std::ostream& os, bool with_identification) const
{
    if (with_identification)
        write_to(os, sequences.cbegin(), sequences.cend(), identifications.cbegin());
    else
        write_to(os, sequences.cbegin(), sequences.cend());
}

void utils::Fasta::_read(std::istream& is)
{
    std::string each_line, each_sequence;

    while (std::getline(is, each_line))
        if (each_line.size() && each_line[0] == '>')
        {
            identifications.push_back(each_line.substr(1));
            break;
        }

    while (std::getline(is, each_line))
    {
        if (each_line.size() == 0) continue;

        if (each_line[0] == '>')
        {
            identifications.push_back(each_line.substr(1));
            sequences.push_back(std::move(each_sequence));
        }
        else
        {
            each_sequence += each_line;
        }
    }

    sequences.push_back(each_sequence);
}
