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
    std::string line, current;

    while (std::getline(is, line))
        if (line.size() && line[0] == '>')
        {
            identifications.push_back(line.substr(1));
            break;
        }

    while (std::getline(is, line))
    {
        if (line.size() == 0) continue;

        if (line[0] == '>')
        {
            identifications.push_back(line.substr(1));
            sequences.push_back(std::move(current));
        }
        else
        {
            current += line;
        }
    }

    sequences.push_back(current);
}
