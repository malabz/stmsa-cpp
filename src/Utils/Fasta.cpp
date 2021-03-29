#include "Fasta.hpp"
#include "Utils.hpp"

#include <fstream>

utils::Fasta::Fasta(const std::string& infile) :
    size(_read(infile))
{}

void utils::Fasta::write_to(std::ostream& os, bool with_identification) const
{
    if (with_identification)
        write_to(os, sequences.cbegin(), sequences.cend(), identifications.cbegin());
    else
        write_to(os, sequences.cbegin(), sequences.cend());
}

size_t utils::Fasta::_read(const std::string& infile)
{
    std::ifstream ifs(infile);

    if (!ifs) err_exit("cannot access file " + infile);
    std::cout << "reading " << infile << "..." << std::flush;

    std::string line, current;

    while (std::getline(ifs, line))
        if (line.size() && line[0] == '>')
        {
            identifications.push_back(line.substr(1));
            break;
        }

    while (std::getline(ifs, line))
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

    std::cout << "finished, " << sequences.size() << " sequences found" << std::endl;
    return sequences.size();
}

std::tuple<std::vector<std::string>, std::vector<std::vector<unsigned char>>>
utils::Fasta::read_to_pseudo(const std::string& infile)
{
    std::ifstream ifs(infile);

    if (!ifs) err_exit("cannot access file " + infile);
    std::cout << "reading " << infile << "..." << std::flush;

    std::vector<std::string> identifications;
    std::vector<std::vector<unsigned char>> sequences;

    std::string line, current;

    while (std::getline(ifs, line))
        if (line.size() && line[0] == '>')
        {
            identifications.push_back(line.substr(1));
            break;
        }

    while (std::getline(ifs, line))
    {
        if (line.size() == 0) continue;

        if (line[0] == '>')
        {
            identifications.push_back(line.substr(1));
            sequences.push_back(to_pseudo(current));
            current.clear();
        }
        else
        {
            current += line;
        }
    }

    sequences.push_back(to_pseudo(current));

    std::cout << "finished, " << sequences.size() << " sequences found" << std::endl;
    return std::make_tuple(std::move(identifications), std::move(sequences));
}
