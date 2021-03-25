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

    std::cout << "finished, " << _sequences.size() << " sequences found" << std::endl;
    return _sequences.size();
}

