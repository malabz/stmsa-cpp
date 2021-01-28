#pragma once

#include <string>
#include <vector>
#include <iostream>

class Fasta
{
private:
    std::vector<std::string> _sequences;
    std::vector<std::string> _identifications;

    const size_t size;

    size_t _read(const std::string& filename);

public:
    const std::vector<std::string>& sequences;
    const std::vector<std::string>& identifications;

    explicit Fasta(const std::string& infile);

    void write_to(std::ostream& os, bool with_idification = true) const;

    template<typename InputIterator>
    static void write_to(std::ostream& os, InputIterator sequence_first, InputIterator sequence_last);

    template<typename InputIterator1, typename InputIterator2>
    static void write_to(std::ostream& os, InputIterator1 sequence_first, InputIterator1 sequence_last,
        InputIterator2 identification_first);
};

