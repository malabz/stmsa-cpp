#pragma once

#include <string>
#include <vector>
#include <iostream>

namespace utils
{

    class Fasta
    {
    private:
        size_t _read(const std::string& filename);

        std::vector<std::string> _sequences;
        std::vector<std::string> _identifications;
        const size_t size;

    public:
        const std::vector<std::string>& sequences;
        const std::vector<std::string>& identifications;

        explicit Fasta(const std::string& infile);

        void write_to(std::ostream& os, bool with_idification = true) const;

        template<typename InputIterator>
        static void write_to(std::ostream& os, InputIterator sequence_first, InputIterator sequence_last)
        {
            if (sequence_first == sequence_last) return;

            using difference_type = decltype(std::distance(sequence_first, sequence_last));
            const difference_type len = std::distance(sequence_first, sequence_last);

            for (difference_type i = 0; i != len - 1; ++sequence_first, ++i)
                os << *sequence_first << '\n';
            os << *sequence_first;
        }

        template<typename InputIterator1, typename InputIterator2>
        static void write_to(std::ostream& os, InputIterator1 sequence_first, InputIterator1 sequence_last,
            InputIterator2 identification_first)
        {
            if (sequence_first == sequence_last) return;

            using difference_type = decltype(std::distance(sequence_first, sequence_last));
            const difference_type len = std::distance(sequence_first, sequence_last);

            for (difference_type i = 0; i != len - 1; ++sequence_first, ++identification_first, ++i)
            {
                os << '<' << *identification_first << '\n';
                os << *sequence_first << '\n';
            }
            os << '<' << *identification_first << '\n';
            os << *sequence_first;
        }

    };

}
