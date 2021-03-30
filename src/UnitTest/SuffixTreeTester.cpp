#include "SuffixTreeTester.hpp"
#include "../StarAlignment/StarAligner.hpp"
#include "../Utils/Fasta.hpp"
#include "../Utils/Profile.hpp"
#include "../Utils/Utils.hpp"
#include "../Utils/Pseudo.hpp"

#include <random>
#include <cassert>
#include <numeric>
#include <iomanip>
#include <cmath>

void suffixtree::SuffixTreeTester::test()
{
    std::random_device rd;
    std::mt19937 generator(rd());
    constexpr size_t max_width = 1 << 4;

    for (size_t width = 1; width != max_width; ++width)
    {
        const size_t max_element = width - 1;
        std::uniform_int_distribution<> dist(0, max_element); // [0, max_element]

        for (size_t length = 1; length != static_cast<size_t>(1) << width; length <<= 1)
        {
            std::vector<unsigned char> word;
            word.reserve(length + 1);
            for (size_t i = 0; i != length; ++i) word.push_back(dist(generator));
            // std::copy(word.cbegin(), word.cend(), std::ostream_iterator<int>(std::cout)); std::cout << '\n';

            auto suffixes = SuffixTree<unsigned char>(word.cbegin(), word.cend(), max_element).get_all_suffixes();
            assert(suffixes.size() == length + 1);

            std::sort(suffixes.begin(), suffixes.end(),
                [](const auto& lhs, const auto& rhs) { return lhs.size() < rhs.size(); });

            word.push_back(max_element + 1);
            for (size_t i = 1; i <= length; ++i)
            {
                const std::vector<unsigned char>& lhs = suffixes[i - 1];
                assert(lhs.size() == i);

                const std::vector<unsigned char> rhs(word.cend() - i, word.cend());
                // std::copy(lhs.cbegin(), lhs.cend(), std::ostream_iterator<unsigned int>(std::cout)); std::cout << '\n';
                // std::copy(rhs.cbegin(), rhs.cend(), std::ostream_iterator<unsigned int>(std::cout)); std::cout << '\n';
                assert(lhs == rhs);
            }
        }
    }
}

void suffixtree::SuffixTreeTester::test_profile(const char* lhs_file_path, const char* rhs_file_path)
{
    using namespace std::chrono;

    // preprocess input
    auto time_point = system_clock::now();
    utils::Fasta lhs(lhs_file_path), rhs(rhs_file_path);
    std::vector<std::vector<unsigned char>> lhs_sequences(lhs.sequences.size()), rhs_sequences(rhs.sequences.size());
    utils::transform_to_pseudo(lhs.sequences.cbegin(), lhs.sequences.cend(), lhs_sequences.begin());
    utils::transform_to_pseudo(rhs.sequences.cbegin(), rhs.sequences.cend(), rhs_sequences.begin());
    utils::print_duration(time_point, "reading and mapping"); std::cout << '\n';
    // std::copy(lhs_sequences[0].cbegin(), lhs_sequences[0].cend(), std::ostream_iterator<int>(std::cout)); std::cout << '\n';
    // std::copy(rhs_sequences[0].cbegin(), rhs_sequences[0].cend(), std::ostream_iterator<int>(std::cout)); std::cout << '\n';

    // classic suffix tree
    // suffixtree::SuffixTree<unsigned char> char_st(lhs_sequences[0].cbegin(), lhs_sequences[0].cend(), nucleic_acid_pseudo::MAX_ELE);
    // auto char_result = char_st.search_for_prefix(rhs_sequences[0].cbegin(), rhs_sequences[0].cend(), 3);
    // std::copy(char_result.cbegin(), char_result.cend(), std::ostream_iterator<size_t>(std::cout, ",")); std::cout << '\n';

    // profiles
    time_point = system_clock::now();
    std::vector<utils::Profile<nucleic_acid_pseudo::MAX_ELE>> lhs_profiles, rhs_profiles;
    lhs_profiles.reserve(lhs_sequences[0].size() + 1); // suffix tree '$'
    rhs_profiles.reserve(rhs_sequences[0].size());
    for (size_t i = 0; i != lhs_sequences[0].size(); ++i) lhs_profiles.emplace_back(lhs_sequences.cbegin(), lhs_sequences.cend(), i);
    for (size_t i = 0; i != rhs_sequences[0].size(); ++i) rhs_profiles.emplace_back(rhs_sequences.cbegin(), rhs_sequences.cend(), i);
    utils::print_duration(time_point, "profiles"); std::cout << '\n';
    // std::copy(lhs_profiles.cbegin(), lhs_profiles.cend(), std::ostream_iterator<unsigned>(std::cout)); std::cout << '\n';
    // std::copy(rhs_profiles.cbegin(), rhs_profiles.cend(), std::ostream_iterator<unsigned>(std::cout)); std::cout << '\n';

    auto artificial_sequence = utils::abstract(lhs_sequences, lhs_sequences.size(), lhs_sequences[0].size(), nucleic_acid_pseudo::MAX_ELE);

    // build suffix tree
    // std::copy((*lhs_sequences.cbegin()).cbegin(), (*lhs_sequences.cbegin()).cend(), std::ostream_iterator<int>(std::cout)); std::cout << '\n';
    time_point = system_clock::now();
    suffixtree::SuffixTree<unsigned char> st(artificial_sequence.cbegin(), artificial_sequence.cend(), nucleic_acid_pseudo::MAX_ELE);
    utils::print_duration(time_point, "suffix tree"); std::cout << '\n';

    // search
    lhs_profiles.emplace_back();
    time_point = system_clock::now();
    auto result = star_alignment::StarAligner::_optimal_path
            (utils::get_similar_substrings(st, lhs_profiles.cbegin(), rhs_profiles.cbegin(), rhs_profiles.cend(), 2., 31));
    utils::print_duration(time_point, "search"); std::cout << '\n';
    std::cout << result.size() << " similar substring(s) found" << '\n';
    // if (result.size()) suffixtree::SuffixTreeTester::print_matrix(result, result.size(), result[0].size());

    // analyse
    if (result.size())
    {
        size_t coverage = 0;
        for (size_t i = 0; i != result.size(); ++i) coverage += result[i][2];
        double relative_coverage = static_cast<double>(coverage * 2) / (lhs_sequences[0].size() + rhs_sequences[0].size());
        std::cout << "         coverage = " << coverage << '\n';
        std::cout << "relative_coverage = " << relative_coverage << '\n';

        std::vector<double> differences; differences.reserve(coverage);
        double difference_max = 0;
        for (size_t i = 0; i != result.size(); ++i)
            for (size_t j = 0; j != result[i][2]; ++j)
            {
                differences.push_back(lhs_profiles[result[i][0] + j].differs_from(rhs_profiles[result[i][1] + j]));
                if (differences.back() > difference_max) difference_max = differences.back();
            }
        std::cout << "average = " << (std::accumulate(differences.cbegin(), differences.cend(), .0) / differences.size()) << '\n';
        std::cout << "    max = " << difference_max << '\n';

        // size_t statistics[100];
        // memset(statistics, 0, sizeof(statistics));
        // for (size_t i = 0; i != differences.size(); ++i)
        //     ++statistics[static_cast<unsigned>(std::floor(differences[i] * 100))];
        // std::copy(statistics, statistics + 100, std::ostream_iterator<size_t>(std::cout, " ")); std::cout << '\n';
    }
}

template<typename MatrixType>
void suffixtree::SuffixTreeTester::print_matrix(const MatrixType& matrix, size_t row, size_t clm)
{
    for (size_t i = 0; i != row; ++i)
    {
        for (size_t j = 0; j != clm; ++j)
            std::cout << std::right << std::setw(8) << matrix[i][j] << ' ';
        std::cout << '\n';
    }
}
