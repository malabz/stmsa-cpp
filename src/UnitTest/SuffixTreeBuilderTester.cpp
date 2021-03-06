#include "../SuffixTree/SuffixTreeBuilder.hpp"
#include "../UnitTest/SuffixTreeBuilderTester.hpp"
#include "../Utils/Fasta.hpp"
#include "../Utils/Profile.hpp"
#include "../Utils/Utils.hpp"
#include "../Utils/Pseudo.hpp"

#include <random>
#include <cassert>
#include <iterator>

void suffixtree::SuffixTreeBuilderTester::test()
{
    std::random_device rd;
    std::mt19937 generator(rd());
    for (size_t width = 1; width != 1 << 3; ++width)
    {
        std::uniform_int_distribution<> dist(0, width - 1); // inclusive
        for (size_t length = 1 << 0; length != static_cast<size_t>(1) << (width + 8); length <<= 1)
            for (int repeat = 0; repeat != (width == 1 ? 1 : 1 << 4); ++repeat)
            {
                std::vector<char> word;
                word.reserve(length);
                for (size_t i = 0; i != length; ++i) word.push_back(dist(generator));
                // std::copy(word.cbegin(), word.cend(), std::ostream_iterator<int>(std::cout)); std::cout << std::endl;

                auto suffixes = SuffixTreeBuilder<char>(word.cbegin(), word.cend(), width, width)._get_all_suffixes();
                assert(suffixes.size() == length + 1);

                std::sort(suffixes.begin(), suffixes.end(),
                    [](const std::vector<char>& lhs, const std::vector<char>& rhs) { return lhs.size() < rhs.size(); });
                word.push_back(width);

                for (size_t i = 1; i != length + 1; ++i)
                {
                    const std::vector<char>& lhs = suffixes[i - 1];
                    assert(lhs.size() == i);

                    const std::vector<char> rhs(word.cend() - i, word.cend());
                    assert(lhs == rhs);
                    // for (size_t j = 0; j != lhs.size(); ++j) assert(lhs[j] == rhs[j]);
                }
            }
    }
}

void suffixtree::SuffixTreeBuilderTester::test_profile(const char* lhs_file_path, const char* rhs_file_path)
{
    using namespace utils;
    using namespace pseudo;

    // preprocess input
    Fasta lhs(lhs_file_path), rhs(rhs_file_path);
    std::vector<std::vector<unsigned char>> lhs_sequences, rhs_sequences;
    lhs_sequences.reserve(lhs.sequences.size());
    rhs_sequences.reserve(rhs.sequences.size());
    transform_to_pseudo(lhs.sequences.cbegin(), lhs.sequences.cend(), std::back_inserter(lhs_sequences));
    transform_to_pseudo(rhs.sequences.cbegin(), rhs.sequences.cend(), std::back_inserter(rhs_sequences));
    std::copy(lhs_sequences[0].cbegin(), lhs_sequences[0].cend(), std::ostream_iterator<int>(std::cout)); std::cout << std::endl;
    std::copy(rhs_sequences[0].cbegin(), rhs_sequences[0].cend(), std::ostream_iterator<int>(std::cout)); std::cout << std::endl;

    // classic suffix tree
    suffixtree::SuffixTreeBuilder<char> char_stb(lhs_sequences[0].cbegin(), lhs_sequences[0].cend(), NUMBER, NUMBER);
    auto char_result = char_stb.search_for_prefix(rhs_sequences[0].cbegin(), rhs_sequences[0].cend(), 3);
    std::copy(char_result.cbegin(), char_result.cend(), std::ostream_iterator<size_t>(std::cout, ",")); std::cout << std::endl;

    // generate profiles
    std::vector<Profile> lhs_profiles, rhs_profiles;
    lhs_profiles.reserve(lhs_sequences[0].size());
    lhs_profiles.reserve(rhs_sequences[0].size());
    for (size_t i = 0; i != lhs_sequences[0].size(); ++i) lhs_profiles.emplace_back(lhs_sequences.cbegin(), lhs_sequences.cend(), i);
    for (size_t i = 0; i != rhs_sequences[0].size(); ++i) rhs_profiles.emplace_back(rhs_sequences.cbegin(), rhs_sequences.cend(), i);
    std::copy(lhs_profiles.cbegin(), lhs_profiles.cend(), std::ostream_iterator<unsigned>(std::cout)); std::cout << std::endl;
    std::copy(rhs_profiles.cbegin(), rhs_profiles.cend(), std::ostream_iterator<unsigned>(std::cout)); std::cout << std::endl;

    // build suffix tree
    std::copy((*lhs_sequences.cbegin()).cbegin(), (*lhs_sequences.cbegin()).cend(), std::ostream_iterator<int>(std::cout)); std::cout << std::endl;
    suffixtree::SuffixTreeBuilder<Profile> profile_stb(lhs_profiles.cbegin(), lhs_profiles.cend(), NUMBER, Profile::end_mark);

    // test by suffixes
    std::vector<std::vector<Profile>> suffixes = profile_stb._get_all_suffixes();
    std::sort(suffixes.begin(), suffixes.end(), [](const std::vector<Profile>& lhs, const std::vector<Profile>& rhs) { return lhs.size() < rhs.size(); });
    for (size_t i = 0; i != suffixes.size(); ++i)
    { std::copy(suffixes[i].cbegin(), suffixes[i].cend(), std::ostream_iterator<unsigned>(std::cout)); std::cout << std::endl; }

    // search
    auto result = profile_stb.search_for_prefix(rhs_profiles.cbegin(), rhs_profiles.cbegin() + 7, 2);
    std::copy(result.cbegin(), result.cend(), std::ostream_iterator<size_t>(std::cout, ",")); std::cout << std::endl;
}

void suffixtree::SuffixTreeBuilderTester::test_substring(const char* lhs_file_path, const char* rhs_file_path)
{
    using namespace utils;
    using namespace pseudo;

    // preprocess input
    Fasta lhs(lhs_file_path), rhs(rhs_file_path);
    std::vector<std::vector<unsigned char>> lhs_sequences, rhs_sequences;
    lhs_sequences.reserve(lhs.sequences.size());
    rhs_sequences.reserve(rhs.sequences.size());
    transform_to_pseudo(lhs.sequences.cbegin(), lhs.sequences.cend(), std::back_inserter(lhs_sequences));
    transform_to_pseudo(rhs.sequences.cbegin(), rhs.sequences.cend(), std::back_inserter(rhs_sequences));
    // std::copy(lhs_sequences[0].cbegin(), lhs_sequences[0].cend(), std::ostream_iterator<int>(std::cout)); std::cout << std::endl;
    // std::copy(rhs_sequences[0].cbegin(), rhs_sequences[0].cend(), std::ostream_iterator<int>(std::cout)); std::cout << std::endl;

    // generate profiles
    std::cout << "generating profile..." << std::flush;
    std::vector<Profile> lhs_profiles, rhs_profiles;
    lhs_profiles.reserve(lhs_sequences[0].size());
    lhs_profiles.reserve(rhs_sequences[0].size());
    for (size_t i = 0; i != lhs_sequences[0].size(); ++i) lhs_profiles.emplace_back(lhs_sequences.cbegin(), lhs_sequences.cend(), i);
    for (size_t i = 0; i != rhs_sequences[0].size(); ++i) rhs_profiles.emplace_back(rhs_sequences.cbegin(), rhs_sequences.cend(), i);
    // std::copy(lhs_profiles.cbegin(), lhs_profiles.cend(), std::ostream_iterator<unsigned>(std::cout)); std::cout << std::endl;
    // std::copy(rhs_profiles.cbegin(), rhs_profiles.cend(), std::ostream_iterator<unsigned>(std::cout)); std::cout << std::endl;
    std::cout << "finished" << std::endl;

    // build suffix tree
    std::cout << "building suffix tree..." << std::flush;
    // std::copy((*lhs_sequences.cbegin()).cbegin(), (*lhs_sequences.cbegin()).cend(), std::ostream_iterator<int>(std::cout)); std::cout << std::endl;
    suffixtree::SuffixTreeBuilder<Profile> profile_stb(lhs_profiles.cbegin(), lhs_profiles.cend(), NUMBER, Profile::end_mark);
    std::cout << "finished" << std::endl;

    // test by suffixes
    // std::vector<std::vector<Profile>> suffixes = profile_stb._get_all_suffixes();
    // std::sort(suffixes.begin(), suffixes.end(), [](const std::vector<Profile>& lhs, const std::vector<Profile>& rhs) { return lhs.size() < rhs.size(); });
    // for (size_t i = 0; i != suffixes.size(); ++i)
    // { std::copy(suffixes[i].cbegin(), suffixes[i].cend(), std::ostream_iterator<unsigned>(std::cout)); std::cout << std::endl; }

    // search
    std::cout << "searching..." << std::flush;
    auto result = profile_stb.get_identical_substring_with(rhs_profiles.cbegin(), rhs_profiles.cend(), 15);
    std::cout << "finished" << std::endl;
    std::cout << result.size() << std::endl;
    for (size_t i = 0; i != result.size(); ++i)
    {
        std::cout << '[';
        std::copy(result[i].cbegin(), result[i].cend(), std::ostream_iterator<size_t>(std::cout, ", "));
        std::cout << ']' << std::endl;
    }
}

template<typename T>
void suffixtree::SuffixTreeBuilderTester::print(const std::vector<std::vector<T>>& vt)
{
    for (size_t i = 0; i != vt.size(); ++i)
    {
        std::copy(vt[i].cbegin(), vt[i].cend(), std::ostream_iterator<T>(std::cout));
        std::cout << std::endl;
    }
}

