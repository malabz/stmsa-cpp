#include "LeftChildRightSiblingBinaryTree.hpp"
#include "../Utils/Utils.hpp"

#include <vector>
#include <array>
#include <iostream>
#include <fstream>

suffix_tree::LeftChildRightSiblingBinaryTree::LeftChildRightSiblingBinaryTree(const SuffixTreeType &st)
    : root(_build(st.root))
    , word(st.word)
    , length(st.length)
{}

auto suffix_tree::LeftChildRightSiblingBinaryTree::_build(const SuffixTreeType::Node *rhs)
    -> Node *
{
    Node *root = new Node();
    root->first = rhs->first;
    root->length = rhs->length;
    root->length_from_root = rhs->length_from_root;

    if (rhs->children)
    {
        std::vector<const SuffixTreeType::Node *> rhs_children;
        rhs_children.reserve(nucleic_acid_pseudo::NUMBER);
        for (unsigned i = 0; i != nucleic_acid_pseudo::NUMBER; ++i)
            if (rhs->children[i])
                rhs_children.push_back(rhs->children[i]);
        if (rhs_children.empty()) return nullptr;

        std::vector<Node *> lhs_children;
        lhs_children.reserve(rhs_children.size());
        for (auto i : rhs_children)
            lhs_children.push_back(_build(i));
        for (size_t i = 1; i != lhs_children.size(); ++i)
            lhs_children[i - 1]->next_sibling = lhs_children[i];

        root->first_child = lhs_children.front();
    }
    else
    {
        root->first_child = nullptr;
    }

    return root;
}

suffix_tree::LeftChildRightSiblingBinaryTree::~LeftChildRightSiblingBinaryTree()
{
    _clear(root);
}

void suffix_tree::LeftChildRightSiblingBinaryTree::_clear(Node *tree)
{
    std::vector<Node *> stack;
    for (auto next = tree->first_child; next; next = next->next_sibling)
        stack.push_back(next);

    for (; stack.size(); stack.pop_back())
        _clear(stack.back());

    delete tree;
}

template <typename RandomAccessIterator>
std::vector<std::array<size_t, 3>>
suffix_tree::LeftChildRightSiblingBinaryTree::get_common_substrings(RandomAccessIterator first, RandomAccessIterator last, size_t threshold) const
{
    std::vector<std::array<size_t, 3>> common_substrings;
    const size_t rhs_len = last - first;

    for (size_t rhs_index = 0; rhs_index < rhs_len; )
    {
        auto found = search_for_prefix(first + rhs_index, last, threshold);

        if (found.empty())
        {
            ++rhs_index;
        }
        else
        {
            for (size_t i = 1; i != found.size(); ++i)
                common_substrings.push_back(std::array<size_t, 3>({ found[i], rhs_index, found[0] }));

            rhs_index += found[0] - threshold + 1;
        }
    }

    return common_substrings;
}

template<typename InputIterator>
std::vector<size_t>
suffix_tree::LeftChildRightSiblingBinaryTree::search_for_prefix(InputIterator first, InputIterator last, size_t threshold) const
{
    size_t common_prefix_length = 0;
    // for (Node *last_node = root, *curr_node = last_node->children[*first]; ; last_node = curr_node, curr_node = curr_node->children[*first])
    for (Node *last_node = root, *curr_node = search_for_child(last_node, *first); ; last_node = curr_node, curr_node = search_for_child(last_node, *first))
    {
        if (curr_node == nullptr)
            return common_prefix_length < threshold ? std::vector<size_t>() : get_all_beginning_with(last_node, common_prefix_length);

        for (const unsigned char *lhs_begin = word + curr_node->first,  *lhs_end = lhs_begin + curr_node->length;
            lhs_begin != lhs_end && first != last; ++lhs_begin, ++first, ++common_prefix_length)
            if (*lhs_begin != *first)
                return common_prefix_length < threshold ? std::vector<size_t>() : get_all_beginning_with(curr_node, common_prefix_length);

        // if (curr_node->children == nullptr || first == last)
        if (curr_node->first_child == nullptr || first == last)
            return common_prefix_length < threshold ? std::vector<size_t>() : get_all_beginning_with(curr_node, common_prefix_length);
    }
}

std::vector<size_t>
suffix_tree::LeftChildRightSiblingBinaryTree::get_all_beginning_with(Node *curr_root, size_t to_add) const
{
    std::vector<size_t> ret{ to_add };
    get_all_beginning_with(curr_root, ret);
    return ret;
}

void suffix_tree::LeftChildRightSiblingBinaryTree::get_all_beginning_with(Node *curr_root, std::vector<size_t> &vsz) const
{
    if (curr_root->first_child == nullptr)
        vsz.push_back(length - curr_root->length_from_root);
    else
        for (auto next = curr_root->first_child; next; next = next->next_sibling)
            get_all_beginning_with(next, vsz);
}

auto suffix_tree::LeftChildRightSiblingBinaryTree::search_for_child(const Node *parent, unsigned char target) const noexcept
    ->Node *
{
    for (auto next = parent->first_child; next; next = next->next_sibling)
        if (word[next->first] == target)
            return next;

    return nullptr;
}

size_t suffix_tree::LeftChildRightSiblingBinaryTree::_count_nodes(const Node *tree) noexcept
{
    size_t count = 1;
    for (auto next = tree->first_child; next; next = next->next_sibling)
        count += _count_nodes(next);
    return count;
}

size_t suffix_tree::LeftChildRightSiblingBinaryTree::count_nodes() const noexcept
{
    return _count_nodes(root);
}

void suffix_tree::experiment(const char *in_file, const char *out_file)
{
    using SuffixTreeType = suffix_tree::LeftChildRightSiblingBinaryTree::SuffixTreeType;

    std::ifstream ifs(in_file);
    if (!ifs)
    {
        std::cout << "cannot access file " << in_file << '\n';
        exit(0);
    }

    const auto sequences = utils::read_to_pseudo(ifs);
    std::cout << sequences.size() << " sequences found\n";
    if (sequences.size() < 2) return;

    std::vector<long long> time_consumtion[2];
    time_consumtion[0].reserve(sequences.size());
    time_consumtion[1].reserve(sequences.size());

    std::vector<size_t> memory_consumtion[2];
    memory_consumtion[0].reserve(sequences.size());
    memory_consumtion[1].reserve(sequences.size());

    for (size_t i = 0; i != sequences.size(); ++i)
    {
        using namespace std::chrono;
        decltype(high_resolution_clock::now()) time_point;

        time_point = high_resolution_clock::now();
        SuffixTreeType st(sequences[i].cbegin(), sequences[i].cend(), nucleic_acid_pseudo::GAP);
        for (const auto &i : sequences) st.get_common_substrings(i.cbegin(), i.cend(), 15);
        time_consumtion[0].push_back(duration_cast<microseconds>(high_resolution_clock::now() - time_point).count());

        time_point = high_resolution_clock::now();
        suffix_tree::LeftChildRightSiblingBinaryTree lcrsbt(st);
        for (const auto &i : sequences) lcrsbt.get_common_substrings(i.cbegin(), i.cend(), 15).size();
        time_consumtion[1].push_back(duration_cast<microseconds>(high_resolution_clock::now() - time_point).count());

        const size_t num_node = lcrsbt.count_nodes();
        memory_consumtion[1].push_back(num_node * sizeof(suffix_tree::LeftChildRightSiblingBinaryTree::Node));
        memory_consumtion[0].push_back(SuffixTreeType::memory_consumtion(num_node, st.length));

        if (i != 0)
        {
            size_t prev_len = std::to_string(i - 1).size();
            for (size_t i = 0; i != prev_len; ++i)
                std::cout << '\b';
        }
        std::cout << std::to_string(i) << std::flush;
    }

    std::ofstream ofs(out_file);
    if (!ofs)
    {
        std::cout << "cannot write to file " << out_file << '\n';
        exit(0);
    }

    for (size_t i = 0; i != sequences.size(); ++i)
        ofs << time_consumtion[0][i] << ',' << memory_consumtion[0][i] << ','
            << time_consumtion[1][i] << ',' << memory_consumtion[1][i] << '\n';
}
