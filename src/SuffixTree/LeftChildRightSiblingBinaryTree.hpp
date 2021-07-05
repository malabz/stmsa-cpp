#include "SuffixTree.hpp"
#include "../Utils/NucleicAcidColumn.hpp"

namespace suffix_tree
{
    class LeftChildRightSiblingBinaryTree
    {
    public:
        using SuffixTreeType = SuffixTree<nucleic_acid_pseudo::NUMBER>;

        struct Node
        {
            size_t first;
            size_t length;
            size_t length_from_root;

            Node *first_child;
            Node *next_sibling;

            LeftChildRightSiblingBinaryTree *tree;

            Node *_suffix_link;
        };

        LeftChildRightSiblingBinaryTree(const SuffixTreeType &st);
        ~LeftChildRightSiblingBinaryTree();

        template <typename RandomAccessIterator>
        std::vector<std::array<size_t, 3>> get_identical_substrings(RandomAccessIterator first, RandomAccessIterator last, size_t threshold) const;

        template<typename InputIterator>
        std::vector<size_t> search_for_prefix(InputIterator first, InputIterator last, size_t threshold) const;

        std::vector<size_t> get_all_beginning_with(Node *curr_root, size_t to_add) const;
        void get_all_beginning_with(Node *curr_root, std::vector<size_t> &vsz) const;

        Node *search_for_child(const Node *parent, unsigned char target) const noexcept;

        size_t count_nodes() const noexcept;

        Node *const root;

        const unsigned char *const word;
        const size_t length;

    private:
        static Node *_build(const SuffixTreeType::Node *rhs);

        static size_t _count_nodes(const Node *tree) noexcept;

        void _clear(Node *tree);
    };

    // time and memory consumption difference
    void experiment(const char *in_file, const char *out_file);
}
