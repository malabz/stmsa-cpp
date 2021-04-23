#pragma once

#include <vector>
#include <iterator>
#include <algorithm>
#include <cstring>
#include <array>
#include <limits>

namespace suffixtree
{

    template<size_t width>
    class SuffixTree
    {
    public:
        class Node;
        friend class Node;

        template<typename InputIterator>
        SuffixTree(InputIterator first, InputIterator last, unsigned char end_mark) :
            root(_new_node(-1, 0, 0, false)),

            length(std::distance(first, last) + 1),
            word(_copy_word(first, last, end_mark)),

            _active_node(root),
            _active_edge(nullptr),
            _active_length(0),

            _remainder(0)
        {
            for (_curr_end = 1; _curr_end <= length; ++_curr_end)
            { // 每次向后处理一个字符
                ++_remainder;
                _previous = nullptr;
                _current = word[_curr_end - 1];
                for (size_t i = _remainder; i != 0; --i)
                {
                    if (_active_length == 0)
                    { // 活动长度为0
                        if (_active_node->children[_current] != nullptr)
                        { // 可隐式插入
                            _active_edge_aux = _current;
                            _active_edge = &_active_edge_aux;
                            _active_length = 1;
                            _link(_active_node);
                            _try_to_walk_down_by_one_step();
                            break;
                        }
                        else
                        { // 不可隐式插入
                            size_t new_edge_length = length - _curr_end + 1;
                            _active_node->children[_current] = _new_node(_curr_end - 1, new_edge_length,
                                _active_node->length_from_root + new_edge_length, true);
                            if (_active_node != root)
                            { // 活动结点为root时无需擦屁股
                                _link(_active_node);
                                _active_node = _active_node->_suffix == nullptr ? root : _active_node->_suffix;
                            }
                            --_remainder;
                        }
                    }
                    else if (word[_active_node->children[*_active_edge]->first + _active_length] == _current)
                    { // 活动长度非0且可隐式插入
                        ++_active_length;
                        _try_to_walk_down_by_one_step();
                        break;
                    }
                    else
                    { // 活动长度非0且不可隐式插入
                        _active_node->_split();
                        size_t new_edge_length = length - _curr_end + 1;
                        _active_node->children[*_active_edge]->children[_current] = _new_node(_curr_end - 1, new_edge_length,
                            _active_node->children[*_active_edge]->length_from_root + new_edge_length, true);
                        _link(_active_node->children[*_active_edge]);
                        if (_active_node == root) // 活动结点为根结点
                            _active_edge = --_active_length == 0 ? nullptr : word + _curr_end - _remainder + 1;
                        else // 活动结点非根结点
                            _active_node = _active_node->_suffix == nullptr ? root : _active_node->_suffix;
                        _try_to_walk_down();
                        --_remainder;
                    }
                }
            }
        }

        ~SuffixTree()
        {
            delete[] word;
        }

        template<typename InputIterator>
        std::vector<size_t> search_for_prefix(InputIterator first, InputIterator last, size_t threshold) const
        {
            size_t common_prefix_length = 0;
            for (Node *last_node = root, *curr_node = last_node->children[*first]; ; last_node = curr_node, curr_node = curr_node->children[*first])
            {
                if (curr_node == nullptr)
                    return common_prefix_length < threshold ? std::vector<size_t>() : get_all_beginning_with(last_node, common_prefix_length);

                for (const unsigned char* lhs_begin = word + curr_node->first, * lhs_end = lhs_begin + curr_node->length;
                    lhs_begin != lhs_end && first != last; ++lhs_begin, ++first, ++common_prefix_length)
                    if (*lhs_begin != *first)
                        return common_prefix_length < threshold ? std::vector<size_t>() : get_all_beginning_with(curr_node, common_prefix_length);

                if (curr_node->children == nullptr || first == last)
                    return common_prefix_length < threshold ? std::vector<size_t>() : get_all_beginning_with(curr_node, common_prefix_length);
            }
        }

        template<typename RandomAccessIterator>
        std::vector<std::array<size_t, 3>> get_identical_substrings(RandomAccessIterator first, RandomAccessIterator last, size_t threshold) const
        {
            std::vector<std::array<size_t, 3>> identical_substrings;
            const size_t rhs_len = last - first;

            // auto set_array = new std::unordered_set<size_t>*[length]();
            for (size_t rhs_index = 0; rhs_index < rhs_len; )
            {
                auto found = search_for_prefix(first + rhs_index, last, threshold);

                if (found.size() == 0)
                {
                    ++rhs_index;
                }
                else
                {
                    // const size_t rhs_bgn = rhs_index;
                    // const size_t rhs_end = rhs_index + found[0];

                    // {
                    //     const size_t lhs_bgn = found[i];
                    //     const size_t lhs_end = found[i] + found[0];

                    //     auto set = set_array[lhs_end];
                    //     if (set == nullptr || set->find(rhs_end) == set->cend())
                    //     {
                    //         identical_substrings.push_back(std::array<size_t, 3>({ lhs_bgn, rhs_bgn, found[0] }));

                    //         if (set == nullptr) set_array[lhs_end] = new std::unordered_set<size_t>();
                    //         set_array[lhs_end]->insert(rhs_end);
                    //     }
                    // }

                    for (size_t i = 1; i != found.size(); ++i)
                        identical_substrings.push_back(std::array<size_t, 3>({ found[i], rhs_index, found[0] }));

                    rhs_index += found[0] - threshold + 1;
                }
            }

            // for (size_t i = 0; i != length; ++i) if (set_array[i]) delete set_array[i];
            // delete[] set_array;
            return identical_substrings;
        }

        std::vector<std::vector<unsigned char>> get_all_suffixes() const
        {
            std::vector<std::vector<unsigned char>> suffixes;
            suffixes.reserve(length);
            std::vector<Node*> stack;
            _get_all_suffixes(root, stack, suffixes);
            return suffixes;
        }

        std::vector<size_t> get_all_beginning_with(Node* curr_root, size_t to_add) const
        {
            std::vector<size_t> ret{ to_add };
            get_all_beginning_with(curr_root, ret);
            return ret;
        }

        void get_all_beginning_with(Node* curr_root, std::vector<size_t>& vsz) const
        {
            if (curr_root->children == nullptr)
                vsz.push_back(length - curr_root->length_from_root);
            else
                for (size_t i = 0; i != width; ++i)
                    if (curr_root->children[i] != nullptr)
                        get_all_beginning_with(curr_root->children[i], vsz);
        }

    private:
        template<typename InputIterator>
        unsigned char* _copy_word(InputIterator first, InputIterator last, unsigned char end_mark) const
        {
            unsigned char* result = new unsigned char[length];
            *std::copy(first, last, result) = end_mark;
            return result;
        }

        void _try_to_walk_down_by_one_step() noexcept
        {
            if (_active_length == _active_node->children[*_active_edge]->length)
            {
                _active_node = _active_node->children[*_active_edge];
                _active_edge = nullptr;
                _active_length = 0;
            }
        }

        void _try_to_walk_down() noexcept
        {
            while (_active_edge != nullptr && _active_length >= _active_node->children[*_active_edge]->length)
            {
                _active_length -= _active_node->children[*_active_edge]->length;
                _active_node = _active_node->children[*_active_edge];
                _active_edge = _active_length == 0 ? nullptr : word + _curr_end - 1 - _active_length;
            }
        }

        void _link(Node* next) noexcept
        {
            if (_previous != nullptr) _previous->_suffix = next;
            _previous = next;
        }

        Node* _new_node(size_t first, size_t len, size_t len_from_root, bool is_leef)
        {
            return new (_allocator.allocate_node()) Node(first, len, len_from_root, is_leef, this);
        }

        void _get_all_suffixes(const Node* root, std::vector<Node *>& stack, std::vector<std::vector<unsigned char>>& suffixes) const
        {
            if (root->children == nullptr)
            {
                std::vector<unsigned char> curr;
                size_t len = 0;
                for (auto np : stack) len += np->length;
                curr.reserve(len);
                for (auto np : stack)
                    for (size_t i = np->first; i != np->first + np->length; ++i)
                        curr.push_back(word[i]);
                suffixes.push_back(std::move(curr));
            }
            else
            {
                for (size_t i = 0; i != width; ++i)
                    if (root->children[i] != nullptr)
                    {
                        stack.push_back(root->children[i]);
                        _get_all_suffixes(root->children[i], stack, suffixes);
                        stack.pop_back();
                    }
            }
        }

        void _delete_tree(Node* tree)
        {
            if (tree && tree->children)
                for (size_t i = 0; i != width; ++i)
                    _delete_tree(tree->children[i]);

            delete tree;
        }

    private:
        class Allocator;
        Allocator _allocator;

    public:
        Node* const root;

        const size_t length;
        const unsigned char* const word;

    private:
        // 活动点
        Node*                   _active_node;       // 活动结点
        const unsigned char*    _active_edge;       // 用字符代表的活动边
        size_t                  _active_length;     // 活动长度, 只有在活动边非null时有效
        unsigned char           _active_edge_aux;

        // 每趟循环临时变量
        size_t                  _curr_end;          // 待显式插入后缀的区间, 注意这里所代表的区间形式为(begin_index, end_index)
        size_t                  _remainder;         // 待插入后缀数量
        Node*                   _previous;          // 最后创建的节点, 后缀链接时使用
        unsigned char           _current;           // 当前待插入后缀的最后一个元素
    };

    template <size_t width>
    class SuffixTree<width>::Node
    {
    public:
        friend class SuffixTree<width>;

        Node(size_t first, size_t len, size_t len_from_root, bool is_leef, SuffixTree<width>* tree) :
            first(first), length(len),
            length_from_root(len_from_root),

            children(is_leef ? nullptr : tree->_allocator.allocate_array()),
            tree(tree),

            _suffix(nullptr)
        {}

    private:
        void _split()
        {
            Node* next = children[*tree->_active_edge];

            auto new_node = new (tree->_allocator.allocate_node())
                Node(next->first, tree->_active_length, length_from_root + tree->_active_length, false, tree);
            new_node->children[tree->word[next->first + tree->_active_length]] = next;

            next->first += tree->_active_length;
            next->length -= tree->_active_length;

            children[*tree->_active_edge] = new_node;
        }

    public:
        size_t first;
        size_t length;
        size_t length_from_root;

        Node** const children;
        SuffixTree<width>* const tree;

    private:
        Node* _suffix;
    };

    template<size_t width>
    class SuffixTree<width>::Allocator
    {
    public:
        Allocator() :
            _node_block_size(_initial_block_size),
            _arr_block_size(_initial_block_size),

            _node_index(_node_block_size),
            _arr_index(_arr_block_size)
        {}

        Node* allocate_node()
        {
            if (_node_index == _node_block_size)
                _allocate_node_block();

            return reinterpret_cast<Node*>(_node_blocks.back() + sizeof(Node) * _node_index++);
        }

        // allocate and set 0
        Node** allocate_array()
        {
            if (_arr_index == _arr_block_size)
                _allocate_arr_block();

            auto mem = _arr_blocks.back() + sizeof(Node*) * width * _arr_index++;
            memset(mem, 0, sizeof(Node*) * width);

            return reinterpret_cast<Node**>(mem);
        }

        ~Allocator()
        {
            for (auto i : _node_blocks) delete[] i;
            for (auto i : _arr_blocks)  delete[] i;
        }

    private:
        void _allocate_node_block()
        {
            _node_index = 0;
            _node_block_size <<= 1;
            _node_blocks.push_back(new char[sizeof(Node) * _node_block_size]);
        }

        void _allocate_arr_block()
        {
            _arr_index = 0;
            _arr_block_size += _arr_blocks.size();
            _arr_blocks.push_back(new char[sizeof(Node*) * width * _arr_block_size]);
        }

        static constexpr size_t _initial_block_size = 16;

        std::vector<char*> _node_blocks;
        std::vector<char*> _arr_blocks;

        size_t _node_block_size, _arr_block_size;
        size_t _node_index, _arr_index;
    };

}
