#pragma once

#ifdef _DEBUG
    #include <iostream>
#endif

#include <vector>
#include <iterator>
#include <algorithm>
#include <cstring>
#include <array>
#include <limits>
// #include <unordered_set>

namespace suffixtree
{

    template<typename T>
    class SuffixTree
    {
    public:
        using value_type = T;

        class Node;
        friend class Node;

        const size_t width;
        Node* const root;

        const size_t length;
        const value_type* const word;

        template<typename InputIterator>
        SuffixTree(InputIterator first, InputIterator last, value_type max_element) :
            width(static_cast<size_t>(max_element) + 2), // max_element + 1 is the end mark '$'
            root(_new_node(-1, 0, 0, false)),

            length(std::distance(first, last) + 1),
            word(_copy_word(first, last)),

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
            delete []word;
            _delete_tree(root);
        }

        template<typename InputIterator>
        std::vector<size_t> search_for_prefix(InputIterator first, InputIterator last, size_t threshold) const
        {
            size_t common_prefix_length = 0;
            for (Node *last_node = root, *curr_node = last_node->children[*first]; ; last_node = curr_node, curr_node = curr_node->children[*first])
            {
                if (curr_node == nullptr)
                    return common_prefix_length < threshold ? std::vector<size_t>() : get_all_beginning_with(last_node, common_prefix_length);

                for (const value_type* lhs_begin = word + curr_node->first, * lhs_end = lhs_begin + curr_node->length;
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
            const size_t len = last - first;

            // auto set_array = new std::unordered_set<size_t>*[length]();
            for (size_t rhs_index = 0; rhs_index < len; )
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

                    for (size_t i = 1; i != found.size(); ++i)
                        identical_substrings.push_back(std::array<size_t, 3>({ found[i], rhs_index, found[0] }));

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

                    rhs_index += found[0] - threshold + 1;
                }
            }

            // for (size_t i = 0; i != length; ++i) if (set_array[i]) delete set_array[i];
            // delete[] set_array;
            return identical_substrings;
        }

        std::vector<std::vector<value_type>> get_all_suffixes() const
        {
            std::vector<std::vector<value_type>> suffixes;
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
        value_type* _copy_word(InputIterator first, InputIterator last) const
        {
            value_type* result = new value_type[length];
            *std::copy(first, last, result) = width - 1;
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

        Node* _new_node(size_t first, size_t len, size_t len_from_root, bool is_leef) const
        {
            return new Node(first, len, len_from_root, is_leef, width, this);
        }

        void _get_all_suffixes(const Node* root, std::vector<Node *>& stack, std::vector<std::vector<value_type>>& suffixes) const
        {
            if (root->children == nullptr)
            {
                std::vector<value_type> curr;
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

        // 活动点
        Node* _active_node;                     // 活动结点
        const value_type* _active_edge;         // 用字符代表的活动边
        size_t _active_length;                  // 活动长度, 只有在活动边非null时有效
        value_type _active_edge_aux;

        // 每趟循环临时变量
        size_t _curr_end;                       // 待显式插入后缀的区间, 注意这里所代表的区间形式为(begin_index, end_index)
        size_t _remainder;                      // 待插入后缀数量
        Node*  _previous;                       // 最后创建的节点, 后缀链接时使用
        value_type      _current;               // 当前待插入后缀的最后一个元素
    };

    template<typename value_type>
    class SuffixTree<value_type>::Node
    {
    public:
        using tree_type = SuffixTree<value_type>;
        friend class tree_type;

        const tree_type* const tree;
        const size_t width;
        Node** children;

        size_t first;
        size_t length;
        size_t length_from_root;

        ~Node()
        {
            if (children) delete children;
        }

    private:
        Node(size_t first, size_t len, size_t len_from_root, bool is_leef, size_t width, const tree_type* tree) :
            tree(tree),
            width(width),

            first(first), length(len),
            length_from_root(len_from_root),

            children(is_leef ? nullptr : new Node * [width]()), _suffix(nullptr)
        {}

        void _split()
        {
            Node* next = children[*tree->_active_edge];

            auto new_node = new Node(next->first, tree->_active_length, length_from_root + tree->_active_length, false, width, tree);
            new_node->children[tree->word[next->first + tree->_active_length]] = next;

            next->first += tree->_active_length;
            next->length -= tree->_active_length;

            children[*tree->_active_edge] = new_node;
        }

        Node* _suffix;
    };

}
