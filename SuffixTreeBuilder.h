#pragma

#include <vector>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <cstring>

namespace suffixtree
{

    template<typename T>
    class SuffixTreeBuilder
    {
    public:
        template<typename InputIterator>
        SuffixTreeBuilder(InputIterator first, InputIterator last, size_t element_count, const T& end_mark) :
            _length(std::distance(first, last) + 1), _word(_copy_word(first, last, end_mark)),
            _width(element_count + 1), _end_mark(end_mark), _root(_new_node(-1, 0, 0, false)),
            _active_node(_root), _active_edge(nullptr), _active_length(0), _remainder(0)
        {
            for (_curr_end = 1; _curr_end != _length + 1; ++_curr_end)
            { // 每次向后处理一个字符
                ++_remainder;
                _previous = nullptr;
                _current = _word[_curr_end - 1];
                for (size_t i = _remainder; i != 0; --i)
                {
                    if (_active_length == 0)
                    { // 活动长度为0
                        if (_active_node->_children[_current] != nullptr)
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
                            size_t new_edge_length = _length - _curr_end + 1;
                            _active_node->_children[_current] = _new_node(_curr_end - 1, new_edge_length,
                                _active_node->_length_from_root + new_edge_length, true);
                            if (_active_node != _root)
                            { // 活动结点为root时无需擦屁股
                                _link(_active_node);
                                _active_node = _active_node->_suffix == nullptr ? _root : _active_node->_suffix;
                            }
                            --_remainder;
                        }
                    }
                    else if (_word[_active_node->_children[*_active_edge]->_first + _active_length] == _current)
                    { // 活动长度非0且可隐式插入
                        ++_active_length;
                        _try_to_walk_down_by_one_step();
                        break;
                    }
                    else
                    { // 活动长度非0且不可隐式插入
                        _active_node->_split();
                        size_t new_edge_length = _length - _curr_end + 1;
                        _active_node->_children[*_active_edge]->_children[_current] = _new_node(_curr_end - 1, new_edge_length,
                            _active_node->_children[*_active_edge]->_length_from_root + new_edge_length, true);
                        _link(_active_node->_children[*_active_edge]);
                        if (_active_node == _root) // 活动结点为根结点
                            _active_edge = --_active_length == 0 ? nullptr : _word + _curr_end - _remainder + 1;
                        else // 活动结点非根结点
                            _active_node = _active_node->_suffix == nullptr ? _root : _active_node->_suffix;
                        _try_to_walk_down();
                        --_remainder;
                    }
                }
            }
        }

        ~SuffixTreeBuilder()
        {
            delete []_word;
        }

        template<typename InputIterator>
        std::vector<size_t> search_for_prefix(InputIterator first, InputIterator last, size_t threshold) const
        {
            size_t common_prefix_length = 0;
            for (_Node *last_node = _root, *curr_node = last_node->_children[*first]; ; last_node = curr_node, curr_node = curr_node->_children[*first])
            {
                if (curr_node == nullptr)
                    return common_prefix_length < threshold ? std::vector<size_t>() : _get_all_beginning_with(last_node, common_prefix_length);

                for (const T* lhs_first = _word + curr_node->_first, *lhs_end = lhs_first + curr_node->_length;
                    lhs_first != lhs_end && first != last; ++lhs_first, ++first, ++common_prefix_length)
                    if (*lhs_first != *first)
                        return common_prefix_length < threshold ? std::vector<size_t>() : _get_all_beginning_with(curr_node, common_prefix_length);

                if (curr_node->_children == nullptr || first == last)
                    return common_prefix_length < threshold ? std::vector<size_t>() : _get_all_beginning_with(curr_node, common_prefix_length);
            }
        }

    private:

        class _Node;

        friend class _Node;
        friend class SuffixTreeBuilderTester;

        template<typename InputIterator>
        T* _copy_word(InputIterator first, InputIterator last, const T& end_mark)
        {
            T* result = new T[_length + 1];
            *std::copy(first, last, result) = end_mark;
            return result;
        }

        void _try_to_walk_down_by_one_step() noexcept
        {
            if (_active_length == _active_node->_children[*_active_edge]->_length)
            {
                _active_node = _active_node->_children[*_active_edge];
                _active_edge = nullptr;
                _active_length = 0;
            }
        }

        void _try_to_walk_down() noexcept
        {
            while (_active_edge != nullptr && _active_length >= _active_node->_children[*_active_edge]->_length)
            {
                _active_length -= _active_node->_children[*_active_edge]->_length;
                _active_node = _active_node->_children[*_active_edge];
                _active_edge = _active_length == 0 ? nullptr : _word + _curr_end - 1 - _active_length;
            }
        }

        void _link(_Node* next) noexcept
        {
            if (_previous != nullptr) _previous->_suffix = next;
            _previous = next;
        }

        inline _Node* _new_node(size_t first, size_t len, size_t len_from_root, bool is_leef) const
        {
            return new _Node(first, len, len_from_root, is_leef, _width, this);
        }

        std::vector<std::vector<T>> _get_all_suffixes() const
        {
            std::vector<std::vector<T>> suffixes;
            suffixes.reserve(_length);
            std::vector<_Node*> stack;
            _get_all_suffixes(_root, stack, suffixes);
            return suffixes;
        }

        void _get_all_suffixes(const _Node* root, std::vector<_Node *>& stack, std::vector<std::vector<T>>& suffixes) const
        {
            if (root->_children == nullptr)
            {
                std::vector<T> curr;
                size_t len = 0;
                for (auto np : stack) len += np->_length;
                curr.reserve(len);
                for (auto np : stack)
                    for (size_t i = np->_first; i != np->_first + np->_length; ++i)
                        curr.push_back(_word[i]);
                suffixes.push_back(std::move(curr));
            }
            else
            {
                for (size_t i = 0; i != _width; ++i)
                    if (root->_children[i] != nullptr)
                    {
                        stack.push_back(root->_children[i]);
                        _get_all_suffixes(root->_children[i], stack, suffixes);
                        stack.pop_back();
                    }
            }
        }

        std::vector<size_t> _get_all_beginning_with(_Node* curr_root, size_t to_add) const
        {
            std::vector<size_t> ret;
            ret.push_back(to_add);
            _get_all_beginning_with(curr_root, ret);
            return ret;
        }

        void _get_all_beginning_with(_Node* curr_root, std::vector<size_t>& vsz) const
        {
            if (curr_root->_children == nullptr)
                vsz.push_back(_length - curr_root->_length_from_root);
            else
                for (size_t i = 0; i != _width; ++i)
                    if (curr_root->_children[i] != nullptr)
                        _get_all_beginning_with(curr_root->_children[i], vsz);
        }

        const size_t   _length;
        const T* const _word;

        const T      _end_mark;
        const size_t _width;

        _Node* const _root;

        // 活动点
        _Node*   _active_node;                       // 活动结点
        const T* _active_edge;                       // 用字符代表的活动边
        size_t   _active_length;                      // 活动长度, 只有在活动边非null时有效
        T        _active_edge_aux;

        // 每趟循环临时变量
        size_t _curr_end;                           // 待显式插入后缀的区间, 注意这里所代表的区间形式为(begin_index, end_index)
        size_t _remainder; // 待插入后缀数量
        _Node* _previous;                               // 最后创建的节点, 后缀链接时使用
        T      _current; // 当前待插入后缀的最后一个元素

    };

    template<typename T>
    class SuffixTreeBuilder<T>::_Node
    {
        using self_type    = _Node;
        using builder_type = SuffixTreeBuilder<T>;

        friend class builder_type;

        _Node(size_t first, size_t len, size_t len_from_root, bool is_leef, size_t width, const builder_type* builder) :
            _builder(builder), _first(first), _length(len), _length_from_root(len_from_root),
            _width(width), _children(is_leef ? nullptr : new self_type * [width]()), _suffix(nullptr)
        {
        }

        void _split()
        {
            _Node* next = _children[*_builder->_active_edge];
            auto new_node = new self_type(next->_first, _builder->_active_length,
                _length_from_root + _builder->_active_length, false, _width, _builder);
            new_node->_children[_builder->_word[next->_first + _builder->_active_length]] = next;
            next->_first += _builder->_active_length;
            next->_length -= _builder->_active_length;
            _children[*_builder->_active_edge] = new_node;
        }

        const builder_type* const _builder;
        const size_t              _width;

        size_t _first;
        size_t _length;

        const size_t _length_from_root;

        _Node** _children;

        _Node* _suffix;
    };

}

