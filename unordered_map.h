#include <iostream>
#include <type_traits>
#include <vector>


template <size_t chunkSize>
class FixedAllocator {
public:
    FixedAllocator() = default;

    static FixedAllocator& get_allocator() {
        static FixedAllocator<chunkSize> new_alloc;
        return new_alloc;
    }

    int8_t* allocate() {
        if (to_ptr.empty()) {
            int8_t* ptr = new int8_t[20 * chunkSize];
            for (int i = 0; i < 20; ++i) {
                to_ptr.push_back(ptr + i * chunkSize);
            }
        }
        int8_t* ptr = to_ptr.back();
        to_ptr.pop_back();
        return ptr;
    }

    void deallocate(int8_t* ptr, size_t) {
        to_ptr.push_back(ptr);
    }
private:
    std::vector<int8_t*> to_ptr;
};

template <typename T>
class FastAllocator {
public:
    FastAllocator() = default;
    template <typename U>
    FastAllocator (const FastAllocator<U>&) {}

    using value_type = T;

    T* allocate(size_t n) {
        switch(n) {
            case 1:
                return reinterpret_cast<T*>(FixedAllocator<sizeof(T)>::get_allocator().allocate());
            case 2:
                return reinterpret_cast<T*>(FixedAllocator<2 * sizeof(T)>::get_allocator().allocate());
            case 3:
                return reinterpret_cast<T*>(FixedAllocator<3 * sizeof(T)>::get_allocator().allocate());
            case 4:
                return reinterpret_cast<T*>(FixedAllocator<4 * sizeof(T)>::get_allocator().allocate());
            case 5:
                return reinterpret_cast<T*>(FixedAllocator<5 * sizeof(T)>::get_allocator().allocate());
            default:
                return reinterpret_cast<T*>(::operator new(n * sizeof(T)));
        }
    }

    template <typename... Args>
    void construct(T* ptr, const Args&... args) {
        new(ptr) T(args...);
    }

    void destroy(T* ptr) {
        ptr->~T();
    }

    void deallocate(T* ptr, size_t n) {
        switch(n) {
            case 1:
                FixedAllocator<sizeof(T)>::get_allocator().deallocate(reinterpret_cast<int8_t*>(ptr), sizeof(T));
                break;
            case 2:
                FixedAllocator<2 * sizeof(T)>::get_allocator().deallocate(reinterpret_cast<int8_t*>(ptr), 2 * sizeof(T));
                break;
            case 3:
                FixedAllocator<3 * sizeof(T)>::get_allocator().deallocate(reinterpret_cast<int8_t*>(ptr), 3 * sizeof(T));
                break;
            case 4:
                FixedAllocator<4 * sizeof(T)>::get_allocator().deallocate(reinterpret_cast<int8_t*>(ptr), 4 * sizeof(T));
                break;
            case 5:
                FixedAllocator<5 * sizeof(T)>::get_allocator().deallocate(reinterpret_cast<int8_t*>(ptr), 5 * sizeof(T));
                break;
            default:
                ::operator delete(ptr);
        }
    }
} ;

template <typename T, typename Allocator = std::allocator<T>>
class List {
protected:
    size_t _size = 0;

    struct Node {
        Node() = default;
        Node(const T& value): value(value) {}
        Node(T&& value) noexcept: value(std::move(value)) {}
        Node(const T& value, Node* prev): value(value), prev(prev) {}
        Node(T&& value, Node* prev) noexcept: value(std::move(value)), prev(prev) {}
        T value = T();
        Node* prev = nullptr;
        Node* next = nullptr;
    };

    Node* _begin = nullptr;
    Node* _end = nullptr;

    typename std::iterator_traits<std::iterator<std::bidirectional_iterator_tag, T>>::iterator_category x;

    template<bool isConst, typename U = T>
    class common_iterator: public std::iterator<std::bidirectional_iterator_tag, U> {
    public:
        common_iterator() = default;
        common_iterator(Node* ptr): ptr(ptr) {}

        common_iterator(const common_iterator& another) {
            ptr = another.ptr;
        }

        common_iterator& operator=(const common_iterator& another) {
            ptr = another.ptr;
            return *this;
        }

        operator common_iterator<true, const T>() const {
            return common_iterator<true, const T>(ptr);
        }

        explicit operator common_iterator<false, T>() {
            return common_iterator<false, T>(ptr);
        }

        std::conditional_t<isConst, const U&, U&> operator*() const {
            return ptr->value;
        }
        std::conditional_t<isConst, const U*, U*> operator->() const {
            return &(ptr->value);
        }

        common_iterator& operator++() {
            ptr = ptr->next;
            return *this;
        }
        common_iterator operator++(int) {
            common_iterator tmp = *this;
            ptr = ptr->next;
            return tmp;
        }

        common_iterator& operator--() {
            ptr = ptr->prev;
            return *this;
        }
        common_iterator operator--(int) {
            common_iterator tmp = *this;
            ptr = ptr->prev;
            return *this;
        }

        bool operator==(const common_iterator& another) const {
            return ptr == another.ptr;
        }
        bool operator!=(const common_iterator& another) const {
            return ptr != another.ptr;
        }

        Node* ptr = nullptr;
    };

    typename std::allocator_traits<Allocator>::template rebind_alloc<Node> allocator;
public:
    using value_type = T;
    using allocator_type = typename std::allocator_traits<Allocator>::template rebind_alloc<Node>;
    using AllocTraits = std::allocator_traits<allocator_type>;
    using size_type = size_t;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = typename std::allocator_traits<Allocator>::pointer;
    using const_pointer = typename std::allocator_traits<Allocator>::const_pointer;

    using iterator = common_iterator<false>;
    using const_iterator = common_iterator<true, const T>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    explicit List() {
        _end = AllocTraits::allocate(allocator, 1);
        _begin = _end;
        _end->prev = _begin;
        _end->next = _begin;
    }

    List(size_t count, const T& value): _size(count) {
        _end = AllocTraits::allocate(allocator, 1);
        _begin = _end;
        _end->prev = _begin;
        _end->next = _begin;
        Node* prev = _begin;
        for (size_t i = 0; i < count; ++i) {
            Node* new_node = AllocTraits::allocate(allocator, 1);
            AllocTraits::construct(allocator, new_node, value);
            new_node->prev = prev;
            prev->next = new_node;
            prev = new_node;
        }
        prev->next = _end;
        _end->prev = prev;
        _begin = _end->next;
    }
    List(size_t count, T&& value): _size(count) {
        _end = AllocTraits::allocate(allocator, 1);
        _begin = _end;
        _end->prev = _begin;
        _end->next = _begin;
        Node* prev = _begin;
        for (size_t i = 0; i < count; ++i) {
            Node* new_node = AllocTraits::allocate(allocator, 1);
            AllocTraits::construct(allocator, new_node, std::move(value));
            new_node->prev = prev;
            prev->next = new_node;
            prev = new_node;
        }
        prev->next = _end;
        _end->prev = prev;
        _begin = _end->next;
    }

    List(size_t count): _size(count) {
        _end = AllocTraits::allocate(allocator, 1);
        _begin = _end;
        _end->prev = _begin;
        _end->next = _begin;
        Node* prev = _begin;
        for (size_t i = 0; i < count; ++i) {
            Node* new_node = AllocTraits::allocate(allocator, 1);
            AllocTraits::construct(allocator, new_node);
            new_node->prev = prev;
            prev->next = new_node;
            prev = new_node;
        }
        prev->next = _end;
        _end->prev = prev;
        _begin = _end->next;
    }

    List(const List& another): allocator(AllocTraits::select_on_container_copy_construction(another.allocator)) {
        _end = AllocTraits::allocate(allocator, 1);
        _begin = _end;
        _end->prev = _begin;
        _end->next = _begin;
        for (const_iterator it = another.begin(); it != another.end(); ++it) {
            push_back(*it);
        }
    }

    List(List&& another) noexcept {
        if constexpr (std::is_base_of_v<std::true_type, typename List<T, Allocator>::AllocTraits::propagate_on_container_move_assignment>) {
            allocator = another.allocator;
        } else {
            allocator = Allocator();
        }
        _size = another._size;
        _end = another._end;
        _begin = another._begin;
        another._end = AllocTraits::allocate(allocator, 1);
        another._begin = another._end;
        another._end->next = another._begin;
        another._begin->prev = another._end;
        another._size = 0;
    }

    List& operator=(const List& another) {
        if (this == &another) return *this;
        iterator it = begin();
        while ( it != end()) {
            erase(it++);
        }
        AllocTraits::deallocate(allocator, _end, 1);
        if constexpr (std::is_base_of_v<std::true_type, typename List<T, Allocator>::AllocTraits::propagate_on_container_copy_assignment>) {
            allocator = another.allocator;
        }
        _end = AllocTraits::allocate(allocator, 1);
        _begin = _end;
        _end->prev = _begin;
        _end->next = _begin;
        _size = 0;
        for (const_iterator i = another.begin(); i != another.end(); ++i) {
            push_back(*i);
        }
        return *this;
    }

    List& operator=(List&& another) noexcept {
        if (this == &another) return *this;
        iterator it = begin();
        while ( it != end()) {
            erase(it++);
        }
        AllocTraits::deallocate(allocator, _end, 1);
        if constexpr (std::is_base_of_v<std::true_type, typename List<T, Allocator>::AllocTraits::propagate_on_container_move_assignment>) {
            allocator = another.allocator;
        }
        _size = another._size;
        _end = another._end;
        _begin = another._begin;
        another._end = AllocTraits::allocate(allocator, 1);
        another._begin = another._end;
        another._end->next = another._begin;
        another._begin->prev = another._end;
        another._size = 0;
        return *this;
    }

    iterator insert(const_iterator iter, const T& value) {
        Node* tmp = AllocTraits::allocate(allocator, 1);
        AllocTraits::construct(allocator, tmp, value);
        tmp->next = iter.ptr;
        tmp->prev = iter.ptr->prev;
        tmp->prev->next = tmp;
        iter.ptr->prev = tmp;
        ++_size;
        return iterator(tmp->prev);
    }
    iterator insert(const_iterator iter, T&& value) {
        Node* tmp = AllocTraits::allocate(allocator, 1);
        AllocTraits::construct(allocator, tmp, std::move(value));
        tmp->next = iter.ptr;
        tmp->prev = iter.ptr->prev;
        tmp->prev->next = tmp;
        iter.ptr->prev = tmp;
        ++_size;
        return iterator(tmp->prev);
    }

    template <typename... Args>
    iterator emplace(const_iterator iter, Args&&... args) {
        Node* tmp = AllocTraits::allocate(allocator, 1);
        AllocTraits::construct(allocator, tmp, std::forward<Args>(args)...);
        tmp->next = iter.ptr;
        tmp->prev = iter.ptr->prev;
        tmp->prev->next = tmp;
        iter.ptr->prev = tmp;
        ++_size;
        return iterator(tmp->prev);
    }

    void erase(const_iterator iter) {
        Node* to_be_begin = _begin;
        if (iter == cbegin()) {
            to_be_begin = iter.ptr->next;
        }
        iter.ptr->prev->next = iter.ptr->next;
        iter.ptr->next->prev = iter.ptr->prev;
        AllocTraits::destroy(allocator, iter.ptr);
        AllocTraits::deallocate(allocator, iter.ptr, 1);
        --_size;
        _begin = to_be_begin;
    }

    void push_back(const T& value) {
        insert(_end, value);
        if (_size == 1) {
            _begin = _end->prev;
        }
    }
    void push_back(T&& value) {
        insert(_end, std::move(value));
        if (_size == 1) {
            _begin = _end->prev;
        }
    }

    void push_front(const T& value) {
        insert(_begin, value);
        _begin = iterator(_begin->prev).ptr;
    }
    void push_front(T&& value) {
        insert(_begin, std::move(value));
        _begin = iterator(_begin->prev).ptr;
    }

    void pop_back() {
        iterator to_erase = _end;
        --to_erase;
        erase(to_erase);
    }
    void pop_front() {
        iterator to_erase = _begin;
        _begin = _begin->next;
        erase(to_erase);
    }

    iterator begin() {
        return _begin;
    }
    const_iterator begin() const {
        return _begin;
    }
    iterator end() {
        return _end;
    }
    const_iterator end() const {
        return _end;
    }
    const_iterator cbegin() const {
        return _begin;
    }
    const_iterator cend() const {
        return _end;
    }

    reverse_iterator rbegin() {
        return reverse_iterator(_end);
    }
    const_reverse_iterator rbegin() const {
        return const_reverse_iterator(_end);
    }
    reverse_iterator rend() {
        return reverse_iterator(_begin);
    }
    const_reverse_iterator rend() const {
        return const_reverse_iterator(_begin);
    }
    const_reverse_iterator crbegin() const {
        return reverse_iterator(_end);
    }
    const_reverse_iterator crend() const {
        return reverse_iterator(_begin);
    }

    typename std::allocator_traits<Allocator>::template rebind_alloc<Node> get_allocator() {
        return allocator;
    }

    size_t size() const {
        return _size;
    }

    ~List() {
        iterator it = begin();
        while ( it != end()) {
            erase(it++);
        }
        AllocTraits::deallocate(allocator, _end, 1);
    }
};

template <
        typename Key,
        typename Value,
        typename Hash=std::hash<Key>,
        typename Equal=std::equal_to<Key>,
        typename Alloc=std::allocator<std::pair<const Key, Value>>
        >
class UnorderedMap {
public:
    using NodeType = std::pair<const Key, Value>;
    using pointer = typename std::allocator_traits<Alloc>::pointer;
    using const_pointer = typename std::allocator_traits<Alloc>::const_pointer;

    UnorderedMap(int64_t size = 100): _bucket_count(size), _list(), _buckets(_bucket_count, _list.end()) {}
    UnorderedMap(const UnorderedMap& another):  _bucket_count(another._bucket_count), _list(), _buckets(_bucket_count, _list.end()),
        _allocator(AllocTraits::select_on_container_copy_construction(another._allocator)), _NTallocator(NTAllocTraits::select_on_container_copy_construction(another._NTallocator)) {
        for (auto it: another._list) {
            insert(it);
        }
    }
    UnorderedMap(UnorderedMap&& another) noexcept: _bucket_count(another._bucket_count), _list(std::move(another._list)),
        _buckets(std::move(another._buckets)) {
        if constexpr (std::is_base_of_v<std::true_type, typename AllocTraits::propagate_on_container_move_assignment>) {
            _allocator = another._allocator;
        } else {
            _allocator = allocator_type();
        }
        if constexpr (std::is_base_of_v<std::true_type, typename NTAllocTraits::propagate_on_container_move_assignment>) {
            _NTallocator = another._NTallocator;
        } else {
            _NTallocator = Alloc();
        }
        another._bucket_count = 0;
    }
    UnorderedMap(int64_t size, double max_l_factor): UnorderedMap(size) {
        _max_load_factor = max_l_factor;
    }

    UnorderedMap& operator=(const UnorderedMap& another) {
        if (this == &another)
            return *this;
        _bucket_count = another._bucket_count;
        _list = another._list;
        _buckets = another._buckets;
        if constexpr (std::is_base_of_v<std::true_type, typename AllocTraits::propagate_on_container_copy_assignment>) {
            _allocator = another._allocator;
        }
        if constexpr (std::is_base_of_v<std::true_type, typename NTAllocTraits::propagate_on_container_copy_assignment>) {
            _NTallocator = another._NTallocator;
        }
        return *this;
    }

    UnorderedMap& operator=(UnorderedMap&& another) noexcept {
        if (this == &another)
            return *this;
        _bucket_count = another._bucket_count;
        another._bucket_count = 0;
        _list = std::move(another._list);
        _buckets = std::move(another._buckets);
        if constexpr (std::is_base_of_v<std::true_type, typename AllocTraits::propagate_on_container_move_assignment>) {
            _allocator = another._allocator;
        }
        if constexpr (std::is_base_of_v<std::true_type, typename NTAllocTraits::propagate_on_container_copy_assignment>) {
            _NTallocator = another._NTallocator;
        }
        return *this;
    }

private:
    struct Node {
        Node(const NodeType& x, uint64_t cached): data(x), cached(cached) {}
        Node(NodeType&& x, uint64_t cached) noexcept: data(std::make_pair(const_cast<Key&&>(std::move(x.first)), std::move(x.second))), cached(cached) {}
        Node(NodeType&& x) noexcept: data(std::make_pair(const_cast<Key&&>(std::move(x.first)), std::move(x.second))) {}
        Node(const Node& x): data(x.data), cached(x.cached) {}
        Node(Node&& x) noexcept: data(std::make_pair(const_cast<Key&&>(std::move(x.data.first)), std::move(x.data.second))), cached(x.cached) {
            x.cached = 0;
        }
        NodeType data;
        uint64_t cached;
    };
    using list_iterator = typename List<Node, typename std::allocator_traits<Alloc>::template rebind_alloc<Node>>::iterator;
    using list_const_iterator = typename List<Node, typename std::allocator_traits<Alloc>::template rebind_alloc<Node>>::const_iterator;
    using allocator_type = typename std::allocator_traits<Alloc>::template rebind_alloc<Node>;
    using AllocTraits = std::allocator_traits<allocator_type>;
    using NTAllocTraits = std::allocator_traits<Alloc>;
    using list_type = List<Node, allocator_type>;
    using buckets_type = std::vector<list_iterator>;
    int64_t _bucket_count;
    list_type _list;
    buckets_type _buckets;
    allocator_type _allocator;
    Alloc _NTallocator;
    double _max_load_factor = 0.9;

public:
    class const_iterator: public list_const_iterator {
    public:
        const_iterator(const list_const_iterator& another): list_const_iterator(another) {}
        const NodeType& operator*() const {
            return list_const_iterator::ptr->value.data;
        }
        const NodeType* operator->() const {
            return &(list_const_iterator::ptr->value.data);
        }

    };
    class iterator: public list_iterator {
    public:
        iterator(const list_iterator& another): list_iterator(another) {}
        NodeType& operator*() {
            return list_iterator::ptr->value.data;
        }
        NodeType* operator->() {
            return &(list_iterator::ptr->value.data);
        }
        operator const_iterator() const {
            return const_iterator(list_iterator::ptr);
        }
    };
    using Iterator = iterator;
    using ConstIterator = const_iterator;
    void rehash(int64_t new_size) {
        UnorderedMap<Key, Value, Hash, Equal, Alloc> new_map(new_size, _max_load_factor);
        for (auto it = _list.begin(); it != _list.end(); ++it) {
            new_map.insert(std::move(*it));
        }
        std::swap(*this, new_map);
    }
    double load_factor() const {
        return double(size()) / double(_buckets.size());
    }
    double max_load_factor() const {
        return _max_load_factor;
    }
    double& max_load_factor() {
        return _max_load_factor;
    }
    void max_load_factor(double new_l_factor) {
        _max_load_factor = new_l_factor;
    }
    int64_t max_size() const {
        return int64_t(_max_load_factor * _buckets.size());
    }
    void reserve(int64_t max_count) {
        rehash(max_count);
    }

    const_iterator find(const Key& key, bool = false) const {
        uint64_t current_hash = Hash()(key) % _buckets.size();
        list_const_iterator it = my_special_find(key, current_hash);
        if (it == _list.cend() || it->cached != current_hash)
            return _list.end();
        else
            return it;
    }

    iterator find(const Key& key) {
        return static_cast<list_iterator>(find(key, false));
    }

    std::pair<iterator, bool> insert(const NodeType& value) {
        if (max_size() < size())
            rehash(_buckets.size());
        uint64_t current_hash = Hash()(value.first) % _buckets.size();
        return insert(Node(value, current_hash));
    }
    std::pair<iterator, bool> insert(NodeType&& value) {
        if (max_size() < size())
            rehash(_buckets.size());
        uint64_t current_hash = Hash()(value.first) % _buckets.size();
        return insert(std::move(Node(std::move(value), current_hash)));
    }
    template <typename InputIterator>
    void insert(const InputIterator& first, const InputIterator& second) {
        for (auto it = first; it != second; ++it) {
            insert(*it);
        }
    }

    template <typename... Args>
    std::pair<iterator, bool> emplace(Args&&... args) {
        Node* x = AllocTraits::allocate(_allocator, 1);
        std::allocator_traits<Alloc>::construct(_NTallocator, &(x->data), std::forward<Args>(args)...);
        x->cached = Hash()(x->data.first) % _buckets.size();
        auto to_return = insert(std::move(*x));
        AllocTraits::destroy(_allocator, x);
        AllocTraits::deallocate(_allocator, x, 1);
        return to_return;
    }

    void erase(const list_iterator& it) {
        list_iterator next = it;
        ++next;
        if (_buckets[it->cached] == static_cast<list_iterator>(it))
            _buckets[it->cached] = next;
        if (next == _list.end() || next->cached != it->cached)
            _buckets[it->cached] = _list.end();
        _list.erase(it);
    }
    void erase(const list_iterator& first, const list_iterator& second) {
        for (auto it = first; it != second; ++it) {
            erase(it);
        }
    }

    Value& operator[](const Key& key) {
        uint64_t current_hash = Hash()(key) % _buckets.size();
        list_iterator it = my_special_find(key, static_cast<int64_t>(current_hash));
        if (it == _list.end() || it->cached != current_hash) {
            if (it == _list.end()) {
                _list.push_back(Node(std::make_pair(key, Value()), current_hash));
                _buckets[current_hash] = --_list.end();
            } else if (it == _list.begin()) {
                _list.push_front(Node(std::make_pair(key, Value()), current_hash));
                _buckets[current_hash] = _list.begin();
            } else
                _buckets[current_hash] = _list.insert(it, Node(std::make_pair(key, Value()), current_hash));
            return iterator(_buckets[current_hash])->second;
        } else
            return iterator(it)->second;
    }
    Value& at(const Key& key) {
        uint64_t current_hash = Hash()(key) % _buckets.size();
        list_iterator it = my_special_find(key, current_hash);
        if (it == _list.end() || it->cached != current_hash) {
            throw std::out_of_range("");
        } else
            return iterator(it)->second;
    }


    iterator begin() {
        return _list.begin();
    }
    const_iterator begin() const {
        return _list.begin();
    }
    const_iterator cbegin() const {
        return _list.cbegin();
    }
    iterator end() {
        return _list.end();
    }
    const_iterator end() const {
        return _list.end();
    }
    const_iterator cend() const {
        return _list.cend();
    }

    int64_t size() const {
        return _list.size();
    }

private:
    std::pair<iterator, bool> insert(const Node& value) {
        uint64_t current_hash = value.cached;
        list_iterator it = my_special_find(value.data.first, current_hash);
        if (it == _list.end() || it->cached != current_hash)
            it = _buckets[current_hash];
        else
            return {it, false};

        if (it == _list.end()) {
            _list.push_back(value);
            _buckets[current_hash] = --_list.end();
        } else if (it == _list.begin()) {
            _list.push_front(value);
            _buckets[current_hash] = _list.begin();
        } else
            _buckets[current_hash] = _list.insert(it, value);
        return {_buckets[current_hash], true};
    }
    std::pair<iterator, bool> insert(Node&& value) {
        uint64_t current_hash = value.cached;
        list_iterator it = my_special_find(value.data.first, current_hash);
        if (it == _list.end() || it->cached != current_hash)
            it = _buckets[current_hash];
        else
            return {it, false};

        if (it == _list.end()) {
            _list.push_back(std::move(value));
            _buckets[current_hash] = --_list.end();
        } else if (it == _list.begin()) {
            _list.push_front(std::move(value));
            _buckets[current_hash] = _list.begin();
        } else
            _buckets[current_hash] = _list.insert(it, std::move(value));
        return {_buckets[current_hash], true};
    }

    list_iterator my_special_find(const Key& key, uint64_t current_hash) const {
        list_const_iterator it = _buckets[current_hash];
        while (it != _list.end() && it->cached == current_hash && !Equal()(it->data.first, key)) {
            ++it;
        }
        return static_cast<list_iterator>(it);
    }
};
