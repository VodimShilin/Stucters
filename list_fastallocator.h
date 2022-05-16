#include <iostream>
#include <type_traits>
#include <list>

template <typename T>
class C {
public:
    C(T) = delete;
};

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
private:
    size_t _size = 0;

    struct Node {
        Node() = default;
        Node(const T& value): value(value) {}
        Node(const T& value, Node* prev): value(value), prev(prev) {}
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

        operator common_iterator<true, const T>() {
            return common_iterator<true, const T>(ptr);
        }

        std::conditional_t<isConst, const U&, U&> operator*() {
            return ptr->value;
        }
        std::conditional_t<isConst, const U*, U*> operator->() {
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
    //Allocator exact_allocator;
public:
    using value_type = T;
    using allocator_type = std::allocator_traits<typename std::allocator_traits<Allocator>::template rebind_alloc<Node>>;
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
        _end = allocator_type::allocate(allocator, 1);
//        _end = allocator.allocate(1);
        _begin = _end;
        _end->prev = _begin;
        _end->next = _begin;
    }

    List(size_t count, const T& value): _size(count) {
        _end = allocator.allocate(1);
        _begin = _end;
        _end->prev = _begin;
        _end->next = _begin;
        Node* prev = _begin;
        for (size_t i = 0; i < count; ++i) {
            Node* new_node = allocator.allocate(1);
            allocator.construct(new_node, value);
            new_node->prev = prev;
            prev->next = new_node;
            prev = new_node;
        }
        prev->next = _end;
        _end->prev = prev;
        _begin = _end->next;
    }

    List(size_t count): _size(count) {
        _end = allocator.allocate(1);
        _begin = _end;
        _end->prev = _begin;
        _end->next = _begin;
        Node* prev = _begin;
        for (size_t i = 0; i < count; ++i) {
            Node* new_node = allocator.allocate(1);
            allocator.construct(new_node);
            new_node->prev = prev;
            prev->next = new_node;
            prev = new_node;
        }
        prev->next = _end;
        _end->prev = prev;
        _begin = _end->next;
    }

    List(const List& another) {
        allocator = allocator_type::select_on_container_copy_construction(another.allocator);
        _end = allocator.allocate(1);
        _begin = _end;
        _end->prev = _begin;
        _end->next = _begin;
        for (const_iterator it = another.begin(); it != another.end(); ++it) {
            push_back(*it);
        }
    }

    List& operator=(const List& another) {
        if (this == &another) return *this;
        iterator it = begin();
        while ( it != end()) {
//            iterator new_it = it;
//            ++new_it;
            erase(it++);
//            it = new_it;
        }
        allocator.deallocate(_end, 1);
        if (std::is_base_of_v<std::true_type, typename List<T, Allocator>::allocator_type::propagate_on_container_copy_assignment>) {
            allocator = another.allocator;
        }
        _end = allocator.allocate(1);
        _begin = _end;
        _end->prev = _begin;
        _end->next = _begin;
        _size = 0;
        for (const_iterator i = another.begin(); i != another.end(); ++i) {
            push_back(*i);
        }
        return *this;
    }

    void insert(const_iterator iter, const T& value) {
        Node* tmp = allocator.allocate(1);
        allocator.construct(tmp, value);
        //tmp->value = value;
        tmp->next = iter.ptr;
        tmp->prev = iter.ptr->prev;
        tmp->prev->next = tmp;
        iter.ptr->prev = tmp;
        ++_size;
    }

    void erase(const_iterator iter) {
        Node* to_be_begin = _begin;
        if (iter == cbegin()) {
            to_be_begin = iter.ptr->next;
        }
        iter.ptr->prev->next = iter.ptr->next;
        iter.ptr->next->prev = iter.ptr->prev;
        allocator.destroy(iter.ptr);
        allocator.deallocate(iter.ptr, 1);
        --_size;
        _begin = to_be_begin;
    }

    void push_back(const T& value) {
        insert(_end, value);
        if (_size == 1) {
            _begin = _end->prev;
        }
    }
    void push_front(const T& value) {
        insert(_begin, value);
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
//            iterator new_it = it;
//            ++new_it;
            erase(it++);
//            it = new_it;
        }
        allocator.deallocate(_end, 1);
    }
};
