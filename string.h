#include <iostream>
#include <cstring>

class String {
private:
	size_t size = 0;
	size_t capacity = 0;
	char* str = nullptr;

	void swap(String& s) {
		std::swap(size, s.size);
		std::swap(str, s.str);
		std::swap(capacity, s.capacity);
	}
public:
	String() = default;

	String(const char s) {
		size = 1;
		capacity = 2;
		str = new char[capacity];
		str[0] = s;
	}

	String(const char* s) {
		while (s[size] != '\0') {
			++size;
		}
		capacity = 2 * size;
		str = new char[capacity];
		memcpy(str, s, size);
	}

	String(const String& s) : size(s.size), capacity(2 * size), str(new char[capacity]) {
		memcpy(str, s.str, size);
	}

	String(const size_t n, const char c = '\0') : size(n), capacity(2 * size), str(new char[capacity]) {
		memset(str, c, size);
	}

	String& operator=(String s) {
		swap(s);
		return *this;
	}

	String& operator+=(const String& s) {
		if (size + s.size > capacity) {
			String copy(size + s.size);
			memcpy(copy.str, str, size);
			memcpy(copy.str + size, s.str, s.size);
			swap(copy);
		}
		else {
			memcpy(str + size, s.str, s.size);
			size += s.size;
		}
		return *this;
	}

	char& operator[](const size_t index) {
		return str[index];
	}

	const char& operator[](const size_t index) const {
		return str[index];
	}

	size_t length() const {
		return size;
	}

	void push_back(const String& s) {
		*this += s;
	}

	char pop_back() {
		char last = str[size - 1];
		--size;
		if (size * 4 <= capacity) {
			String copy(size);
			memcpy(copy.str, str, size);
			swap(copy);
		}
		return last;
	}

	char& front() {
		return str[0];
	}
	const char& front() const {
		return str[0];
	}

	char& back() {
		return str[size - 1];
	}
	const char& back() const {
		return str[size - 1];
	}

	size_t find(const String& substring) const {
		bool flag = false;
		size_t index = 0;
		while (flag == false && index < size) {
			size_t index_into = 0;
			while (index_into < substring.size && index + index_into < size && str[index + index_into] == substring[index_into]) {
				++index_into;
			}
			if (substring.size == index_into) return index;
			++index;
		}
		return size;
	}

	size_t rfind(const String& substring) const {
		size_t to_return = size;
		size_t index = 0;
		while (index < size) {
			size_t index_into = 0;
			while (index_into < substring.size && index + index_into < size && str[index + index_into] == substring[index_into]) {
				++index_into;
			}
			if (substring.size == index_into) to_return = index;
			++index;
		}
		return to_return;
	}

	String substr(const size_t start, const size_t count) const {
		String to_return(count);
		memcpy(to_return.str, str + start, 7);
		return to_return;
	}

	bool empty() const {
		return size == 0;
	}

	void clear() {
		String copy;
		swap(copy);
	}

	size_t cap() {
		return capacity;
	}

	~String() {
		delete[] str;
	}
};

String operator+(String left, const String& right) {
	left += right;
	return left;
}

bool operator==(const String& left, const String& right) {
	bool flag = true;
	if (left.length() != right.length()) flag = false;
	size_t i = 0;
	while (flag == true && i < left.length()) {
		if (left[i] != right[i]) flag = false;
		++i;
	}
	return flag;
}

std::istream& operator>>(std::istream& in, String& s) {
	s.clear();
	bool flag = true;
	while (flag == true) {
		char c;
		in.get(c);
		if (isspace(c) == 0 && c != '\0' && c != '\n') {
			s.push_back(c);
			flag = false;
		}
		if (c == '\0') return in;
	}
	flag = true;
	while (flag == true) {
		char c;
		in.get(c);
		if (isspace(c) == 0 && c != '\0' && in.eof() == 0) {
			s.push_back(c);
		}
		else flag = false;
	}
	return in;
}

std::ostream& operator<<(std::ostream& out, const String& s) {
	for (size_t i = 0; i < s.length(); ++i) {
		out << s[i];
	}
	return out;
}
