#ifndef UTIL_H_
#define UTIL_H_

#include <boost/functional/hash.hpp>
#include <string>
#include <vector>
#include <math.h>
#include <pthread.h>

#define BINARY_THRESHOLD 4

using namespace std;

typedef set<pair<int, bool>> star; //el, dir
typedef set<pair<int, int>> match_star; //uid, vid

enum RBI { NONE, RED, IVORY, BLACK, FIXED };

const int MAX_INT = numeric_limits<int>::max();
const int MIN_INT = numeric_limits<int>::min();
const size_t MAX_SIZE = numeric_limits<size_t>::max();

vector<string> parse(string line, string del);

inline int search(const int* arr, int begin, int end, int target) {
	if (end - begin > BINARY_THRESHOLD) {
		int mid;
		while (begin < end) {
			mid = (begin + end) / 2;
			if (arr[mid] < target) begin = mid + 1;
			else if (arr[mid] > target) end = mid;
			else return mid;
		}
	} else {
		for (int i = begin; i < end; i++) {
            if (target == arr[i])
                return i;
			else if (target < arr[i])
				return -1;
		}
	}
	return -1;
}

inline int search_geq(const int* arr, int begin, int end, int target) {
	if (end - begin > BINARY_THRESHOLD) {
		int mid;
		while (begin != end) {
			mid = (begin + end) / 2;
			if (arr[mid] < target) begin = mid + 1;
			else end = mid;
		}
		return begin;
	} else {
		for (int i = begin; i < end; i++) {
            if (target <= arr[i])
                return i;
		}
		return end;
	}
}


vector<string> tokenize(const string& line, const char* delim);

struct range {
	const int* begin;
	const int* end;

	range(const int* b = NULL, const int* e = NULL) : begin(b), end(e) {}

	bool operator< (const range& other) const {
		int size = this->end - this->begin;
		int osize = other.end - other.begin;
		if (size == osize)
			return this->begin < other.begin;
		else
			return size < osize;
	}

	bool operator== (const range& other) const {
		return this->begin == other.begin && this->end == other.end; 
	}

	bool good() {
		for (int i = 0; i < (end-begin)-1; i++)
			if (begin[i] >= begin[i+1])
				return false;
		return true;
	}

	int size() {
		return end - begin;
	}
};

struct AdjElem {
	int id;
	int label;
	bool dir;

	AdjElem() : id(-1), label(-1), dir(true) {} 
	AdjElem(int _id, int _label, bool _dir) : id(_id), label(_label), dir(_dir) {}
	AdjElem(const AdjElem& other) : id(other.id), label(other.label), dir(other.dir) {}
};

struct meta_range {
	const int* begin;
	const int* end;
	AdjElem edge;

	meta_range(const range& r) : begin(r.begin), end(r.end) {}
	meta_range(const range& r, int id, int label, bool dir) : begin(r.begin), end(r.end), edge(id, label, dir) {}

	bool operator< (const meta_range& other) const {
		int size = this->end - this->begin;
		int osize = other.end - other.begin;
		if (size == osize)
			return this->begin < other.begin;
		else
			return size < osize;
	}

	bool operator== (const meta_range& other) const {
		return this->begin == other.begin && this->end == other.end; 
	}
};

struct Edge {
	int src, dst, el;

	Edge(int s, int d, int e) : src(s), dst(d), el(e) {}
	Edge() {}

	bool operator<(const Edge& other) const {
		if (src < other.src)
			return true;
		else if (src == other.src) {
			if (el < other.el)
				return true;
			else if (el == other.el) {
				return dst < other.dst;
			}
		}
		return false;
	}
	bool operator==(const Edge& other) const {
		return src == other.src && dst == other.dst && el == other.el;
	}
};

struct EdgeHasher {
	std::size_t operator () (const Edge &key) const 
	{
		// The following line is a stright forward implementation. But it can be
		// hard to design a good hash function if KeyData is complex.
		//return (key.id << 32 | key.age); // suppose size_t is 64-bit and int is 32-bit
		// A commonly used way is to use boost
		std::size_t seed = 0;
		boost::hash_combine(seed, boost::hash_value(key.src));
		boost::hash_combine(seed, boost::hash_value(key.dst));
		boost::hash_combine(seed, boost::hash_value(key.el));
		return seed;
	}
};

const size_t MIN_CAP = 4096;

template<class T>
class LargeVector {
private:
	size_t m_capacity;
	size_t m_size;
	T*     m_data;
	
public:
	LargeVector(size_t s = 0) {
		m_capacity = std::max(s, MIN_CAP);
		m_size = s;
		m_data = new T[m_capacity];
	}

	inline size_t size() const {
		return m_size;
	}

	inline size_t capacity() const {
		return m_capacity;
	}

	inline T* data() {
		return m_data;
	}

	inline T* begin() {
		return m_data;
	}

	inline T* end() {
		return m_data + m_size;
	}

	inline void resize(size_t new_s) {
		if (m_capacity < new_s) {
			m_capacity = std::max(2 * m_capacity, new_s);
			T* temp = new T[m_capacity];
			memcpy(temp, m_data, sizeof(T) * m_size);
			delete [] m_data;
			m_data = temp;
		}
		m_size = new_s;
	}

	inline void reserve(size_t new_c) {
		if (m_capacity < new_c) {
			m_capacity = std::max(2 * m_capacity, new_c);
			T* temp = new T[m_capacity];
			memcpy(temp, m_data, sizeof(T) * m_size);
			delete [] m_data;
			m_data = temp;
		}
	}

	inline T& operator[](size_t i) {
		return m_data[i];
	}

	inline T& operator[](size_t i) const {
		return m_data[i];
	}

	inline void push_back(const T& t) {
		if (m_capacity == m_size) {
			m_capacity *= 2;
			T* temp = new T[m_capacity];
			memcpy(temp, m_data, sizeof(T) * m_size);
			delete [] m_data;
			m_data = temp;
		}
		m_data[m_size++] = t;
	}

	inline void erase(T* begin) {
		m_size = begin - m_data;
	}
};

struct Bitmap {
	vector<char> bits;
	Bitmap(int size = 0) : bits(size) {} 

	inline void resize(int size) {
		bits.resize(size, 0);
	}

	inline void clear() {
		bits.clear();
	}

	inline bool get(const int& a, const int& b) {
		return bits[a] & (1 << b);
	}

	inline void set(const int& a, const int& b) {
		bits[a] |= (1 << b);
	}

	inline void unset(const int& a, const int& b) {
		bits[a] ^= (1 << b);
	}

	inline void unsetAll(const int& a) {
		bits[a] = 0;
	}
};


#endif
