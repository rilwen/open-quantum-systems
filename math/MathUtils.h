#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <istream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <Eigen/Core>
#include "MathCore.h"


namespace rql { namespace math {
	// Round function (round to nearest)
	inline double round(double x)
	{
		return floor(x + 0.5);
	}

	// Round function (truncate towards zero)
	inline double trunc(double x)
	{
		return x >= 0 ? std::floor(x) : std::ceil(x);
	}

	// Greatest common denominator
	long gcd(long a, long b);

	// Convert type to a string.
	template <typename T>
	std::string stringify(const T& i)
	{
		std::ostringstream o;
		o << i;
		return o.str();
	}

	// Chech if value of type T is finite.
	template <typename T>
	bool finite(const T& val)
	{
		bool finite = true;
		finite &= !(val != val); // check for NaN
		finite &= (val != std::numeric_limits<T>::infinity()); // check for +Inf
		finite &= (val != -std::numeric_limits<T>::infinity()); // check for -Inf
		return finite;
	}

	// value of pi
	RQL_MATH_API extern const double PI;

	template <class T> class eigen_aligned_data
	{
	public:
		explicit eigen_aligned_data(size_t size);		
		eigen_aligned_data(eigen_aligned_data<T>& other);
		eigen_aligned_data(eigen_aligned_data<T>&& other);
		eigen_aligned_data<T>& operator=(eigen_aligned_data<T>& other);		
		~eigen_aligned_data();
		size_t size() const { return m_size; }
		T* get() { return m_data; }
		const T* get() const { return m_data; }
		typedef T* iterator;
		typedef const T* const_iterator;
		iterator begin() { return m_data; }
		const_iterator begin() const { return m_data; }
		iterator end() { return m_end; }
		const_iterator end() const { return m_end; }
	private:
		void wipe();
		void clear();
	private:		
//		eigen_aligned_data(const eigen_aligned_data<T>& other); // not implemented		
//		eigen_aligned_data<T>& operator=(const eigen_aligned_data<T>& other); // not implemented
		T* m_data;
		size_t m_size;
		T* m_end;
	};

	

	template <class T> eigen_aligned_data<T>::eigen_aligned_data(size_t size)
		: m_data(static_cast<T*>(Eigen::internal::aligned_malloc(sizeof(T)*size))), m_size(size), m_end(m_data + m_size)
	{	
		if (!m_data)
			throw std::bad_alloc();
	}

	template <class T> eigen_aligned_data<T>::eigen_aligned_data(eigen_aligned_data<T>& other)
		: m_data(other.m_data), m_size(other.m_size), m_end(other.m_end)
	{
		other.clear();
	}
	
	template <class T> eigen_aligned_data<T>::eigen_aligned_data(eigen_aligned_data<T>&& other)
		: m_data(other.m_data), m_size(other.m_size), m_end(other.m_end)
	{
		other.clear();
	}

	template <class T> eigen_aligned_data<T>& eigen_aligned_data<T>::operator=(eigen_aligned_data<T>& other)
	{
		wipe();
		m_data = other.m_data;
		m_size = other.m_size;
		m_end = other.m_end;
		other.clear();
		return *this;
	}

	template <class T> eigen_aligned_data<T>::~eigen_aligned_data()
	{
		wipe();
	}

	template <class T> void eigen_aligned_data<T>::wipe()
	{
		if (m_data) {
			Eigen::internal::aligned_free(m_data);
			clear();
		}
	}

	template <class T> void eigen_aligned_data<T>::clear()
	{
		m_data = 0;
		m_size = 0;
		m_end = 0;
	}

	template <class Value> 
	size_t load_matrix_from_stream(std::istream& input, std::vector<std::vector<Value> >& matrix)
	{
		size_t max_cols = 0;
		std::string line;
		std::vector<std::string> tokens;
		while (input.good()) {
			getline(input, line);
			boost::split(tokens, line,	boost::is_any_of(" \t"));
			const size_t n = tokens.size();			
			if (n) {
				std::vector<Value> current_row;
				current_row.reserve(n);
				const std::vector<std::string>::const_iterator end = tokens.end();
				unsigned int valid_cols = 0;
				for (std::vector<std::string>::const_iterator it = tokens.begin(); it != end; ++it) {
					Value v;
					try {
						v = boost::lexical_cast<Value>(*it);
					} catch (boost::bad_lexical_cast&) {
						continue;
					}
					current_row.push_back(v);
					++valid_cols;
				}
				max_cols = std::max(max_cols, static_cast<size_t>(valid_cols));
				if (valid_cols) {
					matrix.push_back(current_row);
				}
			}
		}
		return max_cols;
	}	

	template <class Value, class Matrix> 
	void load_matrix_from_stream(std::istream& input, Matrix& matrix)
	{
		std::vector<std::vector<Value> > data;
		const size_t ncols = load_matrix_from_stream(input, data);
		const size_t nrows = data.size();		
		matrix.setZero(nrows, ncols);
		const typename std::vector<std::vector<Value> >::const_iterator rend = data.end();
		size_t row = 0;
		for (typename std::vector<std::vector<Value> >::const_iterator r = data.begin(); r != rend; ++r) {
			const typename std::vector<Value>::const_iterator cend = r->end();
			size_t col = 0;
			for (typename std::vector<Value>::const_iterator c = r->begin(); c != cend; ++c) {
				matrix(row, col) = *c;
				++col;
			}
			++row;
		}
	}
}}
