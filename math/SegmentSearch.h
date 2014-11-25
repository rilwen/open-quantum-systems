#ifndef __MATH_SEGMENT_SEARCH_H
#define __MATH_SEGMENT_SEARCH_H

#include <stdexcept>

namespace rql {
	struct SegmentSearch
	{
		//! Find such i that data[i] <= x < data[i+1], using binary search
		//! @param data Ordered vector of values
		//! @param size Size of data, non-zero
		//! @tparam V Vector type with operator[]
		//! @return Index i if a suitable segment is found OR -1 if x < data[0]
		template <class V, class T>
		static int binary_search_left_inclusive(const V& data, unsigned int size, const T& x);

		//! Find such i that data[i] <= x < data[i+1], using binary search
		//! @param data Ordered vector of values
		//! @tparam V Vector type with operator[]
		//! @return Index i if a suitable segment is found OR -1 if x < data[0]
		template <class V, class T>
		static int binary_search_left_inclusive(const V& data, const T& x)
		{
			return binary_search_left_inclusive(data, data.size(), x);
		}

		//! Find such i that data[i] < x <= data[i+1], using binary search
		//! @param data Ordered vector of values
		//! @param size Size of data, non-zero.
		//! @tparam V Vector type with operator[] and .size() method
		//! @return Index i if a suitable segment is found OR -1 if x <= data[0]
		template <class V, class T>
		static int binary_search_right_inclusive(const V& data, unsigned int size, const T& x);

		//! Find such i that data[i] < x <= data[i+1], using binary search
		//! @param data Ordered vector of values
		//! @tparam V Vector type with operator[] and .size() method
		//! @return Index i if a suitable segment is found OR -1 if x <= data[0]
		template <class V, class T>
		static int binary_search_right_inclusive(const V& data, const T& x)
		{
			return binary_search_right_inclusive(data, data.size(), x);
		}
	};

	template <class V, class T> int SegmentSearch::binary_search_left_inclusive(const V& data, unsigned int size, const T& x)
	{
		if (size == 0)
			throw std::domain_error("Zero-sized array");
		if (x < data[0])
			return -1;

		// there is a segment which contains x

		// logically, we assume that data[size] == +Infinity
		unsigned int l = 0; // highest i s.t. data[i] <= x
		unsigned int r = size; // lowest i s.t. data[i] > x
		unsigned int d = size;
		while (d > 1u)
		{
			const unsigned int m = l + d/2;
			if (x < data[m])
				r = m;
			else
				l = m;
			d = r - l;
		}
		return l;
	}

	template <class V, class T> int SegmentSearch::binary_search_right_inclusive(const V& data, unsigned int size, const T& x)
	{
		if (size == 0)
			throw std::domain_error("Zero-sized array");
		if (x <= data[0])
			return -1;

		// there is a segment which contains x

		// logically, we assume that data[size] == +Infinity
		unsigned int l = 0; // highest i s.t. data[i] < x
		unsigned int r = size; // lowest i s.t. data[i] >= x
		unsigned int d = size;
		while (d > 1)
		{
			const unsigned int m = l + d/2;
			if (x <= data[m])
				r = m;
			else
				l = m;
			d = r - l;
		}
		return l;
	}
}

#endif // __MATH_SEGMENT_SEARCH_H
