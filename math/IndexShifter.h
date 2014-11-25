#ifndef __RQL_MATH_INDEX_SHIFTER_H
#define __RQL_MATH_INDEX_SHIFTER_H

#include <stdexcept>
#include <cassert>

namespace rql {
	/** Allows indexing from different base than 0.
	* @tparam ValueType Type of the returned value
	*/
	template <class ValueType>
	class IndexShifter
	{
	public:
		IndexShifter(int base);

		//! @tparam Vector Vector type with operator[]
		template <class Vector> ValueType& idx(Vector& vector, int idx) const;

		//! @tparam Vector Vector type with operator[]
		template <class Vector> const ValueType& idx(const Vector& vector, int idx) const;
	private:
		IndexShifter& operator=(const IndexShifter&); // not implemented
		const int m_base;
	};

	template <class VT> IndexShifter<VT>::IndexShifter(int base)
		: m_base(base)
	{
	}

	template <class VT>  template <class V>
	inline VT& IndexShifter<VT>::idx(V& vector, int idx) const
	{
		assert(idx - m_base >= 0);
		return vector[idx - m_base];
	}

	template <class VT> template <class V> 
	inline const VT& IndexShifter<VT>::idx(const V& vector, int idx) const
	{
		assert(idx - m_base >= 0);
		return vector[idx - m_base];
	}



	/** Allows indexing from different base than 0. Stores a copy of the original vector.
	* @tparam Vector Vector with 0-based indexing via operator[]
	* @tparam ValueType Type of the returned value
	*/
	template <class Vector, class ValueType>
	class ShiftedVector
	{
	public:
		ShiftedVector(const Vector& vec, int base);
		ValueType& operator[](int idx);
		const ValueType& operator[](int idx) const;
	private:
		ShiftedVector& operator=(const ShiftedVector&); // not implemented
		Vector m_vec;
		IndexShifter<ValueType> m_shifter;
	};

	template <class Vector, class ValueType> ShiftedVector<Vector,ValueType>::ShiftedVector(const Vector& vec, int base)
		: m_vec(vec), m_shifter(base)
	{
	}

	template <class Vector, class ValueType> ValueType& ShiftedVector<Vector,ValueType>::operator[](int idx)
	{
		return m_shifter.idx(m_vec,idx);
	}

	template <class Vector, class ValueType> const ValueType& ShiftedVector<Vector,ValueType>::operator[](int idx) const
	{
		return m_shifter.idx(m_vec,idx);
	}
}

#endif // __RQL_MATH_INDEX_SHIFTER_H
