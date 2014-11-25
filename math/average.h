#ifndef __MATH_AVERAGE_H
#define __MATH_AVERAGE_H

namespace rql { namespace math {

	template <class T> class Average
	{
	public:
		Average();
		//! Use when averaging non-scalar types, like matrices
		Average(const T& initial_zero);
		void update(const T& x);
		static void update(const T& x, T& average, double& counter);
		T value() const { return m_avg; }
		double counter() const { return m_cntr; }
	private:
		T m_avg;
		double m_cntr;
	};

	template <class T> Average<T>::Average()
		: m_avg(0), m_cntr(0.0)
	{
	}

	template <class T> Average<T>::Average(const T& initialZero)
		: m_avg(initialZero), m_cntr(0.0)
	{
	}

	template <class T> inline void Average<T>::update(const T& x, T& average, double& counter)
	{
		++counter;
		average += (x - average) / counter;
	}

	template <class T> void Average<T>::update(const T& x)
	{
		update(x, m_avg, m_cntr);
	}

}}

#endif // __MATH_AVERAGE_H

