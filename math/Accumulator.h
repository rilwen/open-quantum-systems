#pragma once
#include "MathCore.h"

namespace rql { namespace math {

	/**
	 * Accumulates a sum of numbers trying to avoid numerical errors. Uses Kahan summation.
	 * @tparam T Any floating-point number type, like double or float.
	 */
	template <typename T>
	class Accumulator
	{
		public:
		/**
		 * Default constructor. Initializes the sum to exactly 0.
		 */
		Accumulator();

		/**
		 * Copy constructor. Copies internal state.
		 * @param acc Original accumulator.
		 */
		Accumulator(const Accumulator<T>& acc);

		/**
		 * Destructor.
		 */
		~Accumulator();

		/** Compares the sum stored in the accumulator with another value.
		 * @param value Value to be compared with.
		 */
		bool operator==(const T& value) const;

		/**
		 * Increment the sum by adding new value. Note that adding zero will change the internal state
		 * and the results of subsequent increments.
		 * @return Reference to *this.
		 */
		Accumulator<T> & operator+=(const T& value);

		/**
		 * Assign a new sum. Clears the internal state.
		 * @param value Assigned value.
		 * @return Reference to *this.
		 */
		Accumulator<T> & operator=(const T& value);

		/**
		 * Assign a new accumulator. Copies both the sum and the internal state.
		 * @param acc Original accumulator.
		 * @return Reference to *this.
		 */
		Accumulator<T> & operator=(const Accumulator<T>& acc);

		/**
		 * Convert the accumulator to the value of the sum.
		 * @return Value equal to the result of calling sum().
		 */
		operator T() const;

		/**
		 * Return the accumulated sum.
		 * @return Accumulated sum.
		 */
		T sum() const;

		/**
		 * Pull out the sum from the accumulator. Useful for functors.
		 * @param acc Accumulator.
		 * @return Value equal to the result of calling sum() on the accumulator.
		 */
		static T get_sum(const Accumulator& acc);
	private:
		T c_;
		T y_;
		T sum_;
		T t_;
		inline void accumulate(const T& value);
	};

	template<typename T>
	Accumulator<T>::Accumulator()
	: c_(0.0), y_(0.0), sum_(0.0), t_(0.0)
	{
	}

	template<typename T>
	Accumulator<T>::Accumulator(const Accumulator<T>& acc)
	: c_(acc.c_), y_(acc.y_), sum_(acc.sum_), t_(acc.t_)
	{
	}

	template<typename T>
	Accumulator<T>::~Accumulator()
	{
	}

	#ifdef _MSC_VER
	#pragma optimize("",off)
	#endif
	template<typename T>
	Accumulator<T>& Accumulator<T>::operator+=(const T& value)
	{
		y_ = value - c_;
		t_ = sum_ + y_;
		c_ = (t_ - sum_) - y_;
		sum_ = t_;
		return *this;
	}
	#ifdef _MSC_VER
	#pragma optimize("",on)
	#endif

	template<typename T>
	Accumulator<T>& Accumulator<T>::operator=(const T& value)
	{
		y_ = 0;
		c_ = 0;
		t_ = 0;
		sum_ = value;
		return *this;
	}

	template<typename T>
	Accumulator<T>& Accumulator<T>::operator=(const Accumulator<T>& acc)
	{
		y_ = acc.y_;
		c_ = acc.c_;
		sum_ = acc.sum_;
		t_ = acc.t_;
		return *this;
	}

	template<typename T>
	inline T Accumulator<T>::sum() const
	{
		return sum_;
	}

	template<typename T>
	inline Accumulator<T>::operator T() const
	{
		return sum_;
	}

	template<typename T>
	inline bool Accumulator<T>::operator==(const T& value) const
	{
		return sum_ == value;
	}

	template<typename T>
	inline T Accumulator<T>::get_sum(const Accumulator<T>& acc)
	{
		return acc.sum();
	}

}}