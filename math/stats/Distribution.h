#pragma once

#include <vector>
#include <cassert>
#include <algorithm>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include "../MathCore.h"

namespace rql { namespace math { namespace stats {
	/** Probability distribution for a discrete random variable.
	 * @tparam T Scalar type.
	 */
	template <typename T>
	class Distribution
	{
	public:
		/** Simple constructor.
		 * @param[in] values Array of random variable values in ascending order.
		 * @param[in] probabilities Corresponding probabilities. Must sum to 1 with #PROB_EPS tolerance.
		 * @tparam V1 vector-like class with size(), begin() and end() containing T
		 * @tparam V2 vector-like class with size(), begin() and end() containing double
		 */
		template <class V1, class V2>
		Distribution(const V1& values, const V2& probabilities, typename boost::disable_if<boost::is_scalar<V2> >::type* dummy = 0);

		/** Construct a distribution of multiples of unit, starting from zero.
		 * @param[in] probabilities Array of probabilities for the multiples. Length of the array - 1 is the maximum multiplicity. Must sum to 1 with PROB_EPS tolerance.
		 * @param[in] unit Unit of the random variable values.
		 * @tparam V vector-like class with size(), begin() and end() containing double
		 */
		template <class V>
		Distribution(const V& probabilities, T unit);

		Distribution(const Distribution<T>& other);

		/** Holds a single value with 100% probability.
		 * @param value Single certain value. 0 by default.
		 */
		Distribution(T value = 0);

		virtual ~Distribution(void);

		/** Assingment operator
		 * @param[in] other Other distribution
		 */
		Distribution<T> & operator=(const Distribution<T> other);

		/**
		 * Number of probabilities in the distribution.
		 */
		size_t size() const;

		/** Return a given value.
		 * @param[in] idx Value index, less than size().
		 * @return Copy of the value.
		 */
		T value(size_t idx) const;

		/** Return the probability of a given value.
		 * @param[in] idx Value index, less than size().
		 * @return Probability.
		 */
		double probability(size_t idx) const;

		/**
		 * Maximum value in the distribution.
		 * @return Copy of the maximum.
		 */
		T max() const;

		/**
		 * Minimum value in the distribution.
		 * @return Copy of the minimum.
		 */
		T min() const;

		/**
		 * Expected value of the random variable. Makes sense only for types T which can be cast to double.
		 */
		double expected_value() const;

		/** Calculate expected value of the function.
		 * @tparam F Function type with operator()(T), returning the function value as double.
		 * @param function Function object. Must be calculable for every value in the distribution.
		 * @return Expected value of the function
		 */
		template <typename F>
		double expected_value(F function) const;

		/** Calculate expected value of the function which takes some parameters.
		 * @tparam F Function type with operator()(T, const T&), returning the function value as double.
		 * @tparam P Parameter object type.
		 * @param function Function object. Must be calculable for every value in the distribution.
		 * @param params Parameter object.
		 */
		template <typename F, typename P>
		double expected_value(F function, const P& params) const;

		/** Convolute this distribution with another one. See @link<convolute_discrete_positive(const std::vector<Distribution<T> >& distributions)
		 * for the properties required of this and the other distribution, and the description of the returned value.
		 * @tparam T See @link<convolute_discrete_positive(const std::vector<Distribution<T> >& distributions)
		 * @param[in] distribution Distribution to be convoluted.
		 * @return Convoluted distribution.
		 */
		Distribution<T> convolute_discrete_positive(const Distribution<T>& distribution) const;

		/** Convolute probability distributions, calculating the probability distribution of a random variable
		 * which is the sum of random variables which have the probability distributions stored in the array
		 * passed as argument.
		 * @tparam T Must be an integer type. Does not have to be signed.
		 * @param[in] distributions Array of convoluted distributions. Must have at least 1 element.
		 * @return Convoluted distribution, with probabilities present for every value 0 and max(), inclusive. Some probabilities can be zero.
		 */
		static Distribution<T> convolute_discrete_positive(const std::vector<Distribution<T> >& distributions);

		/** Similar to @link<convolute_discrete_positive(const std::vector<Distribution<T> >& distributions),
		 * but the distribution is saved as an array of probabilities for values from 0 to maximum value of
		 * the sum, inclusive.
		 * @tparam T See @link<convolute_discrete_positive(const std::vector<Distribution<T> >& distributions)
		 * @param[in] distributions Array of convoluted distributions. Must have at least 1 element.
		 * @param[in,out] new_probs Destination array-like object for the computed probabilities which is also used as workspace.
		 * It is assumed that it holds enough space for all probabilities for values from 0 to maximum sum, which were all
		 * initialized with 0. At output, new_probs[i] is the probability that the sum of random variables is equal to i.
		 */
		template <typename Vec>
		static void convolute_discrete_positive(const std::vector<Distribution<T> >& distributions, Vec& new_probs);

		template <typename Vec>
		static void convolute_discrete_positive(const Distribution<T>& distro, Vec& probs);

		/**
		 * @return new curr_max
		 */
		template <typename Vec>
		static T convolute_discrete_positive(const Distribution<T>& distro, Vec& probs, T curr_max);

		/**
		 * Copy values to a container for T with push_back().
		 * @tparam Container Can be std::vector or std::list
		 * @param[in] cont Container reference.
		 */
		template <typename Container>
		void copy_values(Container& cont) const;

		/**
		 * Copy probabilities to a container for double with push_back().
		 * @tparam Container Can be std::vector or std::list
		 * @param[in] cont Container reference.
		 */
		template <typename Container>
		void copy_probabilities(Container& cont) const;

		typedef typename std::vector<T>::const_iterator const_value_iterator;
		typedef typename std::vector<double>::const_iterator const_prob_iterator;

		const_value_iterator values_begin() const
		{
			return values_.begin();
		}

		const_value_iterator values_end() const
		{
			return values_.end();
		}

		const_prob_iterator probabilities_begin() const
		{
			return probabilities_.begin();
		}

		const_prob_iterator probabilities_end() const
		{
			return probabilities_.end();
		}

		//!Tolerance for the sum of probabilities.
		static double prob_epsilon();
	private:
		std::vector<T> values_;
		std::vector<double> probabilities_;
	};

	RQL_MATH_API_DEBUG bool probs_are_ok(const std::vector<double>& probabilities);

	template <typename T>
	Distribution<T>::Distribution(T value)
	: values_(1, value), probabilities_(1, 1.0)
	{
	}

	template <typename T>
	bool values_are_ok(const std::vector<T>& values)
	{
		bool is_ok = true;
		for (size_t i = 1; i < values.size(); ++i) {
			is_ok &= values[i] > values[i - 1];
		}
		return is_ok;
	}

	template <typename T>
	template <class V1, class V2>
	Distribution<T>::Distribution(const V1& values, const V2& probabilities, typename boost::disable_if<boost::is_scalar<V2> >::type* dummy)
	: values_(values.size()), probabilities_(probabilities.size())
	{
		assert(values.size() == probabilities.size());
		std::copy(values.begin(), values.end(), values_.begin());
		std::copy(probabilities.begin(), probabilities.end(), probabilities_.begin());
		assert(probs_are_ok(probabilities_));
		assert(values_are_ok(values_));
	}

	template <typename T>
	std::vector<T> uniform(size_t cnt, T unit)
	{
		std::vector<T> v(cnt);
		for (size_t i = 0; i < cnt; ++i) {
			v[i] = i*unit;
		}
		return v;
	}

	template <typename T>
	template <class V>
	Distribution<T>::Distribution(const V& probabilities, T unit)
	: values_(uniform(probabilities.size(), unit)), probabilities_(probabilities.size())
	{
		assert(unit > 0);
		std::copy(probabilities.begin(), probabilities.end(), probabilities_.begin());
		assert(probs_are_ok(probabilities_));
		assert(values_are_ok(values_));
	}

	template <typename T>
	Distribution<T>::Distribution(const Distribution<T>& other)
	: values_(other.values_), probabilities_(other.probabilities_)
	{
	}

	template <typename T>
	Distribution<T>::~Distribution(void)
	{
	}
	template <typename T>
	Distribution<T>& Distribution<T>::operator=(const Distribution<T> other)
	{
		if (this != &other) {
			values_ = other.values_;
			probabilities_ = other.probabilities_;
		}
		return *this;
	}
	template <typename T>
	inline size_t Distribution<T>::size() const
	{
		return values_.size();
	}
	template <typename T>
	inline T Distribution<T>::value(size_t idx) const
	{
	#ifndef NDEBUG
		return values_.at(idx);
	#else
		return values_[idx];
	#endif
	}
	template <typename T>
	inline double Distribution<T>::probability(size_t idx) const
	{
	#ifndef NDEBUG
		return probabilities_.at(idx);
	#else
		return probabilities_[idx];
	#endif
	}
	template <typename T>
	inline T Distribution<T>::max() const
	{
		return values_.at(values_.size() - 1);
	}
	template <typename T>
	inline T Distribution<T>::min() const
	{
		return values_.at(0);
	}
	template <typename T>
	template <typename F>
	double Distribution<T>::expected_value(F function) const
	{
		double sum = 0;
		for (size_t i = 0; i < values_.size(); ++i) {
			sum += function(values_.at(i)) * probabilities_.at(i);
		}
		return sum;
	}
	template <typename T>
	template <typename F, typename P>
	double Distribution<T>::expected_value(F function, const P& params) const
	{
		double sum = 0;
		for (size_t i = 0; i < values_.size(); ++i) {
			sum += function(values_.at(i), params) * probabilities_.at(i);
		}
		return sum;
	}

	template <typename T>
	template <typename Vec>
	T Distribution<T>::convolute_discrete_positive(const Distribution<T>& d, Vec& new_probs, T curr_max_value)
	{
		assert( probs_are_ok(d.probabilities_) );
		curr_max_value += d.max();

		switch (d.size()) {
			case 1:
			{
				// Shift the convoluted distribution by d->value(0)
				const T value = d.value(0);
				T reverse_k = 0;
				while (reverse_k <= curr_max_value - value) {
					const T k = curr_max_value - reverse_k;
					new_probs[k] = new_probs[k - value];
					++reverse_k;
				}
				while (reverse_k <= curr_max_value) {
					new_probs[curr_max_value - reverse_k] = 0;
					++reverse_k;
				}
			}
				break;
			case 2:
			{
				const T value0 = d.value(0);
				const double probability0 = d.probability(0);
				assert(probability0 >= -Distribution<double>::prob_epsilon());
				assert(probability0 <= 1+Distribution<double>::prob_epsilon());
				const T value1 = d.value(1);
				assert(value0 < value1);
				const double probability1 = d.probability(1);
				assert(probability1 >= -Distribution<double>::prob_epsilon());
				assert(probability1 <= 1+Distribution<double>::prob_epsilon());
				T reverse_k = 0;
				T bnd = curr_max_value - value1;
				T k = curr_max_value;
				while (reverse_k <= bnd) {
					new_probs[k] = new_probs[k - value1] * probability1 + new_probs[k - value0] * probability0;
					++reverse_k;
					--k;
				}
				bnd = curr_max_value - value0;
				while (reverse_k <= bnd) {
					new_probs[k] = new_probs[k - value0] * probability0;
					++reverse_k;
					--k;
				}
				if (reverse_k <= curr_max_value) {
					for (k = 0; k <= curr_max_value - reverse_k; ++k) {
						new_probs[k] = 0;
					}
				}
			}
				break;
			default:
			{
				T top = curr_max_value;
				for (typename std::vector<T>::const_reverse_iterator vr = d.values_.rbegin(); vr != d.values_.rend(); ++vr) {
					const T bottom = *vr;
					const size_t nbr_values = d.values_.rend() - vr;

					assert( top >= bottom );
					for (T reverse_k = 0; reverse_k <= top - bottom; ++reverse_k) {
						const T k = top - reverse_k;
						double sum = 0;
						for (size_t k2 = 0; k2 < nbr_values; ++k2) {
	#ifndef NDEBUG
							sum += d.probability(k2) * new_probs[k - d.values_.at(k2)];
	#else
							sum += d.probability(k2) * new_probs[k - d.values_[k2]];
	#endif
						}
						new_probs[k] = sum;
					}

					top = bottom - 1;
				}
				if (d.value(0) > 0) {
					for (T reverse_k = curr_max_value - d.value(0) + 1; reverse_k <= curr_max_value; ++reverse_k) {
						new_probs[curr_max_value - reverse_k] = 0;
					}
				}
			}
		}
		return curr_max_value;
	}

	template <typename T>
	template <typename Vec>
	void Distribution<T>::convolute_discrete_positive(const std::vector<Distribution<T> >& distributions, Vec& new_probs)
	{
		new_probs[0] = 1;
		T curr_max_value = 0;

		for (typename std::vector<Distribution<T> >::const_iterator d = distributions.begin(); d != distributions.end(); ++d) {
			curr_max_value = convolute_discrete_positive(*d, new_probs, curr_max_value);
		}
	}
	template <typename T>
	Distribution<T> Distribution<T>::convolute_discrete_positive(const std::vector<Distribution<T> >& distributions)
	{
		T max_value = 0;
		for (typename std::vector<Distribution<T> >::const_iterator d = distributions.begin(); d != distributions.end(); ++d) {
			max_value += d->max();
		}
		std::vector<double> new_probs(max_value + 1);

		convolute_discrete_positive(distributions, new_probs);

		std::vector<T> v;
		v.reserve(new_probs.size());
		std::vector<double> p;
		p.reserve(new_probs.size());
		for (T i = 0; i <= max_value; ++i) {
			if (new_probs[i] != 0) {
				v.push_back(i);
				p.push_back(new_probs[i]);
			}
		}
		return Distribution<T > (v, p);
	}

	template <typename T>
	Distribution<T> Distribution<T>::convolute_discrete_positive(const Distribution<T>& distribution) const
	{
		std::vector<Distribution<T> > distros;
		distros.push_back(*this);
		distros.push_back(distribution);
		return Distribution<T>::convolute_discrete_positive(distros);
	}

	template <typename T>
	template <typename Container>
	void Distribution<T>::copy_values(Container& cont) const
	{
		for (typename std::vector<T>::const_iterator i = values_.begin(); i != values_.end(); ++i) {
			cont.push_back(*i);
		}
	}

	template <typename T>
	template <typename Container>
	void Distribution<T>::copy_probabilities(Container& cont) const
	{
		for (typename std::vector<double>::const_iterator i = probabilities_.begin(); i != probabilities_.end(); ++i) {
			cont.push_back(*i);
		}
	}

	template <typename T>
	double Distribution<T>::expected_value() const
	{
		double sum = 0;
		for (size_t i = 0; i < size(); ++i) {
			sum += static_cast<double>(values_.at(i)) * probabilities_.at(i);
		}
		return sum;
	}

	template <>
	template <typename Vec>
	int Distribution<int>::convolute_discrete_positive(const Distribution<int>& d, Vec& new_probs, int curr_max_value)
	{
		curr_max_value += d.max();
		const int value0 = d.value(0);
		switch (d.size()) {
			case 1:
				// Shift the convoluted distribution by d.value(0)
				if (value0 != 0) {
					int k = curr_max_value;
					for (; k >= value0; --k) {
						new_probs[k] = new_probs[k - value0];
					}
					for (; k >= 0; --k) {
						new_probs[k] = 0;
					}
				}
				break;
			case 2:
			{
				const double probability0 = d.probability(0);
				assert(probability0 >= -Distribution<double>::prob_epsilon());
				assert(probability0 <= 1+Distribution<double>::prob_epsilon());
				const int value1 = d.value(1);
				assert(value0 < value1);
				const double probability1 = d.probability(1);
				assert(probability1 >= -Distribution<double>::prob_epsilon());
				assert(probability1 <= 1+Distribution<double>::prob_epsilon());

				int k = curr_max_value;
				for (; k >= value1; --k) {
					new_probs[k] = new_probs[k - value1] * probability1 + new_probs[k - value0] * probability0;
				}
				for (; k >= value0; --k) {
					new_probs[k] = new_probs[k - value0] * probability0;
				}
				if (value0 != 0) {
					std::fill(&new_probs[0], &new_probs[value0], 0.0);
				}
			}
				break;
			default:
			{
				int top = curr_max_value;
				for (typename std::vector<int>::const_reverse_iterator vr = d.values_.rbegin(); vr != d.values_.rend(); ++vr) {
					const int bottom = *vr;
					const size_t nbr_values = d.values_.rend() - vr;

					assert( top >= bottom );
					for (int k = top; k >= bottom; --k) {
						double sum = 0;
						for (size_t k2 = 0; k2 < nbr_values; ++k2) {
	#ifndef NDEBUG
							sum += d.probability(k2) * new_probs[k - d.values_.at(k2)];
	#else
							sum += d.probability(k2) * new_probs[k - d.values_[k2]];
	#endif
						}
						new_probs[k] = sum;
					}

					top = bottom - 1;
				}
				if (value0 > 0) {
					for (int reverse_k = curr_max_value - value0 + 1; reverse_k <= curr_max_value; ++reverse_k) {
						new_probs[curr_max_value - reverse_k] = 0;
					}
				}
			}
		}
		return curr_max_value;
	}

	template <class T> inline double Distribution<T>::prob_epsilon()
	{
		return 1E-8;
	}
}}}