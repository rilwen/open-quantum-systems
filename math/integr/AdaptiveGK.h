#pragma once

#include <stdexcept>
#include <utility>
#include <cassert>
#include "../MathUtils.h"
#include "../MathCore.h"

namespace rql { namespace math { namespace integr {

	struct Normalizator1D
	{
		double operator()(double x) const
		{
			return std::abs(x);
		}
	};

	// size_t -- the the size_t of the GK integrator used (N,2*N+1)
	template <size_t N>
	class AdaptiveGK
	{
	public:
		// max_iter - maximum number of subdivisions
		// tolerance - tolerance for absolute error
		// blow_up -- when asked to return the error with the result, either throw runtime_error on
		// failure to converge (true), or return the result anyway (dflt false)
		AdaptiveGK(unsigned int max_subdiv, double tolerance, bool blow_up = false);

		// throws std::runtime_error if could not converge
		// "function" param maps double to V
		// Norm -- provides the operator()(const V& v) function which measures the norm of a result type
		template <typename F, typename V, typename Norm>
		void integrate(F function, double x0, double x1, V& integral, Norm normalizator) const;
		template <typename F>
		inline double integrate(F function, double x0, double x1) const;	

		// throws std::runtime_error if could not converge
		// "function" param maps double to V
		// Norm -- provides the norm(const V& v) function which measures the norm of a result type
		// returns a pair (result, error measure)
		template <typename F, typename V, typename Norm>
		double integrate_with_error(F function, double x0, double x1, V& integral, Norm normalizator = std::abs<V>) const;
	private:
		template <typename F, typename V, typename Norm>
		double integrate_impl(F function, double x0, double x1, unsigned int level, double current_tolerance, V& result, Norm normalizator) const;
		unsigned int max_subdiv_;
		double tolerance_;
		bool blow_up_;
		RQL_MATH_API static const double start_;
		RQL_MATH_API static const double end_;
		RQL_MATH_API static const double xgk_[2*N+1];
		RQL_MATH_API static const double wgk_[2*N+1];
		RQL_MATH_API static const double xg_[N];
		RQL_MATH_API static const double wg_[N];
	};

	template <size_t N>
	AdaptiveGK<N>::AdaptiveGK(unsigned int max_subdiv, double tolerance, bool blow_up)
	{
		max_subdiv_ = max_subdiv;
		tolerance_ = tolerance;
		blow_up_ = blow_up;
	}

	template <size_t N>
	template <typename F, typename V, typename Norm>
	void AdaptiveGK<N>::integrate(F function, double x0, double x1, V& result, Norm normalizator) const
	{
		const double error = integrate_with_error(function, x0, x1, result, normalizator);
		if (blow_up_ && error > tolerance_) {
			throw std::runtime_error("Could not converge: limit of subdivisions exceeded");
		}
	}

	template <size_t N>
	template <typename F>
	inline double AdaptiveGK<N>::integrate(F function, double x0, double x1) const
	{
		double result = 0.0;
		integrate<F,double,double(*)(double)>(function, x0, x1, result, std::abs);
		return result;
	}

	template <size_t N>
	template <typename F, typename V, typename Norm>
	double AdaptiveGK<N>::integrate_with_error(F function, double x0, double x1, V& integral, Norm normalizator) const
	{
		const double error = integrate_impl<F,V,Norm>(function, x0, x1, 0, tolerance_, integral, normalizator);
		return error;
	}

	template <size_t N>
	template <typename F, typename V, typename Norm>
	double AdaptiveGK<N>::integrate_impl(F function, double x0, double x1, unsigned int level, double current_tolerance, V& result, Norm normalizator) const
	{
		const double a = (x1 - x0) / (end_ - start_);
		assert( a > 0 );
		const double b = x0 - a * start_;
		result = function(b + a*xgk_[0]) * wgk_[0];
		V less_precise(function(b + a*xg_[0]) * wg_[0]);
		for (size_t i = 1; i < N; ++i) {
			less_precise += function(b + a*xg_[i]) * wg_[i];
			result += function(b + a*xgk_[i]) * wgk_[i];
		}
		for (size_t i = N; i <= 2*N; ++i) {
			result += function(b + a*xgk_[i]) * wgk_[i];
		}
		less_precise *= a;
		result *= a;
		const double error = normalizator(result - less_precise);
		if (error <= current_tolerance)
		{
			return error;
		} 
		else 
		{
			if (level < max_subdiv_)
			{
				const double mid_x = 0.5*x0 + 0.5*x1;
				if (mid_x > x0 && mid_x < x1) 
				{
					V left_integral;
					const double left_error = integrate_impl<F,V,Norm>(function, x0, mid_x, level + 1, current_tolerance / 2, left_integral, normalizator);
					V right_integral;
					const double right_error = integrate_impl<F,V,Norm>(function, mid_x, x1, level + 1, current_tolerance / 2, right_integral, normalizator);
					result = left_integral + right_integral;
					return left_error + right_error;
				} 
				else 
				{
					return error;
				}
			}
			else 
			{
				return error;
			}
		}
	}

	typedef AdaptiveGK<7> AdaptiveGK15;
	/*typedef AdaptiveGK<10> AdaptiveGK21;*/
	typedef AdaptiveGK<20> AdaptiveGK41;

}}}