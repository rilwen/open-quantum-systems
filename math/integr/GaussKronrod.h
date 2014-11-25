#pragma once
#include <cassert>
#include <cstdio>
#include <utility>
#include "../MathCore.h"

namespace rql
{
	namespace math
	{
		namespace integr
		{
			// Integrator which uses 2N+1 point Gauss-Kronrod rule
			template <size_t N>
			class GaussKronrod
			{
			public:
				// function maps double to V
				template <typename F, typename V>
				void integrate(F function, double x0, double x1, V& sum) const;

				// double -> double
				template <typename F>
				double integrate(F function, double x0, double x1) const;

				// function maps double to V
				// return pair (result, result - less_precise_result); the latter is a measure of error of eac
				// element of the vector integral
				template <typename F, typename V>
				void integrate_with_error(F function, double x0, double x1, V& integral, V& error) const;

				// double -> double
				template <typename F>
				std::pair<double,double> integrate_with_error(F function, double x0, double x1) const;

				inline size_t size() const
				{
					return 2 * N + 1;
				}

				inline double position(size_t i, double x0, double x1) const
				{
					return x0 + scale(x0, x1)*(xgk_[i] - start_);
				}

				inline double weight(size_t i, double x0, double x1) const
				{
					return scale(x0, x1) * wgk_[i];
				}
			private:

				inline double scale(double x0, double x1) const
				{
					return (x1 - x0) / (end_ - start_);
				}
				RQL_MATH_API_DEBUG static const double start_;
				RQL_MATH_API_DEBUG static const double end_;
				RQL_MATH_API_DEBUG static const double xgk_[2*N+1];
				RQL_MATH_API_DEBUG static const double wgk_[2*N+1];
				RQL_MATH_API_DEBUG static const double xg_[N];
				RQL_MATH_API_DEBUG static const double wg_[N];
			};

			template<size_t N>
			template<typename F, typename V>
			void GaussKronrod<N>::integrate(F f, double x0, double x1, V& sum) const
			{
				assert(x1 >= x0);
				if (x1 == x0) {
					return;
				}

				const double scale = this->scale(x0, x1);
				const size_t len = this->size();
				sum = wgk_[0] * f(x0 + scale * (xgk_[0] - start_));
				for (size_t i = 1; i < len; ++i) 
				{
					sum += wgk_[i] * f(x0 + scale * (xgk_[i] - start_));
				}
				sum *= scale;
			}

			template<size_t N>
			template<typename F>
			double GaussKronrod<N>::integrate(F f, double x0, double x1) const
			{
				double result = 0.0;
				integrate(f, x0, x1, result);
				return result;
			}

			template<size_t N>
			template<typename F, typename V>
			void GaussKronrod<N>::integrate_with_error(F f, double x0, double x1, V& integral, V& error) const
			{
				assert(x1 >= x0);
				if (x1 == x0) {
					return;
				}

				const double scale = this->scale(x0, x1);
				integral = f(x0 + scale * (xgk_[0] - start_))*wgk_[0];
				error = wg_[0] * f(x0 + scale * (xg_[0] - start_));
				for (size_t i = 1; i < N; ++i) {					
					integral += f(x0 + scale * (xgk_[i] - start_))*wgk_[i];
					error += wg_[i] * f(x0 + scale * (xg_[i] - start_));
				}
				for (size_t i = N; i <= 2*N; ++i) {
					integral += wgk_[i] * f(x0 + scale * (xgk_[i] - start_));
				}
				error -= integral;
				integral *= scale;
				error *= -scale;
			}

			template<size_t N>
			template<typename F>
			std::pair<double,double> GaussKronrod<N>::integrate_with_error(F f, double x0, double x1) const
			{
				double integral = 0.0;
				double error = 0.0;
				integrate_with_error(f, x0, x1, integral, error);
				return std::pair<double,double>(integral, error);
			}

			typedef GaussKronrod < 30 > GaussKronrod61;
			typedef GaussKronrod < 20 > GaussKronrod41;
			typedef GaussKronrod < 10 > GaussKronrod21;
		}
	}
}
