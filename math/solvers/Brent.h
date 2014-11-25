#pragma once
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <cmath>
#include "../MathCore.h"
#include <boost/lexical_cast.hpp>

namespace rql
{
	namespace math
	{
		namespace solvers 
		{ 
			// Implementation detail
			inline bool opposite_signs(double f1, double f2)
			{
				return f1 * (f2 / std::abs(f2)) < 0;
			}

			// Implementation detail
			inline double effective_tolerance(double x, double tol)
			{
				// 4*machine_epsilon copied from Netlib's zeroin.f
				return tol + 4*std::numeric_limits<double>::epsilon()*std::abs(x);
			}

			RQL_MATH_API_DEBUG double root_finder_result(double root_position, double value_at_root, double x1, double f1, double x2, double f2, double& lower_bound, double& higher_bound);

			template <typename F>
			inline double safe_evaluate(const F& fun, double x)
			{
				const double f = fun(x);
				if (!finite(f)) {
					throw std::runtime_error("Function returned invalid number");
				}
				return f;
			}

	#ifndef NDEBUG
			struct RootFinderStats
			{
				int bisection_cnt;
				int inverse_quadratic_cnt;
				int linear_cnt;
				int rejection_cnt;

				void clear()
				{
					bisection_cnt = 0;
					rejection_cnt = 0;
					inverse_quadratic_cnt = 0;
					linear_cnt = 0;
				}
			};

			static RootFinderStats root_finder_stats;
	#endif

			// More-or-less Brent
			// Throws std::runtime_error if not converged or supplied bounds do not bracket a root.
			template <typename F>
			double find_root(F f, double low, double high, double xtol, double ytol, double& final_low, double& final_high)
			{
				assert( low < high );
				assert( ytol >= 0 );

	#ifndef NDEBUG
				root_finder_stats.clear();
	#endif

				double other_end = low;
				double best_guess = high;

				double other_value = safe_evaluate(f, other_end);
				double smallest_value = safe_evaluate(f, best_guess);
				if (!opposite_signs(other_value, smallest_value)) {
					throw std::runtime_error(
						(std::string("Bounds do not straddle root: ") +  boost::lexical_cast<std::string>(other_end) + " " + boost::lexical_cast<std::string>(other_value)
						+ " " + boost::lexical_cast<std::string>(best_guess) + " " + boost::lexical_cast<std::string>(smallest_value)).c_str()
						);
				}
				if (std::abs(other_value) <= ytol) {
					return root_finder_result(other_end, other_value, other_end, other_value, best_guess, smallest_value, final_low, final_high);
				}
				if (std::abs(smallest_value) <= ytol) {
					return root_finder_result(best_guess, smallest_value, other_end, other_value, best_guess, smallest_value, final_low, final_high);
				}

				if (std::abs(other_value) < std::abs(smallest_value)) {
					std::swap(other_end, best_guess);
					std::swap(other_value, smallest_value);
				}
				double prev_best_guess = other_end;
				double prev_smallest_value = other_value;
				double dist = other_end - best_guess;
				double abs_dist = std::abs(dist);
				double prev_delta = best_guess - other_end;

				double eff_tol = effective_tolerance(best_guess, xtol);
				while (abs_dist > eff_tol) {
					double delta=0.0;
					if (std::abs(smallest_value) < std::abs(prev_smallest_value) && 2*std::abs(prev_delta) >= eff_tol ) {
						// Try interpolation if the previous step brought improvement in the function value (closer to 0)
	#ifndef NDEBUG
						bool used_inv;
	#endif
						const double s = smallest_value / prev_smallest_value;
						const double m = 0.5*(other_end - best_guess);
						assert( std::abs(s) <= 1 );
						prev_delta = delta;
						if ((prev_smallest_value != other_value) && (prev_smallest_value != smallest_value)) {
							// Inverse quadratic interpolation
							const double q = prev_smallest_value / other_value;
							const double r = smallest_value / other_value;
							const double p = s * (2 * m * q * (q - r) - (best_guess - prev_best_guess) * (r - 1));
							const double t = (q - 1) * (r - 1) * (s - 1);
							delta = - p / t;
	#ifndef NDEBUG
							++root_finder_stats.inverse_quadratic_cnt;
							used_inv = true;
	#endif
						} else {
							// Linear interpolation
							delta = (2 * m * s) / (s - 1);
	#ifndef NDEBUG
							++root_finder_stats.linear_cnt;
							used_inv = false;
	#endif
						}

						// Avoid "oscillations" back and forth: interpolation is effective when "diving in" for the root
						// so we should be concentrating at one side of the bracket
						const double abs_delta = std::abs(delta);
						if ( abs_delta > 0.75 * abs_dist ) {
							// fall back to bisection
							delta = 0.5*dist;
							// this is a cheat to avoid doing too many bisections
							prev_delta = delta;
	#ifndef NDEBUG
							++root_finder_stats.bisection_cnt;
							++root_finder_stats.rejection_cnt;
							if (used_inv) {
								--root_finder_stats.inverse_quadratic_cnt;
							} else {
								--root_finder_stats.linear_cnt;
							}
	#endif
						}

					} else {
	#ifndef NDEBUG
						++root_finder_stats.bisection_cnt;
	#endif
						delta = 0.5*dist;
						// this is a cheat to avoid doing too many bisections
						prev_delta = delta;
					}

					// Step too small: force movement towards the center of the bracket.
					if (2*std::abs(delta) < eff_tol) {
						delta = other_end > best_guess ? eff_tol : -eff_tol;
						delta *= 0.5;
						// don't update the prev_delta -- another cheat!
					}

					prev_best_guess = best_guess;
					prev_smallest_value = smallest_value;
					best_guess += delta;
					smallest_value = safe_evaluate(f, best_guess);
					if (std::abs(smallest_value) <= ytol) {
						return root_finder_result(best_guess, smallest_value, other_end, other_value, prev_best_guess, prev_smallest_value, final_low, final_high);
					}
					assert(opposite_signs(prev_smallest_value, other_value));

					// Choose new bracket
					if (!opposite_signs(other_value, smallest_value)) {
						other_end = prev_best_guess;
						other_value = prev_smallest_value;
					}
					if (std::abs(other_value) < std::abs(smallest_value)) {
						std::swap(other_end, best_guess);
						std::swap(other_value, smallest_value);
					}

					eff_tol = effective_tolerance(best_guess, xtol);
					dist = other_end - best_guess;
					abs_dist = std::abs(dist);
				}
				return root_finder_result(best_guess, smallest_value, other_end, other_value, best_guess, smallest_value, final_low, final_high);
			}

	#ifndef NDEBUG
			// NR algorithm
			template <typename F>
			double brent(F f, double low, double high, double xtol, double ytol)
			{
				if (xtol < 0 || ytol < 0) {
					throw std::runtime_error("Tolerances cannot be negative");
				}

				assert(low < high);

				double a = low;
				double b = high;
				double c;

				double d = b - a; // current distance between root brackets
				double e = b - a;

				double fa;

				double fb;
				fa = safe_evaluate(f, a);
				if (std::abs(fa) <= ytol) {
					return a;
				}
				fb = safe_evaluate(f, b);
				if (std::abs(fb) <= ytol) {
					return b;
				}

				if ((fa < 0 && fb < 0) || (fa > 0 && fb > 0)) {
					throw std::runtime_error("Incorrect initial bounds.");
				}

				c = b;
				double fc = fb;

				for (int i = 0; i < 50; i++) {
					// Invariant: a is the previous best guess.

					// Keep the root between b and c.
					if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
						c = a;
						fc = fa;
						e = b - a;
						d = b - a;
					}

					// Keep the current best guess in b.
					if (std::abs(fc) < std::abs(fb)) {
						a = b;
						b = c;
						c = a;

						fa = fb;
						fb = fc;
						fc = fa;
					}

					double tolerance1 = 2 * std::numeric_limits<double>::epsilon() * std::abs(b) + 0.5 * xtol;

					double xm = 0.5 * (c - b);

					if ((std::abs(xm) < tolerance1) || std::abs(fb) <= ytol) {
						return b;
					}

					if ((std::abs(e) >= tolerance1) && (std::abs(fa) > std::abs(fb))) {
						double p;
						double q;

						const double s = fb / fa;
						if (a == c) {
							p = 2 * xm * s;
							q = 1 - s;
						} else {
							q = fa / fc;
							const double r = fb / fc;
							p = s * ((2 * xm * q * (q - r)) - ((b - a) * (r - 1)));
							q = (q - 1) * (r - 1) * (s - 1);
						}

						if (p > 0) {
							q = -q;
						}

						p = std::abs(p);

						const double min1 = 3 * xm * q - std::abs(tolerance1 * q);
						const double min2 = std::abs(e * q);

						if (2 * p < std::min(min1, min2)) {
							e = d;
							d = p / q;
						} else {
							d = xm;
							e = xm;
						}
					} else {
						d = xm;
						e = xm;
					}

					a = b;
					fa = fb;

					if (std::abs(d) > tolerance1) {
						b += d;
					} else {
						b += (xm >= 0 ? tolerance1 : -tolerance1);
					}

					fb = safe_evaluate(f, b);
					if (std::abs(fb) <= ytol) {
						return b;
					}
				}
				throw std::runtime_error("No solution found after MAXIMUM_ITERATIONS.");
			}
	#endif // NDEBUG

			// Throws std::runtime_error if not converged or supplied bounds do not bracket a root.
			template <typename F>
			double find_root(const F& f, double xa, double xb, double xtol, double ytol)
			{
				double final_xb, final_xa;
				return find_root(f, xa, xb, xtol, ytol, final_xa, final_xb);
			}
		}
	}
}
