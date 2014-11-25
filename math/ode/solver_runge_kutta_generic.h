#ifndef __SOLVER_RUNGE_KUTTA_GENERIC_H
#define __SOLVER_RUNGE_KUTTA_GENERIC_H

#include <cassert>
#include <stdexcept>
#include <utility>
#include <vector>
#include <boost/array.hpp>
#include <Eigen/Core>

namespace rql { 
	namespace math {
		namespace ode {

			template <class V, bool Adaptive> struct SolverRungeKuttaGenericMidPointWorkspace
			{
			};

			template <class V, bool Adaptive> struct SolverRungeKuttaGenericNorm
			{
				static double norm(const V& v);
			};

			template <class V> struct SolverRungeKuttaGenericNorm<V,false>
			{
				static double norm(const V& v) { return 0; }
			};			

			template <> struct SolverRungeKuttaGenericNorm<double,true>
			{
				static double norm(double v) { return std::abs(v); }
			};

			template <> struct SolverRungeKuttaGenericNorm<std::complex<double>,true>
			{
				static double norm(const std::complex<double>& v) { return std::abs(v); }
			};

			template<typename S, int R, int C, int O, int MR, int MC> struct SolverRungeKuttaGenericNorm<Eigen::Matrix<S,R,C,O,MR,MC>, true>
			{
				static double norm(const Eigen::Matrix<S,R,C,O,MR,MC>& m)
				{
					return m.norm();
				}
			};

			template <class V> struct SolverRungeKuttaGenericNormalizer
			{
				inline void normalize(V& state) const
				{
					; // do nothing
				}
			};

			//! Solves an equation dy/dt = f(y,t) using 4th order RK method			
			//! @tparam F Provides operator()(double t, const V& y, V& f) which calculates f(y,t)
			template <class V, bool Adaptive> class SolverRungeKuttaGeneric
			{
			public:
				typedef V y_type;
				SolverRungeKuttaGeneric();
				//! Initialize V variables to right dimensions
				SolverRungeKuttaGeneric(const V& zero);
				template <class F> void solve(F& ode, V& state, double prev_t, double next_t)
				{
					solve(ode, m_null_normalizer, state, prev_t, next_t);
				}
				template <class F, class Normalizer> void solve(F& ode, Normalizer& normalizer, V& state, double prev_t, double next_t);
				template <class F> void solve(F& ode, const V& prev_it, V& next_it, double prev_t, double next_t)
				{
					solve(ode, m_null_normalizer, prev_it, next_it, prev_t, next_t);
				}
				template <class F, class Normalizer> void solve(F& ode, Normalizer& normalizer, const V& prev_it, V& next_it, double prev_t, double next_t);
				void set_norm_range(double new_norm_range);
				static const unsigned int MAX_NBR_SUBDIVISIONS = 8;
			private:
				//! Init common params
				void init();
				//! Save the change in m_tmp;
				template <class F, bool CheckNorms> bool calc_change(F& ode, const V& prev, double prev_t, double next_t);
				template <class F, class Normalizer> void solve_impl(F& ode, Normalizer& normalizer, const V& prev_it, V& next_it, double prev_t, double next_t, unsigned int nbr_allowed_subdivisions);
				bool norms_diverged(double min_norm, double max_norm)
				{
					return min_norm > 0 && max_norm/min_norm > m_norm_range;
				}
				template <class Normalizer> void add_change(const V& prev, V& next, Normalizer& normalizer) const
				{
					next = prev;
					next += m_tmp;
					normalizer.normalize(next);
				}
				double m_norm_range;
				SolverRungeKuttaGenericMidPointWorkspace<V,Adaptive> m_mid_pts;
				V m_k1;
				V m_k2;
				V m_k3;
				V m_k4;
				V m_tmp;
				const SolverRungeKuttaGenericNormalizer<V> m_null_normalizer;
			};

			template <class V> struct SolverRungeKuttaGenericMidPointWorkspace<V,true>
			{
				boost::array<V, SolverRungeKuttaGeneric<V,true>::MAX_NBR_SUBDIVISIONS> m_points;
			};

			template <class V> struct SolverRungeKuttaGenericMidPointWorkspace<V,false>
			{
				boost::array<V, 0u> m_points;
			};

			template <class V, bool Adaptive> SolverRungeKuttaGeneric<V,Adaptive>::SolverRungeKuttaGeneric()
				: m_null_normalizer()
			{
				init();
			}

			template <class V, bool Adaptive> SolverRungeKuttaGeneric<V,Adaptive>::SolverRungeKuttaGeneric(const V& zero)
				: m_k1(zero), m_k2(zero), m_k3(zero), m_k4(zero), m_tmp(zero), m_null_normalizer()
			{
				init();
			}

			template <class V, bool Adaptive> void SolverRungeKuttaGeneric<V,Adaptive>::init()
			{
				set_norm_range(1.2);
			}

			template <class V, bool Adaptive> void SolverRungeKuttaGeneric<V,Adaptive>::set_norm_range(double new_norm_range)
			{
				m_norm_range = new_norm_range;
			}

			template <class V, bool Adaptive> template <class F, class Normalizer> void SolverRungeKuttaGeneric<V,Adaptive>::solve(F& ode, Normalizer& normalizer, V& state, double prev_t, double next_t)
			{
				solve(ode, normalizer, state, state, prev_t, next_t);
			}

			template <class V, bool Adaptive> template <class F, class Normalizer> void SolverRungeKuttaGeneric<V,Adaptive>::solve(F& ode, Normalizer& normalizer, const V& prev, V& next, const double prev_t, const double next_t)
			{
				if (Adaptive) {
					solve_impl(ode, normalizer, prev, next, prev_t, next_t, MAX_NBR_SUBDIVISIONS);
				} else {
					calc_change<F,false>(ode, prev, prev_t, next_t);
					add_change(prev, next, normalizer);
				}
			}

			template <class V, bool Adaptive> template <class F, class Normalizer> void SolverRungeKuttaGeneric<V,Adaptive>::solve_impl(F& ode, Normalizer& normalizer, const V& prev, V& next, const double prev_t, const double next_t, unsigned int nbr_allowed_subdivisions)
			{
				if (Adaptive && nbr_allowed_subdivisions > 0) {
					const bool step_ok = calc_change<F,true>(ode, prev, prev_t, next_t);
					if (step_ok) {
						add_change(prev, next, normalizer);
					} else {
						const double mid_t = 0.5*prev_t + 0.5*next_t;
						if (mid_t <= prev_t || mid_t >= next_t) {
							throw std::runtime_error("Time step too small");
						} else {
							const unsigned int new_nbr_subdiv = nbr_allowed_subdivisions - 1;
							V& mid_point = m_mid_pts.m_points[new_nbr_subdiv];
							//std::cout << "going to level " << new_nbr_subdiv << std::endl;
							solve_impl(ode, normalizer, prev, mid_point, prev_t, mid_t, new_nbr_subdiv);
							solve_impl(ode, normalizer, mid_point, next, mid_t, next_t, new_nbr_subdiv);							
						}
					}
				} else {
					calc_change<F,false>(ode, prev, prev_t, next_t);
					add_change(prev, next, normalizer);
				}
			}

			template <class V, bool Adaptive> static void update_norms(const V& v, double& min_norm, double& max_norm)
			{
				const double norm = SolverRungeKuttaGenericNorm<V,Adaptive>::norm(v);
				min_norm = std::min(min_norm, norm);
				max_norm = std::max(max_norm, norm);
			}			

			template <class V, bool Adaptive> template <class F, bool CheckNorms> bool SolverRungeKuttaGeneric<V,Adaptive>::calc_change(F& ode, const V& prev, double prev_t, double next_t)
			{
				const double h = next_t - prev_t;
				assert(h > 0);
				double min_norm;
				double max_norm;

				ode(prev_t, prev, m_k1);
				m_k1 *= h;
				if (Adaptive && CheckNorms) {
					min_norm = SolverRungeKuttaGenericNorm<V,Adaptive>::norm(m_k1);
					max_norm = min_norm;
				}

				m_tmp = m_k1;
				m_tmp *= 0.5;
				m_tmp += prev;

				ode(prev_t + 0.5*h, m_tmp, m_k2);
				m_k2 *= h;
				if (Adaptive && CheckNorms) {
					update_norms<V,Adaptive>(m_k2, min_norm, max_norm);
					if (norms_diverged(min_norm, max_norm))
						return false;
				}

				m_tmp = m_k2;
				m_tmp *= 0.5;
				m_tmp += prev;	

				ode(prev_t + 0.5*h, m_tmp, m_k3);
				m_k3 *= h;
				if (Adaptive && CheckNorms) {
					update_norms<V,Adaptive>(m_k3, min_norm, max_norm);
					if (norms_diverged(min_norm, max_norm))
						return false;
				}

				m_tmp = prev;
				m_tmp += m_k3;

				ode(next_t, m_tmp, m_k4);
				m_k4 *= h;
				if (Adaptive && CheckNorms) {
					update_norms<V,Adaptive>(m_k4, min_norm, max_norm);
					if (norms_diverged(min_norm, max_norm))
						return false;
				}

				m_tmp = m_k2;
				m_tmp += m_k3;
				m_tmp *= 2;
				m_tmp += m_k1;
				m_tmp += m_k4;
				m_tmp /= 6.0;

				if (Adaptive && CheckNorms) {
					const double delta_norm = SolverRungeKuttaGenericNorm<V,Adaptive>::norm(m_tmp);
					const double prev_norm = SolverRungeKuttaGenericNorm<V,Adaptive>::norm(prev);
					if (prev_norm > 0 && delta_norm / prev_norm > (m_norm_range - 1.0))
						return false;
				}

				return true;
			}
		}
	}
}

#endif // __SOLVER_RUNGE_KUTTA_GENERIC_H
