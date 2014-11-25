#ifndef __CFA_SIMULATOR_H
#define __CFA_SIMULATOR_H

#include <Eigen/Core>
#include <math/ode/solver_runge_kutta_generic.h>

namespace cfa {
	template <class Method> class Simulator
	{
	public:
		typedef Method method_type;
		//! H0 is NxN
		//! g is KxN
		//! omega is Kx1
		Simulator(const Eigen::MatrixXcd& H0, const Eigen::MatrixXcd& g, const Eigen::VectorXd& omega, double kBT = 0, bool fixed_bath = false);
		Simulator(const Method& method);
		void reset_H0(const Eigen::MatrixXcd& H0) { m_method.reset_H0(H0); }
		template <class F> void simulate(const Eigen::VectorXcd& initial_state, double dt, size_t nbr_steps, F& processor);
	private:
		Method m_method;
	};

	template <class M> Simulator<M>::Simulator(const Eigen::MatrixXcd& H0, const Eigen::MatrixXcd& g, const Eigen::VectorXd& omega, double kBT, bool fixed_bath)
		: m_method(H0, g, omega, kBT, fixed_bath)
	{
	}

	template <class M> Simulator<M>::Simulator(const M& method)
		: m_method(method)
	{
	}

	template <class M> template <class F> void Simulator<M>::simulate(const Eigen::VectorXcd& initial_state, double dt, size_t nbr_steps, F& processor)
	{
		typedef typename M::data_type data_type;
		data_type y;
		m_method.init_from_pure_state(initial_state, y);
		rql::math::ode::SolverRungeKuttaGeneric<data_type,true> solver(y);
		solver.set_norm_range(1.2);
		processor(m_method, y, 0, 0.0);
		typename M::normalizer_type normalizer(m_method);
		for (size_t n = 1; n <= nbr_steps; ++n) {
			const double t = n*dt;
			solver.solve(m_method, normalizer, y, (n-1)*dt, t);
			m_method.cfa_normalize(y, t);
			processor(m_method, y, n, t);
		}
	}
}

#endif // __CFA_TRANSPORT_SIMULATOR_H
