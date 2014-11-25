#ifndef __CFA_ODE_H
#define __CFA_ODE_H

#include <Eigen/Core>
#include "core.h"

namespace cfa {
	class ClassicalFieldApproximation
	{
	public:
		typedef Eigen::VectorXcd data_type;
	public:
		//! H0 is NxN
		//! g is KxN
		//! omega is Kx1
		ClassicalFieldApproximation(const Eigen::MatrixXcd& H0, const Eigen::MatrixXcd& g, const Eigen::VectorXd& omega, double kBT = 0);
		//! Calculate dy_dt = f(y, t)
		void operator()(double t, const data_type& y, data_type& dy_dt);
		size_t dim() const { return m_dim; }
		void init_from_pure_state(const Eigen::VectorXcd& exciton, data_type& y) const;
		void init_from_pure_state(const Eigen::VectorXcd& exciton, const Eigen::VectorXcd& bath, data_type& y) const;
		//! int_0^t <Psi(0)|H(t)|Psi(0)> dt
		double energy_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const;
		void exciton_density_matrix(const data_type& y, const Eigen::VectorXcd& initial_pure_state, Eigen::MatrixXcd& rho) const;
		std::complex<double> exciton_density_matrix_trace(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const;
		void reset_H0(const Eigen::MatrixXcd& H0);
		struct Normalizer
		{
			Normalizer(const ClassicalFieldApproximation&) {}
			inline static void normalize(data_type& state)
			{
				; // do nothing
			}
		};
		typedef Normalizer normalizer_type;
		static const bool WRAP_RHO_EIGENVALUES = false;
		static const double RHO_EIGENVALUE_TOLERANCE;
		//! Normalize the state (does nothing)
		void cfa_normalize(data_type&, double) const {}
		static const char* name() { return "CFA"; }
	private:
		void init_exciton_part(const Eigen::VectorXcd& exciton, data_type& y) const;
		size_t bath_index(size_t k) const { return m_nbr_exciton_values + k; }
	private:
		Eigen::MatrixXcd m_H0;
		Eigen::MatrixXcd m_g;
		Eigen::VectorXd m_omega;
		double m_kBT;
		Eigen::VectorXcd m_exp_w;
		Eigen::VectorXcd m_work_g;
		size_t m_dim;
		size_t m_nbr_bath_modes;
		size_t m_nbr_sites;
		size_t m_nbr_exciton_values;	
	};
}

#endif // __CFA_ODE_H
