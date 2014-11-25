#ifndef __MFA_ODE_H
#define __MFA_ODE_H

#include <Eigen/Core>
#include <vector>
#include "core.h"

namespace cfa {
	class MatrixFieldApproximation
	{
	public:
		typedef Eigen::VectorXcd data_type;
	public:
		//! H0 is NxN
		//! g is KxN
		//! omega is Kx1
		CFA_API MatrixFieldApproximation(const Eigen::MatrixXcd& H0, const Eigen::MatrixXcd& g, const Eigen::VectorXd& omega, double kBT = 0, bool fixed_bath = false);
		//! Calculate dy_dt = f(y, t)
		CFA_API void operator()(double t, const data_type& y, data_type& dy_dt);
		CFA_API size_t dim() const { return m_dim; }
		CFA_API size_t operator_dim() const { return m_operator_dim; }
		CFA_API size_t nbr_sites() const { return m_nbr_sites; }
		CFA_API void init_from_pure_state(const Eigen::VectorXcd& exciton, data_type& y) const;
		CFA_API void init_from_pure_state(const Eigen::VectorXcd& exciton, const Eigen::VectorXcd& bath, data_type& y) const;
		CFA_API void exciton_density_matrix(const data_type& y, const Eigen::VectorXcd& initial_pure_state, Eigen::MatrixXcd& rho);
		CFA_API std::complex<double> exciton_density_matrix_trace(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const;
		CFA_API void reset_H0(const Eigen::MatrixXcd& H0);
		//! int_0^t <Psi(0)|H(t)|Psi(0)> dt
		CFA_API double energy_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const;
		CFA_API std::complex<double> evolution_operator_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const;
		struct Normalizer
		{
			CFA_API Normalizer(const MatrixFieldApproximation&) {}
			CFA_API inline void normalize(data_type& state) const
			{
				; // do nothing
			}
		};
		typedef Normalizer normalizer_type;
		CFA_API static const bool WRAP_RHO_EIGENVALUES = false;
		CFA_API static const double RHO_EIGENVALUE_TOLERANCE;
		//! Normalize the state (does nothing)
		CFA_API void cfa_normalize(data_type&, double) const {}
		static const char* name() { return "MFA"; }
	private:
		typedef Eigen::Map<Eigen::MatrixXcd> operator_type;
		typedef Eigen::Map<const Eigen::MatrixXcd> const_operator_type;
		operator_type get_operator(data_type& y, size_t offset) const
		{
			return operator_type(y.data() + offset, m_operator_dim, m_operator_dim);
		}
		const_operator_type get_operator(const data_type& y, size_t offset) const
		{
			return const_operator_type(y.data() + offset, m_operator_dim, m_operator_dim);
		}
		operator_type exciton_operator(data_type& y, size_t m, size_t n) const;
		const_operator_type exciton_operator(const data_type& y, size_t m, size_t n) const;
		operator_type bath_operator(data_type& y, size_t k) const;
		const_operator_type bath_operator(const data_type& y, size_t k) const;
		//! int_0^t H(t) dt
		operator_type hamiltonian_operator(data_type& y) const;
		//! int_0^t H(t) dt
		const_operator_type hamiltonian_operator(const data_type& y) const;
		size_t bath_index(size_t k) const { return m_nbr_exciton_values + k*m_nbr_operator_elements; }
		void init_exciton_part(data_type& y) const;
		template <class M> std::complex<double> mean_value(const M& op, const Eigen::VectorXcd& initial_pure_state) const;
	private:
		Eigen::MatrixXcd m_H0;
		Eigen::MatrixXcd m_g;
		Eigen::VectorXd m_omega;
		double m_kBT;
		bool m_fixed_bath;
		Eigen::VectorXcd m_exp_w;
		std::vector<Eigen::MatrixXcd> m_work_g;
		size_t m_dim;
		size_t m_nbr_bath_modes;
		size_t m_nbr_sites;		
		size_t m_operator_dim;
		size_t m_nbr_operator_elements;
		size_t m_nbr_exciton_values;
		size_t m_hamiltonian_offset;
	};
}

#endif // __MFA_ODE_H
