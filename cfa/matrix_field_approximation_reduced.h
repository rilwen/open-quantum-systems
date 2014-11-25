#ifndef __MATRIX_FIELD_APPROXIMATION_REDUCED_H
#define __MATRIX_FIELD_APPROXIMATION_REDUCED_H

#include <Eigen/Core>
#include <Eigen/QR>
#include <vector>
#include <boost/array.hpp>
#include <cassert>
#include "core.h"

namespace cfa {
	class MatrixFieldApproximationReduced
	{
	public:
		typedef Eigen::VectorXcd data_type;
	public:
		//! H0 is NxN
		//! g is KxN
		//! omega is Kx1
		CFA_API MatrixFieldApproximationReduced(const Eigen::MatrixXcd& H0, const Eigen::MatrixXcd& g, const Eigen::VectorXd& omega, double kBT = 0, bool normalize = true, size_t pivot = 0);
		//! Calculate dy_dt = f(y, t)
		CFA_API void operator()(double t, const data_type& y, data_type& dy_dt);
		CFA_API size_t dim() const { return m_dim; }
		CFA_API void init_from_pure_state(const Eigen::VectorXcd& exciton, data_type& y) const;
		//void init_from_pure_state(const Eigen::VectorXcd& exciton, const Eigen::VectorXcd& bath, data_type& y) const;  // not sure if done right
		CFA_API void exciton_density_matrix(const data_type& y, const Eigen::VectorXcd& initial_pure_state, Eigen::MatrixXcd& rho);
		CFA_API std::complex<double> exciton_density_matrix_trace(const data_type& y, const Eigen::VectorXcd& initial_pure_state);
		CFA_API void reset_H0(const Eigen::MatrixXcd& H0) { m_H0 = H0; }
		CFA_API std::complex<double> evolution_operator_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const;
		//! int_0^t <Psi(0)|H(t)|Psi(0)> dt
		CFA_API double energy_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const;
		class IdentityNormalizer
		{
		public:
			CFA_API IdentityNormalizer(MatrixFieldApproximationReduced& ode);
			CFA_API void normalize(data_type& state);
		private:
			MatrixFieldApproximationReduced& m_ode;
		};
		typedef IdentityNormalizer normalizer_type;
		CFA_API static const bool WRAP_RHO_EIGENVALUES = false;
		CFA_API static const double RHO_EIGENVALUE_TOLERANCE;
		//! Normalize the state (does nothing)
		CFA_API void cfa_normalize(data_type&, double) const {}
		static const char* name() { return "MFAReduced"; }
	private:
		typedef Eigen::Map<Eigen::MatrixXcd> operator_type;
		typedef Eigen::Map<const Eigen::MatrixXcd> const_operator_type;		
		operator_type get_operator(data_type& y, size_t offset) const
		{
			assert(m_dim == y.size());
			assert(offset + m_nbr_sites*m_nbr_sites <= m_dim);
			return operator_type(y.data() + offset, m_nbr_sites, m_nbr_sites);
		}
		const_operator_type get_operator(const data_type& y, size_t offset) const
		{
			assert(m_dim == y.size());
			assert(offset + m_nbr_sites*m_nbr_sites <= m_dim);
			return const_operator_type(y.data() + offset, m_nbr_sites, m_nbr_sites);
		}
		//! return |m><m_pivot|
		operator_type reduced_exciton_operator(data_type& y, size_t m) const
		{
			assert(m < m_nbr_sites);
			return get_operator(y, m*m_nbr_operator_elements);
		}
		//! return |m><m_pivot|
		const_operator_type reduced_exciton_operator(const data_type& y, size_t m) const
		{
			assert(m < m_nbr_sites);
			return get_operator(y, m*m_nbr_operator_elements);
		}
		operator_type bath_operator(data_type& y, size_t k) const
		{
			assert(k < m_nbr_bath_modes);
			return get_operator(y, bath_index(k));
		}
		const_operator_type bath_operator(const data_type& y, size_t k) const
		{
			assert(k < m_nbr_bath_modes);
			return get_operator(y, bath_index(k));
		}
		//! int_0^t H(t) dt
		operator_type hamiltonian_operator(data_type& y) const
		{
			return get_operator(y, m_hamiltonian_offset);
		}
		//! int_0^t H(t) dt
		const_operator_type hamiltonian_operator(const data_type& y) const
		{
			return get_operator(y, m_hamiltonian_offset);
		}
		size_t bath_index(size_t k) const 
		{ 
			assert(k < m_nbr_bath_modes);
			return (m_nbr_sites + k) * m_nbr_operator_elements; 
		}
		void init_exciton_part(data_type& y) const;
		const Eigen::MatrixXcd& exciton_operator(const data_type& y, size_t m, size_t n, size_t tmpidx);		
	private:
		Eigen::MatrixXcd m_H0;
		Eigen::MatrixXcd m_g;
		Eigen::VectorXd m_omega;
		double m_kBT;
		size_t m_pivot;
		Eigen::VectorXcd m_exp_w;
		std::vector<Eigen::MatrixXcd> m_work_g;
		static const size_t NBR_TMP_MATRICES = 2;
		boost::array<Eigen::MatrixXcd, NBR_TMP_MATRICES> m_tmp_mat;
		Eigen::MatrixXd m_tmp_mat_real;
		static const size_t NBR_TMP_VECTORS = 2;
		boost::array<Eigen::VectorXd, NBR_TMP_VECTORS> m_tmp_vec;		
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> m_qr;
		size_t m_dim;
		size_t m_nbr_bath_modes;
		size_t m_nbr_sites;
		size_t m_operator_dim;
		size_t m_nbr_operator_elements;
		size_t m_hamiltonian_offset;
		double m_sqrt_nbr_sites;
		bool m_normalize;
	};
}

#endif // __MATRIX_FIELD_APPROXIMATION_REDUCED_H
