#ifndef __CFA_ROA_LHO_FULL_H
#define __CFA_ROA_LHO_FULL_H

#include <Eigen/Core>
#include <Eigen/SVD>
#include <vector>
#include <cassert>
#include <boost/shared_ptr.hpp>
#include "core.h"

class CorrelationFunctionLorentzians;

namespace cfa {
	class ReducedOperatorApproximationLorentziansHigherOrder
	{
	public:
		typedef Eigen::VectorXcd data_type;
	public:
		CFA_API ReducedOperatorApproximationLorentziansHigherOrder(const Eigen::MatrixXcd& H0, boost::shared_ptr<const CorrelationFunctionLorentzians> alpha/*, double kBT, bool fixed_bath*/);
		//! Calculate dy_dt = f(y, t)
		CFA_API void operator()(double t, const data_type& y, data_type& dy_dt);
		CFA_API size_t dim() const { return m_dim; }
		CFA_API size_t operator_dim() const { return m_operator_dim; }
		CFA_API size_t nbr_sites() const { return m_nbr_sites; }
		CFA_API void init_from_pure_state(const Eigen::VectorXcd& exciton, data_type& y) const;
		CFA_API void exciton_density_matrix(const data_type& y, const Eigen::VectorXcd& initial_pure_state, Eigen::MatrixXcd& rho);
		CFA_API std::complex<double> exciton_density_matrix_trace(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const;
		CFA_API void reset_H0(const Eigen::MatrixXcd& H0);
		//! int_0^t <Psi(0)|H(t)|Psi(0)> dt
		//! not implemented
		CFA_API double energy_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const { return 0; }
		struct Normalizer
		{
			CFA_API Normalizer(const ReducedOperatorApproximationLorentziansHigherOrder&) {}
			CFA_API inline void normalize(data_type& state) const
			{
				; // do nothing
			}
		};
		typedef Normalizer normalizer_type;
		CFA_API static const bool WRAP_RHO_EIGENVALUES = false;
		CFA_API static const double RHO_EIGENVALUE_TOLERANCE;
		//! Normalize the state y at time t
		CFA_API void cfa_normalize(data_type& y, double t);
		static const char* name() { return "MFALorHO"; }
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
		operator_type bath_operator(data_type& y, size_t site_idx, size_t peak_idx) const 
		{ 
			return get_operator(y, bath_index(site_idx, peak_idx)); 
		}
		const_operator_type bath_operator(const data_type& y, size_t site_idx, size_t peak_idx) const 
		{ 
			return get_operator(y, bath_index(site_idx, peak_idx)); 
		}
		//! t_{mn} \bar{t}_{m'j}
		operator_type exciton_bath_operator(data_type& y, size_t m, size_t n, size_t m2, size_t j) const
		{
			return get_operator(y, exciton_bath_index(m, n, m2, j));
		}
		const_operator_type exciton_bath_operator(const data_type& y, size_t m, size_t n, size_t m2, size_t j) const
		{
			return get_operator(y, exciton_bath_index(m, n, m2, j));
		}
		operator_type work_operator(size_t idx) 
		{ 
			assert(idx < NBR_WORK_OPERATORS);
			return get_operator(m_work, idx*m_nbr_operator_elements);
		}
		size_t bath_index(size_t site_idx, size_t peak_idx) const 
		{ 
			assert(site_idx < m_nbr_sites);
			assert(peak_idx < m_nbr_peaks);
			return m_nbr_exciton_values + (site_idx*m_nbr_peaks + peak_idx)*m_nbr_operator_elements; 
		}
		size_t exciton_bath_index(size_t m, size_t n, size_t m2, size_t j) const 
		{ 
			assert(m < m_nbr_sites);
			assert(n < m_nbr_sites);
			assert(m2 < m_nbr_sites);
			assert(j < m_nbr_peaks);
			return m_exciton_bath_offset + (((m*m_nbr_sites + n)*m_nbr_sites + m2)*m_nbr_peaks + j)*m_nbr_operator_elements;
		}
		void init_exciton_part(data_type& y) const;
		template <class M> std::complex<double> mean_value(const M& op, const Eigen::VectorXcd& initial_pure_state) const;
		//! sum of bath operators \bar{t}_{nj}(t) over sites n is equal to identity operator times this constant
		std::complex<double> sum_bath_constant(size_t j, double t) const;		
	private:
		Eigen::MatrixXcd m_iH0;
		/*double m_kBT;
		bool m_fixed_bath;*/
		Eigen::VectorXd m_scales;
		Eigen::VectorXd m_hwhms;
		Eigen::VectorXd m_omegas;
		static const size_t NBR_WORK_OPERATORS = 2;
		Eigen::VectorXcd m_work;
		Eigen::VectorXcd m_multipliers;
		struct cfa_normalizer_workspace
		{
			cfa_normalizer_workspace(size_t nbr_sites);
			Eigen::MatrixXcd m_A;
			Eigen::VectorXcd m_b;
			Eigen::VectorXcd m_x;
			Eigen::JacobiSVD<Eigen::MatrixXcd> m_solver;
		};
		cfa_normalizer_workspace m_cfa_norm_wksp;
		size_t m_operator_dim;
		size_t m_nbr_sites;
		size_t m_dim;
		size_t m_nbr_operator_elements;
		size_t m_nbr_exciton_values;
		size_t m_nbr_peaks;		
		size_t m_nbr_bath_operators;
		size_t m_exciton_bath_offset;
	};
}

#endif // __CFA_ROA_LHO_FULL_H
