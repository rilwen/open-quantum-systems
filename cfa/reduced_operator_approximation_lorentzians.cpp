#include "reduced_operator_approximation_lorentzians.h"
#include "cfa_utils.h"
#include <protein_chain/correlation_function_lorentzians.h>
#include <cassert>
#include <stdexcept>

namespace cfa {
	ReducedOperatorApproximationLorentzians::ReducedOperatorApproximationLorentzians(const Eigen::MatrixXcd& H0, boost::shared_ptr<const CorrelationFunctionLorentzians> alpha/*, double kBT, bool fixed_bath*/)
		: m_H0(H0)/*, m_kBT(kBT), m_fixed_bath(fixed_bath)*/, m_scales(alpha->scales()), m_hwhms(alpha->gammas()), m_omegas(alpha->omegas()), m_multipliers(alpha->scales().size())
	{
		m_nbr_sites = H0.rows();
		if (!m_nbr_sites) {
			throw std::domain_error("No sites");
		}
		m_operator_dim = m_nbr_sites;
		assert(H0.cols() == H0.rows());
		m_nbr_operator_elements = m_operator_dim*m_operator_dim;
		m_nbr_exciton_values = m_nbr_operator_elements*(m_nbr_sites*(m_nbr_sites + 1))/2;
		m_nbr_peaks = alpha->gammas().size();
		m_nbr_bath_operators = m_nbr_sites * m_nbr_peaks;
		m_dim = m_nbr_exciton_values + m_nbr_bath_operators*m_nbr_operator_elements;
		m_work.setZero((m_nbr_sites + 1)*m_nbr_operator_elements);
		for (size_t i = 0; i < m_nbr_peaks; ++i) {
			m_multipliers[i] = std::complex<double>(-m_hwhms[i], -m_omegas[i]);
		}
	}

	ReducedOperatorApproximationLorentzians::operator_type ReducedOperatorApproximationLorentzians::exciton_operator(data_type& y, size_t m, size_t n) const
	{
		return get_operator(y, m_nbr_operator_elements*exciton_index(m, n));
	}

	ReducedOperatorApproximationLorentzians::const_operator_type ReducedOperatorApproximationLorentzians::exciton_operator(const data_type& y, size_t m, size_t n) const
	{
		return get_operator(y, m_nbr_operator_elements*exciton_index(m, n));
	}	

	void ReducedOperatorApproximationLorentzians::operator()(double t, const data_type& y, data_type& dy_dt)
	{
		m_work.setZero();		
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			operator_type work(work_operator(m));
			const_operator_type tmm(exciton_operator(y, m, m));
//#ifndef NDEBUG
//			std::cout << t << " " << m << " " << (tmm - tmm.adjoint()).norm() << std::endl;
//#endif
			for (size_t j = 0; j < m_nbr_peaks; ++j) {
				work += bath_operator(y, m, j) * m_scales[j];
				operator_type op(bath_operator(dy_dt, m, j));
				op = bath_operator(y, m, j) * m_multipliers[j];
				op += tmm;
			}
			const Eigen::MatrixXcd::ColXpr H0m(m_H0.col(m));
			for (size_t n = 0; n <= m; ++n) {
				operator_type dtmn_dt(exciton_operator(dy_dt, m, n));								
				if (n < m) {
					operator_type work2(work_operator(m_nbr_sites));
					work2 = work;
					work2 -= work_operator(n);
					cfa::subtract_adjoint_in_place(work2, m_operator_dim);
					assert(work2 == -work2.adjoint());
					// symmetrized product, because bath and exciton matrices do not commute in MFA approximation
					dtmn_dt.noalias() = exciton_operator(y, m, n) * work2;
					dtmn_dt.noalias() += work2 * exciton_operator(y, m, n);
					dtmn_dt *= -0.5*IMAGINARY;
				} else {
					dtmn_dt.setZero();
				}
				const Eigen::MatrixXcd::RowXpr H0n(m_H0.row(n));
				size_t m2 = 0;
				assert(n <= m);
				for (; m2 < n; ++m2) {
					assert(m2 < n);
					assert(m2 < m);
					dtmn_dt += exciton_operator(y, n, m2).adjoint() * H0m[m2];
					dtmn_dt -= exciton_operator(y, m, m2) * H0n[m2];
				}
				for (; m2 < m; ++m2) {
					assert(m2 >= n);
					assert(m2 < m);
					dtmn_dt += exciton_operator(y, m2, n) * H0m[m2];
					dtmn_dt -= exciton_operator(y, m, m2) * H0n[m2];
				}
				for (; m2 < m_nbr_sites; ++m2) {
					assert(m2 >= n);
					assert(m2 >= m);
					dtmn_dt += exciton_operator(y, m2, n) * H0m[m2];
					dtmn_dt -= exciton_operator(y, m2, m).adjoint() * H0n[m2];
				}				
				dtmn_dt *= IMAGINARY;
			}
		}
	}

	void ReducedOperatorApproximationLorentzians::init_exciton_part(data_type& y) const
	{
		y.resize(m_dim);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			for (size_t n = 0; n <= m; ++n) {
				operator_type op(exciton_operator(y, m, n));
				op.setZero();				
				op(m, n) = 1.0;
			}
		}
	}

	void ReducedOperatorApproximationLorentzians::init_from_pure_state(const Eigen::VectorXcd&, data_type& y) const
	{		
		y.setZero(m_dim);
		init_exciton_part(y);
		assert(y.size() == m_dim);	
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			for (size_t j = 0; j < m_nbr_peaks; ++j) {
				bath_operator(y, m, j).setZero();
			}
		}
	}	

	template <class M> std::complex<double> ReducedOperatorApproximationLorentzians::mean_value(const M& op, const Eigen::VectorXcd& initial_pure_state) const
	{
		return initial_pure_state.dot(op*initial_pure_state);
	}

	void ReducedOperatorApproximationLorentzians::exciton_density_matrix(const data_type& y, const Eigen::VectorXcd& initial_pure_state, Eigen::MatrixXcd& rho)
	{
		/*rho.resize(m_nbr_sites, m_nbr_sites);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			for (size_t n = 0; n <= m; ++n) {
				const_operator_type op(exciton_operator(y, m, n));
				rho(m, n) = initial_pure_state.dot(op*initial_pure_state);
				rho(n, m) = conj(rho(m, n));
			}
		}*/
		rho.setZero(m_nbr_sites, m_nbr_sites);
		operator_type tmp(work_operator(0));
		for (size_t i = 0; i < m_nbr_sites; ++i) {
			for (size_t m = 0; m < m_nbr_sites; ++m) {
				for (size_t n = 0; n < m_nbr_sites; ++n) {
					if (m >= i) {
						if (n > i) {
							tmp.noalias() = exciton_operator(y, m, i) * exciton_operator(y, n, i).adjoint();
						} else {
							tmp.noalias() = exciton_operator(y, m, i) * exciton_operator(y, i, n);
						}
					} else {
						if (n > i) {
							tmp.noalias() = exciton_operator(y, i, m).adjoint() * exciton_operator(y, n, i).adjoint();
						} else {
							tmp.noalias() = exciton_operator(y, i, m).adjoint() * exciton_operator(y, i, n);
						}
					}
					rho(m, n) += mean_value(tmp, initial_pure_state);
				}
			}
		}
		rho /= m_nbr_sites;
	}

	std::complex<double> ReducedOperatorApproximationLorentzians::exciton_density_matrix_trace(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{
		std::complex<double> sum(0.0);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			sum += mean_value(exciton_operator(y, m, m), initial_pure_state);
		}
		return sum;
	}

	std::complex<double> ReducedOperatorApproximationLorentzians::sum_bath_constant(size_t j, double t) const
	{		
		if (m_multipliers[j].real() == 0.0 && m_multipliers[j].imag() == 0.0) {
			return t * m_scales[j];
		} else {
			std::complex<double> result(0.0);
			result = exp(m_multipliers[j] * t);
			result -= 1.0;
			result *= m_scales[j];
			result /= m_multipliers[j];
			return result;
		}
	}

	const double ReducedOperatorApproximationLorentzians::RHO_EIGENVALUE_TOLERANCE = 1E-3;
}
