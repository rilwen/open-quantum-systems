#include <Eigen/Dense>
#include "reduced_operator_approximation_lorentzians_higher_order.h"
#include "cfa_utils.h"
#include <protein_chain/correlation_function_lorentzians.h>
#include <stdexcept>

namespace cfa {
	ReducedOperatorApproximationLorentziansHigherOrder::ReducedOperatorApproximationLorentziansHigherOrder(const Eigen::MatrixXcd& H0, boost::shared_ptr<const CorrelationFunctionLorentzians> alpha/*, double kBT, bool fixed_bath*/)
		: m_iH0(H0)/*, m_kBT(kBT), m_fixed_bath(fixed_bath)*/, m_scales(alpha->scales()), m_hwhms(alpha->gammas()), m_omegas(alpha->omegas()), m_multipliers(alpha->scales().size()), m_cfa_norm_wksp(H0.rows())
	{
		m_iH0 *= IMAGINARY;
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
		m_exciton_bath_offset = m_nbr_exciton_values + m_nbr_bath_operators*m_nbr_operator_elements;
		const size_t nbr_bathexciton_set_values = m_nbr_sites*m_nbr_sites*m_nbr_bath_operators*m_nbr_operator_elements;
		m_dim = m_exciton_bath_offset + nbr_bathexciton_set_values;
		m_work.setZero(m_nbr_operator_elements * NBR_WORK_OPERATORS);
		for (size_t i = 0; i < m_nbr_peaks; ++i) {
			m_multipliers[i] = std::complex<double>(-m_hwhms[i], -m_omegas[i]);
		}
	}

	ReducedOperatorApproximationLorentziansHigherOrder::operator_type ReducedOperatorApproximationLorentziansHigherOrder::exciton_operator(data_type& y, size_t m, size_t n) const
	{
		assert(m < m_nbr_peaks);
		assert(n <= m);
		return get_operator(y, m_nbr_operator_elements*exciton_index(m, n));
	}

	ReducedOperatorApproximationLorentziansHigherOrder::const_operator_type ReducedOperatorApproximationLorentziansHigherOrder::exciton_operator(const data_type& y, size_t m, size_t n) const
	{
		assert(m < m_nbr_peaks);
		assert(n <= m);
		return get_operator(y, m_nbr_operator_elements*exciton_index(m, n));
	}	

	void ReducedOperatorApproximationLorentziansHigherOrder::operator()(double t, const data_type& y, data_type& dy_dt)
	{
		m_work.setZero();
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			const_operator_type tmm(exciton_operator(y, m, m));
			// do bath operators
			for (size_t j = 0; j < m_nbr_peaks; ++j) {
				operator_type op(bath_operator(dy_dt, m, j));
				op = bath_operator(y, m, j) * m_multipliers[j];
				op += tmm * m_scales[j];
			}
			const Eigen::MatrixXcd::ColXpr iH0m(m_iH0.col(m));

			exciton_operator(dy_dt, m, m).setZero();

			for (size_t n = 0; n < m_nbr_sites; ++n) {
				const Eigen::MatrixXcd::RowXpr iH0n(m_iH0.row(n));

				for (size_t m1 = 0; m1 < m_nbr_sites; ++m1) {
					for (size_t j = 0; j < m_nbr_peaks; ++j) {
						operator_type dsmnm1j_dt(exciton_bath_operator(dy_dt, m, n, m1, j));

						// oscillatory part
						dsmnm1j_dt = exciton_bath_operator(y, m, n, m1, j) * m_multipliers[j]; 

						// system hamiltonian part
						for (size_t m2 = 0; m2 < m_nbr_sites; ++m2) {
							dsmnm1j_dt -= exciton_bath_operator(y, m, m2, m1, j) * iH0n[m2];
							dsmnm1j_dt += exciton_bath_operator(y, m2, n, m1, j) * iH0m[m2];
						} // m2
					} // j
				} // m1

				if (n != m) {
					operator_type work0(work_operator(0));
					work0.setZero();
					for (size_t j = 0; j < m_nbr_peaks; ++j) {
						work0 += exciton_bath_operator(y, m, n, m, j);
						work0 += exciton_bath_operator(y, n, m, n, j).adjoint();
						work0 -= exciton_bath_operator(y, m, n, n, j);
						work0 -= exciton_bath_operator(y, n, m, m, j).adjoint();
					} // j
					if (n < m) {
						exciton_operator(dy_dt, m, n) = work0;
						const_operator_type tmn(exciton_operator(y, m, n));
						for (size_t j = 0; j < m_nbr_peaks; ++j) {
							exciton_bath_operator(dy_dt, m, n, n, j) += tmn * (m_scales[j]);
						} // j
					} /* n < m */ else {
						assert(n < m);
						for (size_t j = 0; j < m_nbr_peaks; ++j) {
							const_operator_type tnm(exciton_operator(y, n, m));
							exciton_bath_operator(dy_dt, m, n, n, j) += tnm.adjoint() * (m_scales[j]);
						} // j
					} // n > m

					work0 /= 2.0; // divide by 2 due to symmetrisation of the choice of bath operators
					operator_type work1(work_operator(1));
					for (size_t m1 = 0; m1 < m_nbr_sites; ++m1) {
						for (size_t j = 0; j < m_nbr_peaks; ++j) {
							work1.noalias() = work0 * bath_operator(y, m1, j);
							exciton_bath_operator(dy_dt, m, n, m1, j) += work1;
						} // j
					} // m1

					work0.setZero();
					for (size_t j = 0; j < m_nbr_peaks; ++j) {
						work0 += bath_operator(y, m, j);
						work0 -= bath_operator(y, n, j);
					} // j
					work0 /= 2.0; // divide by 2 due to symmetrisation of the choice of bath operators
					for (size_t m1 = 0; m1 < m_nbr_sites; ++m1) {
						for (size_t j = 0; j < m_nbr_peaks; ++j) {
							work1.noalias() = work0 * exciton_bath_operator(y, m, n, m1, j);
							exciton_bath_operator(dy_dt, m, n, m1, j) += work1;
							work1.noalias() = work0.adjoint() * exciton_bath_operator(y, m, n, m1, j);
							exciton_bath_operator(dy_dt, m, n, m1, j) -= work1;
						} // j
					} // m1
				} // n != m
				else {					
					for (size_t j = 0; j < m_nbr_peaks; ++j) {
						exciton_bath_operator(dy_dt, m, m, m, j) += tmm * (m_scales[j]);
					} // j
				} // n == m
			} // n

			// t_{mn} operators
			for (size_t n = 0; n <= m; ++n) {
				operator_type dtmn_dt(exciton_operator(dy_dt, m, n));								
				const Eigen::MatrixXcd::RowXpr iH0n(m_iH0.row(n));
				size_t m2 = 0;
				assert(n <= m);
				for (; m2 < n; ++m2) {
					assert(m2 < n);
					assert(m2 < m);
					dtmn_dt += exciton_operator(y, n, m2).adjoint() * iH0m[m2];
					dtmn_dt -= exciton_operator(y, m, m2) * iH0n[m2];
				}
				for (; m2 < m; ++m2) {
					assert(m2 >= n);
					assert(m2 < m);
					dtmn_dt += exciton_operator(y, m2, n) * iH0m[m2];
					dtmn_dt -= exciton_operator(y, m, m2) * iH0n[m2];
				}
				for (; m2 < m_nbr_sites; ++m2) {
					assert(m2 >= n);
					assert(m2 >= m);
					dtmn_dt += exciton_operator(y, m2, n) * iH0m[m2];
					dtmn_dt -= exciton_operator(y, m2, m).adjoint() * iH0n[m2];
				}
			}
		}
	}

	void ReducedOperatorApproximationLorentziansHigherOrder::init_exciton_part(data_type& y) const
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

	void ReducedOperatorApproximationLorentziansHigherOrder::init_from_pure_state(const Eigen::VectorXcd&, data_type& y) const
	{		
		init_exciton_part(y);
		assert(y.size() == m_dim);	
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			for (size_t j = 0; j < m_nbr_peaks; ++j) {
				bath_operator(y, m, j).setZero();
				for (size_t n = 0; n < m_nbr_sites; ++n) {
					for (size_t m2 = 0; m2 < m_nbr_sites; ++m2) {
						exciton_bath_operator(y, m, n, m2, j).setZero();
					}
				}
			}
		}
	}	

	template <class M> std::complex<double> ReducedOperatorApproximationLorentziansHigherOrder::mean_value(const M& op, const Eigen::VectorXcd& initial_pure_state) const
	{
		return initial_pure_state.dot(op*initial_pure_state);
	}

	void ReducedOperatorApproximationLorentziansHigherOrder::exciton_density_matrix(const data_type& y, const Eigen::VectorXcd& initial_pure_state, Eigen::MatrixXcd& rho)
	{
		rho.setZero(m_nbr_sites, m_nbr_sites);
		operator_type tmp(work_operator(0));
		for (size_t i = 0; i < m_nbr_sites; ++i) {
			for (size_t m = 0; m < m_nbr_sites; ++m) {
				for (size_t n = 0; n < m_nbr_sites; ++n) {
					if (m >= i) {
						if (n >= i) {
							tmp.noalias() = exciton_operator(y, m, i) * exciton_operator(y, n, i).adjoint();
						} else {
							tmp.noalias() = exciton_operator(y, m, i) * exciton_operator(y, i, n);
						}
					} else {
						if (n >= i) {
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

	std::complex<double> ReducedOperatorApproximationLorentziansHigherOrder::exciton_density_matrix_trace(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{
		std::complex<double> sum(0.0);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			sum += mean_value(exciton_operator(y, m, m), initial_pure_state);
		}
		return sum;
	}

	std::complex<double> ReducedOperatorApproximationLorentziansHigherOrder::sum_bath_constant(size_t j, double t) const
	{		
		if (cfa::is_zero(m_multipliers[j])) {
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

	void ReducedOperatorApproximationLorentziansHigherOrder::cfa_normalize(data_type& y, double t)
	{
		Eigen::VectorXcd& b = m_cfa_norm_wksp.m_b;
		Eigen::MatrixXcd& A = m_cfa_norm_wksp.m_A;
		Eigen::VectorXcd& x = m_cfa_norm_wksp.m_x;

		// normalize t_nj
		for (size_t j = 0; j < m_nbr_peaks; ++j) {
			const std::complex<double> expected_constant = sum_bath_constant(j, t);
			if (!cfa::is_zero(expected_constant)) {
				b.fill(expected_constant);
				for (size_t n1 = 0; n1 < m_nbr_sites; ++n1) {
					const_operator_type tn1j(bath_operator(static_cast<const data_type&>(y), n1, j));
					for (size_t n2 = 0; n2 < m_nbr_sites; ++n2) {
						A(n2, n1) = tn1j(n2, n2);
					}
				}
				//std::cout << t << "\t A == [" << A << "]\t b == " << m_cfa_norm_wksp.m_b << std::endl;
				m_cfa_norm_wksp.m_solver.compute(A, Eigen::ComputeThinU | Eigen::ComputeThinV);				
				x = m_cfa_norm_wksp.m_solver.solve(b);
				//std::cout << t << "\t A == [" << A << "]\t x == " << x << std::endl;
				for (size_t n1 = 0; n1 < m_nbr_sites; ++n1) {
					bath_operator(y, n1, j) *= x[n1];
				}
			}
		}
	}

	ReducedOperatorApproximationLorentziansHigherOrder::cfa_normalizer_workspace::cfa_normalizer_workspace(size_t nbr_sites)
		: m_A(nbr_sites, nbr_sites), m_b(nbr_sites), m_x(nbr_sites), m_solver(nbr_sites, nbr_sites)
	{
	}

	const double ReducedOperatorApproximationLorentziansHigherOrder::RHO_EIGENVALUE_TOLERANCE = 1E-3;
}

