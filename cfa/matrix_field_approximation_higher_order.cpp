#include "matrix_field_approximation_higher_order.h"
#include "cfa_utils.h"
#include <cassert>
#include <stdexcept>
#include <unsupported/Eigen/MatrixFunctions>
#include <protein_chain/thermal_states.h>

namespace cfa {
	MatrixFieldApproximationHigherOrder::MatrixFieldApproximationHigherOrder(const Eigen::MatrixXcd& H0, const Eigen::MatrixXcd& g, const Eigen::VectorXd& omega, double kBT, bool fixed_bath)
		: m_H0(H0), m_g(g), m_omega(omega), m_kBT(kBT), m_fixed_bath(fixed_bath)
	{
		m_nbr_sites = H0.rows();
		if (!m_nbr_sites) {
			throw std::domain_error("No sites");
		}
		m_operator_dim = m_nbr_sites;
		m_work.resize(m_operator_dim, m_operator_dim);
		m_work2.resize(m_operator_dim, m_operator_dim);
		assert(H0.cols() == H0.rows());
		m_nbr_bath_modes = omega.size();
		if (g.rows() != m_nbr_bath_modes) {
			throw std::domain_error("g.rows() != m_nbr_bath_modes");
		}
		if (g.cols() != m_nbr_sites) {
			throw std::domain_error("g.cols() != m_nbr_sites");
		}
		m_nbr_operator_elements = m_operator_dim*m_operator_dim;
		const size_t nbr_exciton_operators = (m_nbr_sites*(m_nbr_sites + 1))/2;
		m_nbr_exciton_values = m_nbr_operator_elements*nbr_exciton_operators;
		m_nbr_exciton_values_both_halves = m_nbr_sites*m_nbr_sites*m_nbr_operator_elements;
		m_exciton_bath_product_offset = m_nbr_exciton_values + m_nbr_bath_modes*m_nbr_operator_elements;
		m_hamiltonian_offset = m_exciton_bath_product_offset + m_nbr_bath_modes*m_nbr_exciton_values_both_halves;
		m_dim = m_hamiltonian_offset + m_nbr_operator_elements;
	}

	MatrixFieldApproximationHigherOrder::operator_type MatrixFieldApproximationHigherOrder::exciton_operator(data_type& y, size_t m, size_t n) const
	{
		return get_operator(y, m_nbr_operator_elements*exciton_index(m, n));
	}

	MatrixFieldApproximationHigherOrder::const_operator_type MatrixFieldApproximationHigherOrder::exciton_operator(const data_type& y, size_t m, size_t n) const
	{
		return get_operator(y, m_nbr_operator_elements*exciton_index(m, n));
	}

	MatrixFieldApproximationHigherOrder::operator_type MatrixFieldApproximationHigherOrder::bath_operator(data_type& y, size_t k) const
	{
		return get_operator(y, bath_index(k));
	}

	MatrixFieldApproximationHigherOrder::const_operator_type MatrixFieldApproximationHigherOrder::bath_operator(const data_type& y, size_t k) const
	{
		return get_operator(y, bath_index(k));
	}

	MatrixFieldApproximationHigherOrder::operator_type MatrixFieldApproximationHigherOrder::hamiltonian_operator(data_type& y) const
	{
		return get_operator(y, m_hamiltonian_offset);
	}

	MatrixFieldApproximationHigherOrder::const_operator_type MatrixFieldApproximationHigherOrder::hamiltonian_operator(const data_type& y) const
	{
		return get_operator(y, m_hamiltonian_offset);
	}

	size_t MatrixFieldApproximationHigherOrder::exciton_bath_index(size_t m, size_t n, size_t k) const
	{
		return m_exciton_bath_product_offset + k*m_nbr_exciton_values_both_halves + m_nbr_operator_elements*(m*m_nbr_sites + n);
	}

	MatrixFieldApproximationHigherOrder::operator_type MatrixFieldApproximationHigherOrder::exciton_bath_operator(data_type& y, size_t m, size_t n, size_t k) const
	{
		return get_operator(y, exciton_bath_index(m, n, k));
	}

	MatrixFieldApproximationHigherOrder::const_operator_type MatrixFieldApproximationHigherOrder::exciton_bath_operator(const data_type& y, size_t m, size_t n, size_t k) const
	{
		return get_operator(y, exciton_bath_index(m, n, k));
	}

	void MatrixFieldApproximationHigherOrder::init_exciton_part(data_type& y) const
	{
		y.resize(m_dim);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			for (size_t n = 0; n <= m; ++n) {
				operator_type op(exciton_operator(y, m, n));
				op.setZero();
				op(m, n) = 1.0;
			}
		}
		hamiltonian_operator(y).setZero();
	}

	void MatrixFieldApproximationHigherOrder::init_from_pure_state(const Eigen::VectorXcd&, data_type& y) const
	{		
		init_exciton_part(y);
		assert(y.size() == m_dim);
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			operator_type bk(bath_operator(y, k));
			const double sqrt_avgocc = m_kBT != 0 ? sqrt(ThermalStates::average_mode_occupation(m_omega[k], m_kBT)) : 0.0;
			if (m_kBT == 0) {
				bk.setZero();
			} else {				
				bk.setIdentity();
				bk *= sqrt_avgocc;
			}
			for (size_t m = 0; m < m_nbr_sites; ++m) {
				for (size_t n = 0; n < m_nbr_sites; ++n) {
					operator_type smnk(exciton_bath_operator(y, m, n, k));
					smnk.setZero();
					if (m_kBT != 0) {
						smnk(m, n) = sqrt_avgocc;
					}
				}
			}
		}
	}

	void MatrixFieldApproximationHigherOrder::init_from_pure_state(const Eigen::VectorXcd&, const Eigen::VectorXcd& bath, data_type& y) const
	{
		init_exciton_part(y);
		assert(y.size() == m_dim);
		const unsigned int ratio = m_operator_dim / m_nbr_sites;
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			operator_type op(bath_operator(y, k));
			op.setIdentity();
			op *= bath[k];
			for (size_t m = 0; m < m_nbr_sites; ++m) {
				for (size_t n = 0; n <= m; ++n) {
					exciton_bath_operator(y, m, n, k) = exciton_operator(y, m, n) * op;
				}
				for (size_t n = m + 1; n < m_nbr_sites; ++n) {
					exciton_bath_operator(y, m, n, k) = exciton_operator(y, m, n).adjoint() * op;
				}
			}
		}
	}

	template <class M> std::complex<double> MatrixFieldApproximationHigherOrder::mean_value(const M& op, const Eigen::VectorXcd& initial_pure_state) const
	{
		return initial_pure_state.dot(op*initial_pure_state);
	}

	void MatrixFieldApproximationHigherOrder::exciton_density_matrix(const data_type& y, const Eigen::VectorXcd& initial_pure_state, Eigen::MatrixXcd& rho)
	{
		rho.setZero(m_nbr_sites, m_nbr_sites);
		for (size_t i = 0; i < m_nbr_sites; ++i) {
			for (size_t m = 0; m < m_nbr_sites; ++m) {
				for (size_t n = 0; n < m_nbr_sites; ++n) {
					if (m >= i) {
						if (n >= i) {
							m_work.noalias() = exciton_operator(y, m, i) * exciton_operator(y, n, i).adjoint();
						} else {
							m_work.noalias() = exciton_operator(y, m, i) * exciton_operator(y, i, n);
						}
					} else {
						if (n >= i) {
							m_work.noalias() = exciton_operator(y, i, m).adjoint() * exciton_operator(y, n, i).adjoint();
						} else {
							m_work.noalias() = exciton_operator(y, i, m).adjoint() * exciton_operator(y, i, n);
						}
					}
					rho(m, n) += mean_value(m_work, initial_pure_state); // TODO: make it work for T > 0 or simulations with bath out of ground state at the beginning					
				}
			}
		}
		rho /= m_nbr_sites;
	}

	std::complex<double> MatrixFieldApproximationHigherOrder::exciton_density_matrix_trace(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{
		std::complex<double> sum(0.0);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			sum += mean_value(exciton_operator(y, m, m), initial_pure_state);
		}
		return sum;
	}

	void MatrixFieldApproximationHigherOrder::operator()(double t, const data_type& y, data_type& dy_dt)
	{
		assert(y.size() == m_dim);
		assert(dy_dt.size() == m_dim);
		dy_dt.setZero();
		
		if (m_fixed_bath) {
			for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
				bath_operator(dy_dt, k).setZero();
			}
		} else {
			// calculate db_k(t)/dt
			for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
				operator_type sum(bath_operator(dy_dt, k));
				sum = bath_operator(y, k);
				sum *= m_omega[k];
				const Eigen::MatrixXcd::RowXpr gk(m_g.row(k));
				for (size_t m = 0; m < m_nbr_sites; ++m) {
					if (gk[m] != 0.0) {
						sum.noalias() += exciton_operator(y, m, m) * gk[m];
					}
				}
				sum *= -IMAGINARY;
			}
		}

		// calculate d|m><n|(t)/dt
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			const Eigen::MatrixXcd::ColXpr H0m(m_H0.col(m));
			const Eigen::MatrixXcd::ColXpr gm(m_g.col(m));
			for (size_t n = 0; n <= m; ++n) {
				operator_type sum(exciton_operator(dy_dt, m, n));
				sum.setZero();				
				if (m != n) {					
					const Eigen::MatrixXcd::ColXpr gn(m_g.col(n));
					for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
						const std::complex<double> dg = gm[k] - gn[k];
						if (dg.real() != 0 || dg.imag() != 0) {
							sum += exciton_bath_operator(y, m, n, k) * conj(dg);
							sum += exciton_bath_operator(y, n, m, k).adjoint() * dg;
						}
					}
				}				
				const Eigen::MatrixXcd::RowXpr H0n(m_H0.row(n));
				size_t m2 = 0;
				assert(n <= m);
				for (; m2 < n; ++m2) {
					assert(m2 < n);
					assert(m2 < m);
					sum += exciton_operator(y, n, m2).adjoint() * H0m[m2];
					sum -= exciton_operator(y, m, m2) * H0n[m2];
				}
				for (; m2 < m; ++m2) {
					assert(m2 >= n);
					assert(m2 < m);
					sum += exciton_operator(y, m2, n) * H0m[m2];
					sum -= exciton_operator(y, m, m2) * H0n[m2];
				}
				for (; m2 < m_nbr_sites; ++m2) {
					assert(m2 >= n);
					assert(m2 >= m);
					sum += exciton_operator(y, m2, n) * H0m[m2];
					sum -= exciton_operator(y, m2, m).adjoint() * H0n[m2];
				}				
				sum *= IMAGINARY;
				assert(m != n || cfa::is_hermitian(sum));
			}
		}

		// calculate ds_{mnk}(t)/dt
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			const std::complex<double>& w_k = m_omega[k];
			const Eigen::MatrixXcd::RowXpr gk(m_g.row(k));
			const_operator_type bk(bath_operator(y, k));
			for (size_t m = 0; m < m_nbr_sites; ++m) {
				const Eigen::MatrixXcd::ColXpr H0m(m_H0.col(m));
				const Eigen::MatrixXcd::ColXpr gm(m_g.col(m));
				for (size_t n = 0; n < m_nbr_sites; ++n) {
					operator_type sum(exciton_bath_operator(dy_dt, m, n, k));
					if (m_fixed_bath) {
						sum.setZero();
					} else {
						sum = exciton_bath_operator(y, m, n, k) * (-w_k);
						const std::complex<double> a = -gk[n];
						if (a.real() != 0.0 || a.imag() != 0.0) {
							if (n <= m) {
								sum += exciton_operator(y, m, n) * a;
							} else {
								sum += exciton_operator(y, n, m).adjoint() * a;
							}
						}
					}

					if (m != n) {
						const Eigen::MatrixXcd::ColXpr gn(m_g.col(n));

						// symmetrise indices

						m_work.setZero();
						for (size_t k2 = 0; k2 < m_nbr_bath_modes; ++k2) {
							const std::complex<double> dgh = (gm[k2] - gn[k2]) / 2.0;
							if (dgh.real() != 0 || dgh.imag() != 0) {
								m_work += exciton_bath_operator(y, m, n, k2) * conj(dgh);
								m_work += exciton_bath_operator(y, n, m, k2).adjoint() * dgh;
							}
						}
						sum.noalias() += m_work * bk;

						m_work.setZero();
						for (size_t k2 = 0; k2 < m_nbr_bath_modes; ++k2) {
							const std::complex<double> dgh = (gm[k2] - gn[k2]) / 2.0;
							if (dgh.real() != 0 || dgh.imag() != 0) {
								m_work += bath_operator(y, k2) * conj(dgh);
							}
						}
						sum.noalias() += m_work.adjoint() * exciton_bath_operator(y, m, n, k);
						sum.noalias() += exciton_bath_operator(y, m, n, k) * m_work;
					}

					const Eigen::MatrixXcd::RowXpr H0n(m_H0.row(n));
					for (size_t m2 = 0; m2 < m_nbr_sites; ++m2) {
						sum += H0m[m2] * exciton_bath_operator(y, m2, n, k);
						sum -= H0n[m2] * exciton_bath_operator(y, m, m2, k);
					}

					sum *= IMAGINARY;
				}
			}
		}

		operator_type dham(dy_dt.data() + m_hamiltonian_offset, m_operator_dim, m_operator_dim);
		dham.setZero();
		//// calculate sum_{m,n} H0(m,n)|m><n| + sum_m |m><m| sum_k ( \tilde{b}_k(t) e^{-i \omega_k t} \cc{g}_{km} + c.c. )
		//for (size_t m = 0; m < m_nbr_sites; ++m) {
		//	const Eigen::MatrixXcd::RowXpr H0m(m_H0.row(m));
		//	const_operator_type op_mm(exciton_operator(y, m, m));
		//	assert(cfa::is_hermitian(op_mm));
		//	assert(H0m[m].imag() == 0);
		//	dham += op_mm * H0m[m];
		//	assert(cfa::is_hermitian(dham));
		//	const Eigen::MatrixXcd::ColXpr gm(m_g.col(m));
		//	for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
		//		if (gm[k] != 0.0) {
		//			dham += exciton_bath_operator(y, m, m, k) * conj(gm[k]);
		//			dham += exciton_bath_operator(y, m, m, k).adjoint() * gm[k];
		//		}
		//	}
		//	cfa::make_hermitian(dham); // this is necessary because op_mm and m_work_g[m] do not have to commute in the MFA approximation
		//	for (size_t n = 0; n < m; ++n) {
		//		const_operator_type op(exciton_operator(y, m, n));
		//		dham += op * H0m[n];
		//		dham += op.adjoint() * conj(H0m[n]);
		//	}			
		//	assert(cfa::is_hermitian(dham));			
		//}
		//// calculate sum_k omega_k b^dagger_k b_k
		//for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
		//	const_operator_type bk(bath_operator(y, k));
		//	m_work.noalias() = bk.adjoint() * bk;
		//	dham += m_work * m_omega[k];
		//}
		//assert(cfa::is_hermitian(dham));
	}

	void MatrixFieldApproximationHigherOrder::reset_H0(const Eigen::MatrixXcd& H0)
	{
		m_H0 = H0;
	}

	double MatrixFieldApproximationHigherOrder::energy_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{
		const_operator_type H(hamiltonian_operator(y));
		return mean_value(H, initial_pure_state).real();
	}

	std::complex<double> MatrixFieldApproximationHigherOrder::evolution_operator_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{		
		Eigen::MatrixXcd minus_hamiltonian_operator(m_operator_dim, m_operator_dim);
		minus_hamiltonian_operator = hamiltonian_operator(y);
		minus_hamiltonian_operator *= -IMAGINARY;
		Eigen::MatrixExponential<Eigen::MatrixXcd> exponential(minus_hamiltonian_operator); // evolution operator
		Eigen::MatrixXcd U(m_operator_dim, m_operator_dim);
		exponential.compute(U);
		return initial_pure_state.dot(U*initial_pure_state);
		return mean_value(U, initial_pure_state);
	}

	const double MatrixFieldApproximationHigherOrder::RHO_EIGENVALUE_TOLERANCE = 1E-3;
}

