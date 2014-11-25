#include "matrix_field_approximation.h"
#include "cfa_utils.h"
#include <cassert>
#include <stdexcept>
#include <unsupported/Eigen/MatrixFunctions>
#include <protein_chain/thermal_states.h>

namespace cfa {
	MatrixFieldApproximation::MatrixFieldApproximation(const Eigen::MatrixXcd& H0, const Eigen::MatrixXcd& g, const Eigen::VectorXd& omega, double kBT, bool fixed_bath)
		: m_H0(H0), m_g(g), m_omega(omega), m_exp_w(omega.size()), m_work_g(H0.rows()), m_kBT(kBT), m_fixed_bath(fixed_bath)
	{
		m_nbr_sites = H0.rows();
		if (!m_nbr_sites) {
			throw std::domain_error("No sites");
		}
		m_operator_dim = m_nbr_sites;
		assert(H0.cols() == H0.rows());
		m_nbr_bath_modes = omega.size();
		if (g.rows() != m_nbr_bath_modes) {
			throw std::domain_error("g.rows() != m_nbr_bath_modes");
		}
		if (g.cols() != m_nbr_sites) {
			throw std::domain_error("g.cols() != m_nbr_sites");
		}
		m_nbr_operator_elements = m_operator_dim*m_operator_dim;
		m_nbr_exciton_values = m_nbr_operator_elements*(m_nbr_sites*(m_nbr_sites + 1))/2;
		m_hamiltonian_offset = m_nbr_exciton_values + m_nbr_bath_modes*m_nbr_operator_elements;
		m_dim = m_hamiltonian_offset + m_nbr_operator_elements;
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			m_work_g[m].setZero(m_operator_dim, m_operator_dim);
		}
	}

	MatrixFieldApproximation::operator_type MatrixFieldApproximation::exciton_operator(data_type& y, size_t m, size_t n) const
	{
		return get_operator(y, m_nbr_operator_elements*exciton_index(m, n));
	}

	MatrixFieldApproximation::const_operator_type MatrixFieldApproximation::exciton_operator(const data_type& y, size_t m, size_t n) const
	{
		return get_operator(y, m_nbr_operator_elements*exciton_index(m, n));
	}

	MatrixFieldApproximation::operator_type MatrixFieldApproximation::bath_operator(data_type& y, size_t k) const
	{
		return get_operator(y, bath_index(k));
	}

	MatrixFieldApproximation::const_operator_type MatrixFieldApproximation::bath_operator(const data_type& y, size_t k) const
	{
		return get_operator(y, bath_index(k));
	}

	MatrixFieldApproximation::operator_type MatrixFieldApproximation::hamiltonian_operator(data_type& y) const
	{
		return get_operator(y, m_hamiltonian_offset);
	}

	MatrixFieldApproximation::const_operator_type MatrixFieldApproximation::hamiltonian_operator(const data_type& y) const
	{
		return get_operator(y, m_hamiltonian_offset);
	}

	void MatrixFieldApproximation::init_exciton_part(data_type& y) const
	{
		y.resize(m_dim);
		//const unsigned int ratio = m_operator_dim / m_nbr_sites;
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			for (size_t n = 0; n <= m; ++n) {
				operator_type op(exciton_operator(y, m, n));
				op.setZero();				
				op(m, n) = 1.0;
				/*for (unsigned int r = 0; r < ratio; ++r) {
					op(m + r*m_nbr_sites, n + r*m_nbr_sites) = 1.0;
				}*/
			}
		}
		hamiltonian_operator(y).setZero();
	}

	void MatrixFieldApproximation::init_from_pure_state(const Eigen::VectorXcd&, data_type& y) const
	{		
		init_exciton_part(y);
		assert(y.size() == m_dim);
		/*const unsigned int ratio = m_operator_dim / m_nbr_sites;*/
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			operator_type bk(bath_operator(y, k));
			if (m_kBT == 0) {
				bk.setZero();
			} else {
				bk.setIdentity();
				bk *= sqrt(ThermalStates::average_mode_occupation(m_omega[k], m_kBT));
			}
			/*for (size_t r = 1; r < ratio; ++r) {
				bk.block((r-1)*m_nbr_sites, r*m_nbr_sites, m_nbr_sites, m_nbr_sites).setIdentity();
			}*/
		}
	}

	void MatrixFieldApproximation::init_from_pure_state(const Eigen::VectorXcd&, const Eigen::VectorXcd& bath, data_type& y) const
	{
		init_exciton_part(y);
		assert(y.size() == m_dim);
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			operator_type op(bath_operator(y, k));
			op.setIdentity();
			op *= bath[k];
		}
	}

	template <class M> std::complex<double> MatrixFieldApproximation::mean_value(const M& op, const Eigen::VectorXcd& initial_pure_state) const
	{
		return initial_pure_state.dot(op*initial_pure_state);
		//return initial_pure_state.dot(op.block(0, 0, m_nbr_sites, m_nbr_sites) * initial_pure_state);
		/*std::complex<double> sum(0.0);
		const unsigned int ratio = m_operator_dim / m_nbr_sites;
		for (size_t r = 0; r < ratio; ++r) {
			sum += initial_pure_state.dot(op.block(r*m_nbr_sites, r*m_nbr_sites, m_nbr_sites, m_nbr_sites) * initial_pure_state);
		}
		return sum;*/
	}

	void MatrixFieldApproximation::exciton_density_matrix(const data_type& y, const Eigen::VectorXcd& initial_pure_state, Eigen::MatrixXcd& rho)
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
		Eigen::MatrixXcd& tmp = m_work_g[0];
		const unsigned int ratio = m_operator_dim / m_nbr_sites;
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
					rho(m, n) += mean_value(tmp, initial_pure_state); // TODO: make it work for T > 0 or simulations with bath out of ground state at the beginning
				}
			}
		}
		rho /= m_nbr_sites;
	}

	void MatrixFieldApproximation::operator()(double t, const data_type& y, data_type& dy_dt)
	{
		assert(y.size() == m_dim);
		assert(dy_dt.size() == m_dim);
		dy_dt.setZero();
		
		if (m_fixed_bath) {
			for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
				bath_operator(dy_dt, k).setZero();
			}
		} else {
			// calculate d\tilde{b}_k(t)/dt
			for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
				const std::complex<double> expwk(exp(std::complex<double>(0.0, m_omega[k]*t)));
				m_exp_w[k] = expwk;
				operator_type sum(bath_operator(dy_dt, k));
				sum.setZero();
				const Eigen::MatrixXcd::RowXpr gk(m_g.row(k));
				for (size_t m = 0; m < m_nbr_sites; ++m) {
					sum.noalias() += exciton_operator(y, m, m) * (gk[m] * expwk);
				}
				sum *= -IMAGINARY;
			}
		}

		// set m_work_g[m] to sum_k ( \tilde{b}_k(t) e^{-i \omega_k t} \cc{g}_{km} + c.c. )
		// uses some memory from dy_dt as temporary storage
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			operator_type sum(dy_dt.data(), m_nbr_sites, m_nbr_sites);			
			sum.setZero();
			const Eigen::MatrixXcd::ColXpr gm(m_g.col(m));
			for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
				const std::complex<double>& gmk = gm[k];
				if (gmk != 0.0) { // TODO: check if it's more efficient to have this if, or not
					sum.noalias() += bath_operator(y, k).adjoint() * (gmk*m_exp_w[k]);
				}
			}
			m_work_g[m] = sum;
			m_work_g[m].noalias() += sum.adjoint();			
			assert(cfa::is_hermitian(m_work_g[m]));
		}

		for (size_t m = 0; m < m_nbr_sites; ++m) {
			const Eigen::MatrixXcd::ColXpr H0m(m_H0.col(m));
			for (size_t n = 0; n <= m; ++n) {
				operator_type sum(exciton_operator(dy_dt, m, n));
				sum.setZero();
				if (m != n) {					
					const_operator_type y_op(exciton_operator(y, m, n));
					// symmetrized product, because bath and exciton matrices do not commute in MFA approximation
					sum.noalias() = y_op * m_work_g[m];
					sum.noalias() -= y_op * m_work_g[n];
					sum.noalias() += m_work_g[m] * y_op;
					sum.noalias() -= m_work_g[n] * y_op;
					sum /= 2;
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

		operator_type dham(dy_dt.data() + m_hamiltonian_offset, m_operator_dim, m_operator_dim);
		dham.setZero();
		// calculate sum_{m,n} H0(m,n)|m><n| + sum_m |m><m| sum_k ( \tilde{b}_k(t) e^{-i \omega_k t} \cc{g}_{km} + c.c. )
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			const Eigen::MatrixXcd::RowXpr H0m(m_H0.row(m));
			const_operator_type op_mm(exciton_operator(y, m, m));
			assert(cfa::is_hermitian(op_mm));
			assert(H0m[m].imag() == 0);
			dham += op_mm * H0m[m];
			assert(cfa::is_hermitian(dham));
			assert(cfa::is_hermitian(m_work_g[m]));
			dham.noalias() += op_mm * m_work_g[m];
			cfa::make_hermitian(dham); // this is necessary because op_mm and m_work_g[m] do not have to commute in the MFA approximation
			for (size_t n = 0; n < m; ++n) {
				const_operator_type op(exciton_operator(y, m, n));
				dham += op * H0m[n];
				dham += op.adjoint() * conj(H0m[n]);
			}			
			assert(cfa::is_hermitian(dham));			
		}
		// previous contents of m_work_g are not needed anymore
		// we will reuse m_work_g[0]
		// calculate sum_k omega_k b^dagger_k b_k
		assert(!m_work_g.empty());
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			const_operator_type bk(bath_operator(y, k));
			m_work_g.front().noalias() = bk.adjoint() * bk;
			dham += m_work_g.front() * m_omega[k];
		}
		assert(cfa::is_hermitian(dham));
	}

	void MatrixFieldApproximation::reset_H0(const Eigen::MatrixXcd& H0)
	{
		m_H0 = H0;
	}

	double MatrixFieldApproximation::energy_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{
		const_operator_type H(hamiltonian_operator(y));
		return mean_value(H, initial_pure_state).real();
	}

	std::complex<double> MatrixFieldApproximation::evolution_operator_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{		
		Eigen::MatrixXcd minus_hamiltonian_operator(m_operator_dim, m_operator_dim);
		minus_hamiltonian_operator = hamiltonian_operator(y);
		minus_hamiltonian_operator *= -IMAGINARY;
		Eigen::MatrixExponential<Eigen::MatrixXcd> exponential(minus_hamiltonian_operator); // evolution operator
		Eigen::MatrixXcd U(m_operator_dim, m_operator_dim);
		exponential.compute(U);
		return mean_value(U, initial_pure_state);
	}

	std::complex<double> MatrixFieldApproximation::exciton_density_matrix_trace(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{
		std::complex<double> sum(0.0);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			sum += mean_value(exciton_operator(y, m, m), initial_pure_state);
		}
		return sum;
	}

	const double MatrixFieldApproximation::RHO_EIGENVALUE_TOLERANCE = 1E-3;
}

