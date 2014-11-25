#include "classical_field_approximation.h"
#include "cfa_utils.h"
#include <cassert>
#include <complex>

namespace cfa {
	ClassicalFieldApproximation::ClassicalFieldApproximation(const Eigen::MatrixXcd& H0, const Eigen::MatrixXcd& g, const Eigen::VectorXd& omega, double kBT)
		: m_H0(H0), m_g(g), m_omega(omega), m_exp_w(omega.size()), m_work_g(H0.rows()), m_kBT(kBT)
	{
		m_nbr_sites = H0.rows();
		assert(H0.cols() == H0.rows());
		m_nbr_bath_modes = omega.size();
		assert(g.rows() == m_nbr_bath_modes);
		assert(g.cols() == m_nbr_sites);
		m_nbr_exciton_values = (m_nbr_sites*(m_nbr_sites + 1))/2;
		m_dim = m_nbr_exciton_values + m_nbr_bath_modes;
	}

	void ClassicalFieldApproximation::init_exciton_part(const Eigen::VectorXcd& exciton, data_type& y) const
	{
		y.resize(m_dim);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			const std::complex<double> e_m = conj(exciton[m]);
			for (size_t n = 0; n <= m; ++n)
				y[exciton_index(m, n)] = e_m * exciton[n];
		}
	}

	void ClassicalFieldApproximation::init_from_pure_state(const Eigen::VectorXcd& exciton, data_type& y) const
	{
		init_exciton_part(exciton, y);
		assert(y.size() == m_dim);
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			y[bath_index(k)] = 0.0;
		}
	}

	void ClassicalFieldApproximation::init_from_pure_state(const Eigen::VectorXcd& exciton, const Eigen::VectorXcd& bath, data_type& y) const
	{
		init_exciton_part(exciton, y);		
		assert(y.size() == m_dim);
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			y[bath_index(k)] = bath[k];
		}
	}

	void ClassicalFieldApproximation::exciton_density_matrix(const data_type& y, const Eigen::VectorXcd&, Eigen::MatrixXcd& rho) const
	{
		rho.resize(m_nbr_sites, m_nbr_sites);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			for (size_t n = 0; n <= m; ++n) {
				rho(m, n) = y[exciton_index(m, n)];
				rho(n, m) = conj(rho(m, n));
			}
		}		
	}

	std::complex<double> ClassicalFieldApproximation::exciton_density_matrix_trace(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{
		std::complex<double> sum(0.0);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			sum += y[exciton_index(m, m)];
		}
		return sum;
	}

	void ClassicalFieldApproximation::operator()(double t, const data_type& y, data_type& dy_dt)
	{
		assert(y.size() == m_dim);
		assert(dy_dt.size() == m_dim);
		dy_dt.setZero();
		
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			const std::complex<double> expwk(exp(std::complex<double>(0.0, m_omega[k]*t)));
			m_exp_w[k] = expwk;
			std::complex<double> sum(0.0,0.0);
			const Eigen::MatrixXcd::RowXpr gk(m_g.row(k));
			for (size_t m = 0; m < m_nbr_sites; ++m) {
				sum += gk[m] * expwk * y[exciton_index(m, m)];
			}
			dy_dt[bath_index(k)] = -(IMAGINARY * sum);
		}
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			std::complex<double> sum(0.0, 0.0);
			const Eigen::MatrixXcd::ColXpr gm(m_g.col(m));
			for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
				sum += gm[k] * conj(y[bath_index(k)]) * m_exp_w[k];
			}
			m_work_g[m] = sum.real();
		}
		m_work_g *= IMAGINARY * 2.0;

		for (size_t m = 0; m < m_nbr_sites; ++m) {
			const Eigen::MatrixXcd::ColXpr H0m(m_H0.col(m));
			for (size_t n = 0; n <= m; ++n) {				
				std::complex<double> sum(0.0,0.0);
				const Eigen::MatrixXcd::RowXpr H0n(m_H0.row(n));
				size_t m2 = 0;
				for (; m2 < n; ++m2) {
					assert(m2 < n);
					assert(m2 < m);
					sum += H0m[m2]*conj(y[exciton_index(n,m2)]) - H0n[m2]*y[exciton_index(m,m2)];
				}
				for (; m2 < m; ++m2) {
					assert(m2 >= n);
					assert(m2 < m);
					sum += H0m[m2]*y[exciton_index(m2,n)] - H0n[m2]*y[exciton_index(m,m2)];
				}
				for (; m2 < m_nbr_sites; ++m2) {
					assert(m2 >= n);
					assert(m2 >= m);
					sum += H0m[m2]*y[exciton_index(m2,n)] - H0n[m2]*conj(y[exciton_index(m2,m)]);
				}
				sum *= IMAGINARY;
				const size_t idx = exciton_index(m, n);
				if (m != n) {					
					sum += y[idx] * (m_work_g[m] - m_work_g[n]);
				}
				dy_dt[idx] = sum;
			}
		}
	}

	void ClassicalFieldApproximation::reset_H0(const Eigen::MatrixXcd& H0)
	{
		m_H0 = H0;
	}

	double ClassicalFieldApproximation::energy_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{
		return 0;
	}

	const double ClassicalFieldApproximation::RHO_EIGENVALUE_TOLERANCE = 1E-3;
}
