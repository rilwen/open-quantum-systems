#include "matrix_field_approximation_reduced.h"
#include "cfa_utils.h"
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <math/linalg/MatrixFunctions.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/special_functions/fpclassify.hpp>

namespace cfa {
	MatrixFieldApproximationReduced::MatrixFieldApproximationReduced(const Eigen::MatrixXcd& H0, const Eigen::MatrixXcd& g, const Eigen::VectorXd& omega, double kBT, bool normalize, size_t pivot)
		: m_H0(H0), m_g(g), m_omega(omega), m_exp_w(omega.size()), m_work_g(H0.rows()), m_pivot(pivot), m_qr(H0.rows(), H0.cols()), m_normalize(normalize), m_kBT(kBT)
	{		
		m_nbr_sites = H0.rows();
		m_operator_dim = m_nbr_sites;
		if (!m_nbr_sites) {
			throw std::domain_error("No sites");
		}
		if (pivot >= m_nbr_sites) {
			throw std::domain_error("");
		}
		assert(H0.cols() == H0.rows());
		m_nbr_bath_modes = omega.size();
		assert(g.rows() == m_nbr_bath_modes);
		assert(g.cols() == m_nbr_sites);
		m_nbr_operator_elements = m_operator_dim*m_operator_dim;
		m_hamiltonian_offset = (m_nbr_sites + m_nbr_bath_modes) * m_nbr_operator_elements;
		m_dim = m_hamiltonian_offset + m_nbr_operator_elements;
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			m_work_g[m].setZero(m_operator_dim, m_operator_dim);
		}
		for (size_t i = 0; i < NBR_TMP_MATRICES; ++i) {
			m_tmp_mat[i].setZero(m_operator_dim, m_operator_dim);
		}
		m_tmp_mat_real.setZero(m_nbr_sites, m_nbr_sites);
		for (size_t i = 0; i < NBR_TMP_VECTORS; ++i) {
			m_tmp_vec[i].setZero(m_nbr_sites);
		}
		m_sqrt_nbr_sites = sqrt(static_cast<double>(m_nbr_sites));
	}

	void MatrixFieldApproximationReduced::operator()(double t, const data_type& y, data_type& dy_dt)
	{
		assert(y.size() == m_dim);
		assert(dy_dt.size() == m_dim);
		dy_dt.setZero();
		
		double trace = 0;
		// calculate d\tilde{b}_k(t)/dt and Tr(sum_m |m><m|)
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			m_exp_w[k] = exp(std::complex<double>(0.0, m_omega[k]*t));
			bath_operator(dy_dt, k).setZero();
		}
		// by looping over m first and then over k, we calculate exciton_operator(y, m, m, 0) only once, minimizing the number of matrix*matrix multiplications
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			const Eigen::MatrixXcd& opmm = exciton_operator(y, m, m, 0);
			trace += opmm.trace().real();
			const Eigen::MatrixXcd::ColXpr gm(m_g.col(m));
			for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
				bath_operator(dy_dt, k).noalias() += exciton_operator(y, m, m, 0) * (gm[k] * m_exp_w[k]);
			}
		}
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			bath_operator(dy_dt, k) *= -IMAGINARY;
		}
		assert(trace >= 0.0);
		/*for (size_t k = 0; k < m_nbr_bath_modes; ++k) { // old version (does not calculate the trace)
			const std::complex<double> expwk(exp(std::complex<double>(0.0, m_omega[k]*t)));
			m_exp_w[k] = expwk;
			operator_type sum(bath_operator(dy_dt, k));
			sum.setZero();
			const Eigen::MatrixXcd::RowXpr gk(m_g.row(k));
			for (size_t m = 0; m < m_nbr_sites; ++m) {
				sum.noalias() += exciton_operator(y, m, m, 0) * (gk[m] * expwk);
			}
			sum *= -IMAGINARY;
		}*/

		// set m_work_g[m] to sum_k ( \tilde{b}_k(t) e^{-i \omega_k t} \cc{g}_{km} + c.c. )
		// uses some memory from dy_dt as temporary storage
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			operator_type sum(dy_dt.data(), m_operator_dim, m_operator_dim);
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
			const Eigen::MatrixXcd::RowXpr H0pivot(m_H0.row(m_pivot));
			operator_type sum(reduced_exciton_operator(dy_dt, m));
			sum.setZero();
			if (m != m_pivot) {
				const Eigen::MatrixXcd& y_op = exciton_operator(y, m, m_pivot, 0);
				//// symmetrized product, because bath and exciton matrices do not commute in MFA approximation
				sum.noalias() = y_op * m_work_g[m];
				sum.noalias() -= y_op * m_work_g[m_pivot];
				sum.noalias() += m_work_g[m] * y_op;
				sum.noalias() -= m_work_g[m_pivot] * y_op;
				sum /= 2;
			}				
				
			for (size_t m2 = 0; m2 < m_nbr_sites; ++m2) {
				sum += exciton_operator(y, m2, m_pivot, 0) * H0m[m2];
				sum -= exciton_operator(y, m, m2, 0) * H0pivot[m2];
			}
			sum *= IMAGINARY;
			assert(m != m_pivot || cfa::is_hermitian(sum));
		}

		if (m_normalize) {
			double trace_y_dy_dt = 0;
			for (size_t m = 0; m < m_nbr_sites; ++m) {
				const const_operator_type om(reduced_exciton_operator(y, m));
				m_tmp_mat[0].noalias() = om.adjoint() * reduced_exciton_operator(dy_dt, m);
				trace_y_dy_dt += m_tmp_mat[0].trace().real();
			}
			//assert(trace_y_dy_dt >= 0.0);???
			const double trace_ratio = trace_y_dy_dt / trace;
			const double sqrt_trace = sqrt(trace);
			for (size_t m = 0; m < m_nbr_sites; ++m) {
				operator_type om(reduced_exciton_operator(dy_dt, m));
				om -= reduced_exciton_operator(y, m) * trace_ratio;
				om *= m_sqrt_nbr_sites / sqrt_trace;
			}
		}

		operator_type dham(hamiltonian_operator(dy_dt));
		dham.setZero();
		// calculate sum_{m,n} H0(m,n)|m><n| + sum_m |m><m| sum_k ( \tilde{b}_k(t) e^{-i \omega_k t} \cc{g}_{km} + c.c. )
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			const Eigen::MatrixXcd::RowXpr H0m(m_H0.row(m));
			const Eigen::MatrixXcd& op_mm(exciton_operator(y, m, m, 0));
			assert(cfa::is_hermitian(op_mm));
			assert(H0m[m].imag() == 0);
			dham += op_mm * H0m[m];
			assert(cfa::is_hermitian(dham));
			assert(cfa::is_hermitian(m_work_g[m]));
			dham.noalias() += op_mm * m_work_g[m];
			cfa::make_hermitian(dham); // this is necessary because op_mm and m_work_g[m] do not have to commute in the MFA approximation
			for (size_t n = 0; n < m; ++n) {
				const Eigen::MatrixXcd& op(exciton_operator(y, m, n, 1));
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
		assert(dham.norm() != 0);
		assert(boost::math::isfinite(dham.norm()));
	}

	void MatrixFieldApproximationReduced::init_from_pure_state(const Eigen::VectorXcd& exciton, data_type& y) const
	{
		init_exciton_part(y);
		assert(y.size() == m_dim);
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			operator_type op(bath_operator(y, k));
			op.setZero();
			//op.topRightCorner(m_nbr_sites, m_nbr_sites).setIdentity();
		}
	}

	/*void MatrixFieldApproximationReduced::init_from_pure_state(const Eigen::VectorXcd& exciton, const Eigen::VectorXcd& bath, data_type& y) const
	{
		init_exciton_part(y);
		assert(y.size() == m_dim);
		for (size_t k = 0; k < m_nbr_bath_modes; ++k) {
			operator_type op(bath_operator(y, k));
			op.setIdentity();			
			op *= bath[k];
			op.topRightCorner(m_nbr_sites, m_nbr_sites).setIdentity();
		}
	}*/

	void MatrixFieldApproximationReduced::exciton_density_matrix(const data_type& y, const Eigen::VectorXcd& initial_pure_state, Eigen::MatrixXcd& rho)
	{
		rho.resize(m_nbr_sites, m_nbr_sites);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			for (size_t n = 0; n <= m; ++n) {
				const Eigen::MatrixXcd& op(exciton_operator(y, m, n, 0));
				rho(m, n) = initial_pure_state.dot(op.topLeftCorner(m_nbr_sites, m_nbr_sites)*initial_pure_state);
				rho(n, m) = conj(rho(m, n));
			}
		}		
	}

	std::complex<double> MatrixFieldApproximationReduced::exciton_density_matrix_trace(const data_type& y, const Eigen::VectorXcd& initial_pure_state)
	{
		std::complex<double> sum(0.0);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			const Eigen::MatrixXcd& op(exciton_operator(y, m, m, 0));
			sum += initial_pure_state.dot(op.topLeftCorner(m_nbr_sites, m_nbr_sites)*initial_pure_state);
		}
		return sum;
	}

	double MatrixFieldApproximationReduced::energy_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{
		const_operator_type H(hamiltonian_operator(y));
		return initial_pure_state.dot(H*initial_pure_state).real();
	}

	std::complex<double> MatrixFieldApproximationReduced::evolution_operator_mean(const data_type& y, const Eigen::VectorXcd& initial_pure_state) const
	{
		Eigen::MatrixXcd minus_hamiltonian_operator(m_nbr_sites, m_nbr_sites);
		minus_hamiltonian_operator = hamiltonian_operator(y);
		minus_hamiltonian_operator *= -IMAGINARY;
		assert(boost::math::isfinite(minus_hamiltonian_operator.norm()));
		Eigen::MatrixXcd U;
		rql::math::linalg::MatrixFunctions::hermitian_matrix_function(minus_hamiltonian_operator, static_cast<double(*)(double)>(exp), U);
		assert(U.norm() != 0);
		assert(boost::math::isfinite(U.norm()));
		const std::complex<double> sp = initial_pure_state.dot(U*initial_pure_state);
		assert(boost::math::isfinite(std::abs(sp)));
		if (!boost::math::isfinite(sp.real()) || !boost::math::isfinite(sp.imag())) {
			const std::complex<double> H_norm = minus_hamiltonian_operator.norm();
			const std::complex<double> U_norm = U.norm();
			const double abs_sp = std::abs(sp);
			std::stringstream ss;
			ss << "Numerical errors: " << H_norm << " " << U_norm << " " << abs_sp;
			throw std::runtime_error(ss.str().c_str());
		}
		return sp;
	}

	void MatrixFieldApproximationReduced::init_exciton_part(data_type& y) const
	{
		y.resize(m_dim);
		for (size_t m = 0; m < m_nbr_sites; ++m) {
			operator_type op(reduced_exciton_operator(y, m));
			op.setZero();
			op(m, m_pivot) = 1.0;
		}
		hamiltonian_operator(y).setZero();
	}	

	const Eigen::MatrixXcd& MatrixFieldApproximationReduced::exciton_operator(const data_type& y, size_t m, size_t n, size_t tmpidx)
	{
		assert(tmpidx < NBR_TMP_MATRICES);
		Eigen::MatrixXcd& op = m_tmp_mat[tmpidx];
		op.noalias() = reduced_exciton_operator(y, m) * reduced_exciton_operator(y, n).adjoint();
		return op;
	}

	MatrixFieldApproximationReduced::IdentityNormalizer::IdentityNormalizer(MatrixFieldApproximationReduced& ode)
		: m_ode(ode)
	{
	}

	void MatrixFieldApproximationReduced::IdentityNormalizer::normalize(data_type& state)
	{
		if (m_ode.m_normalize) { // commented out -- not a good normalisation method
	//		Eigen::MatrixXd& A = m_ode.m_tmp_mat_real;
	//		Eigen::VectorXd& b = m_ode.m_tmp_vec[0];
	//		Eigen::VectorXd& x = m_ode.m_tmp_vec[1];
	//		b.setConstant(1.0);
	//		for (size_t m = 0; m < m_ode.m_nbr_sites; ++m) {
	//			const Eigen::MatrixXcd& op_mm = m_ode.exciton_operator(state, m, m, 0);
	//			Eigen::MatrixXd::ColXpr Acol(A.col(m));
	//			for (size_t k = 0; k < m_ode.m_nbr_sites; ++k) {
	//				Acol[k] = op_mm(k, k).real();// + op_mm(k + m_ode.m_nbr_sites, k + m_ode.m_nbr_sites).real(); // A(k, m) = op_mm(k, k);
	//			}			
	//		}
	//		m_ode.m_qr.compute(A);
	//		x = m_ode.m_qr.solve(b);
	//#ifndef NDEBUG
	//		const double residual = (b - A*x).norm();
	//		assert(residual < 1E-8*A.norm());
	//#endif // NDEBUG
	//		for (size_t m = 0; m < m_ode.m_nbr_sites; ++m) {
	//			m_ode.reduced_exciton_operator(state, m) *= sqrt(std::abs(x[m]));
	//		}
		}
	}

	const double MatrixFieldApproximationReduced::RHO_EIGENVALUE_TOLERANCE = 1E-8;
}
