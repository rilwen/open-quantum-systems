#ifndef __CFA_PROCESSOR_DENSITY_MATRIX_H
#define __CFA_PROCESSOR_DENSITY_MATRIX_H

#include <Eigen/Core>
#include <vector>

namespace cfa {
	template <class Method> class ProcessorDensityMatrix
	{
	public:
		typedef Method method_type;
		ProcessorDensityMatrix(const Eigen::VectorXcd& initial_state, const size_t nbr_steps);
		void operator()(Method& method, const typename Method::data_type& state, size_t step_idx, double t);
		const Eigen::MatrixXcd& rho(size_t step_idx) const { return m_rhos[step_idx]; }

		//! Theoretical value of the rho trace
		std::complex<double> rho_trace(size_t step_idx) const { return m_rho_traces[step_idx]; }
	private:
		Eigen::VectorXcd m_initial_state;
		std::vector<Eigen::MatrixXcd> m_rhos;
		std::vector<std::complex<double> > m_rho_traces;
	};

	template <class M> ProcessorDensityMatrix<M>::ProcessorDensityMatrix(const Eigen::VectorXcd& initial_state, const size_t nbr_steps)
		: m_initial_state(initial_state), m_rhos(nbr_steps + 1), m_rho_traces(nbr_steps + 1)
	{
	}

	template <class M> void ProcessorDensityMatrix<M>::operator()(M& method, const typename M::data_type& state, size_t step_idx, double)
	{
		method.exciton_density_matrix(state, m_initial_state, m_rhos[step_idx]);
		m_rho_traces[step_idx] = method.exciton_density_matrix_trace(state, m_initial_state);
	}
}

#endif // __CFA_PROCESSOR_DENSITY_MATRIX_H
