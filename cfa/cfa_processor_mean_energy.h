#ifndef __CFA_PROCESSOR_MEAN_ENERGY_H
#define __CFA_PROCESSOR_MEAN_ENERGY_H

#include <Eigen/Core>
#include <vector>

namespace cfa {
	//! Measures int_0^t <Psi(0)|H(t)|Psi(0)> dt
	template <class Method> class ProcessorMeanEnergy
	{
	public:
		typedef Method method_type;
		ProcessorMeanEnergy(const Eigen::VectorXcd& initial_state, const size_t nbr_steps);
		void operator()(Method& method, const typename Method::data_type& state, size_t step_idx, double t);
		double energy(size_t step_idx) const { return m_energies[step_idx]; }
	private:
		Eigen::VectorXcd m_initial_state;
		std::vector<double> m_energies;
	};

	template <class M> ProcessorMeanEnergy<M>::ProcessorMeanEnergy(const Eigen::VectorXcd& initial_state, const size_t nbr_steps)
		: m_initial_state(initial_state), m_energies(nbr_steps + 1)
	{
	}

	template <class M> void ProcessorMeanEnergy<M>::operator()(M& method, const typename M::data_type& state, size_t step_idx, double)
	{
		m_energies[step_idx] = method.energy_mean(state, m_initial_state);
	}
}

#endif // __CFA_PROCESSOR_MEAN_ENERGY_H
