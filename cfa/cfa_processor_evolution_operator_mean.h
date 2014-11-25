#ifndef __CFA_PROCESSOR_EVOL_OP_MEAN_H
#define __CFA_PROCESSOR_EVOL_OP_MEAN_H

#include <Eigen/Core>
#include <complex>
#include <vector>
#include "classical_field_approximation.h"

namespace cfa {
	template <class Method> class ProcessorEvolutionOperatorMean
	{
	public:
		typedef Method method_type;
		ProcessorEvolutionOperatorMean(const Eigen::VectorXcd& initial_state, const size_t nbr_steps);
		void operator()(const Method& method, const typename Method::data_type& state, size_t step_idx, double t);
		const std::vector<std::complex<double>>& mean_values() const { return m_mvals; }
	private:
		Eigen::VectorXcd m_initial_state;
		std::vector<std::complex<double>> m_mvals;
	};

	template <class M> ProcessorEvolutionOperatorMean<M>::ProcessorEvolutionOperatorMean(const Eigen::VectorXcd& initial_state, const size_t nbr_steps)
		: m_initial_state(initial_state), m_mvals(nbr_steps + 1)
	{
	}

	template <class M> struct caller
	{
		static std::complex<double> call(const M& method, const typename M::data_type& state, const Eigen::VectorXcd& initial_state)
		{
			return method.evolution_operator_mean(state, initial_state);
		}
	};

	class ClassicalFieldApproximation;
	template <> struct caller<ClassicalFieldApproximation>
	{
		static std::complex<double> call(const ClassicalFieldApproximation&, const ClassicalFieldApproximation::data_type&, const Eigen::VectorXcd&);
	};

	template <class M> void ProcessorEvolutionOperatorMean<M>::operator()(const M& method, const typename M::data_type& state, size_t step_idx, double)
	{
		m_mvals[step_idx] = caller<M>::call(method, state, m_initial_state);
	}
}

#endif // __CFA_PROCESSOR_EVOL_OP_MEAN_H
