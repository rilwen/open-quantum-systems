#ifndef __CFA_PROCESSOR_PAIR_H
#define __CFA_PROCESSOR_PAIR_H

namespace cfa {
	template <class FirstProcessor, class SecondProcessor> class ProcessorPair
	{
	public:
		typedef typename FirstProcessor::method_type method_type;
		ProcessorPair(FirstProcessor& first, SecondProcessor& second)
			: m_first(first), m_second(second)
		{
		}
		void operator()(method_type& method, const typename method_type::data_type& state, size_t step_idx, double t)
		{
			m_first(method, state, step_idx, t);
			m_second(method, state, step_idx, t);
		}
		const FirstProcessor& first() const { return m_first; }
		const SecondProcessor& second() const { return m_second; }		
	private:
		FirstProcessor& m_first;
		SecondProcessor& m_second;
	};	

	template <class FP, class SP> ProcessorPair<FP, SP> join_processors(FP& first, SP& second)
	{
		return ProcessorPair<FP,SP>(first, second);
	}
}

#endif // __CFA_PROCESSOR_PAIR_H
