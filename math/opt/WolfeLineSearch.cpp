#include "WolfeLineSearch.h"

#include <stdexcept>

namespace rql {
	WolfeLineSearch::WolfeLineSearch(double c1, double c2, unsigned int maxIter)
		: m_c1(c1), m_c2(c2), m_max_iter(maxIter)
	{
		if (m_c1 <= 0.0)
			throw std::domain_error("Required: c1 > 0");
		if (m_c2 <= m_c1)
			throw std::domain_error("Required: c1 < c2");
		if (m_c2 >= 1.0)
			throw std::domain_error("Required: c2 < 1");
		if (!m_max_iter)
			throw std::domain_error("Required: maxIter > 0");
	}

	double WolfeLineSearch::search(NonlinearConjugateGradient::LineSearchFunction& /*lsf*/, const NonlinearConjugateGradient::Workspace& /*wksp*/)
	{
		return 0;
	}

	const double WolfeLineSearch::DFLT_C1 = 1E-4;
	const double WolfeLineSearch::DFLT_C2 = 0.9;
}