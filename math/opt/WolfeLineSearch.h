#include "NonlinearConjugateGradient.h"

namespace rql {
	class WolfeLineSearch: public NonlinearConjugateGradient::LineSearch
	{
	public:
		//! 0 < c1 < c2 < 1
		//! maxIter > 0
		WolfeLineSearch(double c1, double c2, unsigned int maxIter);
		virtual double search(NonlinearConjugateGradient::LineSearchFunction& lsf, const NonlinearConjugateGradient::Workspace& wksp);
		static const double DFLT_C1;
		static const double DFLT_C2;
	private:
		double zoom(NonlinearConjugateGradient::LineSearchFunction& lsf, const NonlinearConjugateGradient::Workspace& wksp, double alph_lo, double alph_hi);
		WolfeLineSearch& operator=(const WolfeLineSearch&); // not implemented
	private:
		const double m_c1;
		const double m_c2;
		const unsigned int m_max_iter;		
	};
}