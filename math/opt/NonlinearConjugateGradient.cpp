#include "NonlinearConjugateGradient.h"

namespace rql {
	NonlinearConjugateGradient::Function::Function(unsigned int dim)
		: m_dim(dim)
	{
	}

	NonlinearConjugateGradient::Workspace::Workspace(unsigned int dim, const Eigen::VectorXd& x0)
		: dimension(dim), nbrIterations(0), nbrFunctionEvaluations(0), nbrGradientEvaluations(0), x(x0), p(dim), g(dim)
	{		
	}

	NonlinearConjugateGradient::LineSearchFunction::LineSearchFunction(NonlinearConjugateGradient::Function& multivariate, NonlinearConjugateGradient::Workspace& wksp, const Eigen::VectorXd& x0, const Eigen::VectorXd& p)
		: m_fun(multivariate), m_wksp(wksp), m_x0(x0), m_p(p)
	{
	}

	void NonlinearConjugateGradient::optimize(NonlinearConjugateGradient::Function& f, const Eigen::VectorXd& x0, NonlinearConjugateGradient::Result& result) const
	{
		Workspace wksp(f.dimension(), x0);
		do_optimize(wksp, f, x0, result);
	}

	void NonlinearConjugateGradient::do_optimize(NonlinearConjugateGradient::Workspace& wksp, NonlinearConjugateGradient::Function& f, const Eigen::VectorXd& x0, NonlinearConjugateGradient::Result& result) const
	{
		call_vag(wksp, f);
		wksp.p = -wksp.g;
		Status status;
		while ((status = m_convergence_checker->status(wksp)) == CONTINUE)
		{
			wksp.nbrIterations++;
			LineSearchFunction lsf(f, wksp, wksp.x, wksp.p);
			wksp.prev_alpha = wksp.alpha;
			wksp.alpha = m_line_search->search(lsf, wksp);
			wksp.prev_x = wksp.x;
			wksp.prev_g = wksp.g;
			wksp.prev_f = wksp.f;
			wksp.x += wksp.alpha*wksp.p;			
			call_vag(wksp, f);
			m_update_strategy->update(wksp, wksp.p);
		}
		result.converged = status == CONVERGED;
		result.finalGradient = wksp.g;
		result.finalPosition = wksp.x;
		result.finalValue = wksp.f;
		result.nbrFunctionEvaluations = wksp.nbrFunctionEvaluations;
		result.nbrGradientEvaluations = wksp.nbrGradientEvaluations;
		result.nbrIterations = wksp.nbrIterations;
	}
}