#ifndef __RQL_NONLINEAR_CONJUGATE_GRADIENT_H
#define __RQL_NONLINEAR_CONJUGATE_GRADIENT_H

#include <Eigen/Core>
#include <memory>

namespace rql {
	class NonlinearConjugateGradient
	{
	public:
		class Function
		{
		public:
			virtual ~Function() {}
			Function(unsigned int dim);
			unsigned int dimension() const { return m_dim; }
			virtual double value(const Eigen::VectorXd& x) = 0;
			virtual void gradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient) = 0;
			virtual double valueAndGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient) = 0;
		private:
			unsigned int m_dim;
		};
		struct Workspace
		{
			Workspace(unsigned int dim, const Eigen::VectorXd& x0);
			const unsigned int dimension;
			unsigned int nbrIterations;
			unsigned int nbrFunctionEvaluations;
			unsigned int nbrGradientEvaluations;
			Eigen::VectorXd x;
			Eigen::VectorXd p;
			Eigen::VectorXd g;
			double f;
			double alpha;
			Eigen::VectorXd prev_x;
			Eigen::VectorXd prev_g;
			Eigen::VectorXd prev_p;
			double prev_f;
			double prev_alpha;
		private:
			Workspace& operator=(const Workspace&); // not implemented
		};
		class LineSearchFunction
		{
		public:
			LineSearchFunction(Function& multivariate, Workspace& wksp, const Eigen::VectorXd& x0, const Eigen::VectorXd& p);
			double value(double alpha)
			{
				const double v = m_fun.value(m_x0 + m_p*alpha);
				m_wksp.nbrFunctionEvaluations++;
				return v;
			}
		private:
			LineSearchFunction& operator=(const LineSearchFunction&); // not implemented
			Function& m_fun;
			Workspace& m_wksp;
			const Eigen::VectorXd& m_x0;
			const Eigen::VectorXd& m_p;
		};
		class LineSearch
		{
		public:
			virtual ~LineSearch() {}
			virtual double search(LineSearchFunction& lsf, const Workspace& wksp) = 0;
		};
		class UpdateStrategy
		{
		public:
			virtual ~UpdateStrategy() {}
			virtual void update(const Workspace& wksp, Eigen::VectorXd& new_p) = 0;
		};
		class ResetStrategy
		{
		public:
			virtual ~ResetStrategy() {}
			virtual bool reset(const Workspace& wksp) = 0;
		};
		typedef enum { CONTINUE, CONVERGED, FAIL } Status;
		class ConvergenceChecker
		{
		public:
			virtual ~ConvergenceChecker() {}
			virtual Status status(const Workspace& wksp) = 0;
		};
		struct Result
		{
			bool converged;
			Eigen::VectorXd finalPosition;
			double finalValue;
			Eigen::VectorXd finalGradient;
			unsigned int nbrIterations;
			unsigned int nbrFunctionEvaluations;
			unsigned int nbrGradientEvaluations;
		};

		void optimize(Function& f, const Eigen::VectorXd& x0, Result& result) const;
	private:
		void do_optimize(Workspace& wksp, Function& f, const Eigen::VectorXd& x0, Result& result) const;
		void call_vag(Workspace& wksp, Function& fun) const
		{
			wksp.f = fun.valueAndGradient(wksp.x, wksp.g);
			wksp.nbrFunctionEvaluations++;
			wksp.nbrGradientEvaluations++;
		}		
	private:		
		std::shared_ptr<LineSearch> m_line_search;
		std::shared_ptr<UpdateStrategy> m_update_strategy;
		std::shared_ptr<ResetStrategy> m_reset_strategy;
		std::shared_ptr<ConvergenceChecker> m_convergence_checker;
	};
}

#endif // __RQL_NONLINEAR_CONJUGATE_GRADIENT_H
