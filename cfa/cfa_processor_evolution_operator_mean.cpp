#include "cfa_processor_evolution_operator_mean.h"

namespace cfa {
	std::complex<double> caller<ClassicalFieldApproximation>::call(const ClassicalFieldApproximation&, const ClassicalFieldApproximation::data_type&, const Eigen::VectorXcd&)
	{
		return std::complex<double>(1.0, 0.0);
	}
}
