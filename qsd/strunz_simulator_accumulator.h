#ifndef __STRUNZ_SIMULATOR_ACCUMULATOR_H
#define __STRUNZ_SIMULATOR_ACCUMULATOR_H

#include <vector>
#include <Eigen/Core>

class StrunzSimulatorAccumulator
{
public:
	StrunzSimulatorAccumulator(const Eigen::VectorXcd& initial, std::vector<double>& scalarProductTimes, std::vector<std::complex<double> >& scalarProducts);
	void setBaseTime(double baseTime);
	void operator()(double delta, const Eigen::VectorXcd& state);
private:
	const Eigen::VectorXcd& m_initial;
	std::vector<double>& m_scalarProductTimes;
	std::vector<std::complex<double> >& m_scalarProducts;
	double m_base_time;
};

#endif // __STRUNZ_SIMULATOR_ACCUMULATOR_H
