#include "strunz_simulator_accumulator.h"

StrunzSimulatorAccumulator::StrunzSimulatorAccumulator(const Eigen::VectorXcd& initial, std::vector<double>& scalarProductTimes, std::vector<std::complex<double> >& scalarProducts)
	: m_initial(initial), m_scalarProductTimes(scalarProductTimes), m_scalarProducts(scalarProducts), m_base_time(0.0)
{
}

void StrunzSimulatorAccumulator::setBaseTime(double baseTime)
{
	m_base_time = baseTime;
}

void StrunzSimulatorAccumulator::operator()(double delta, const Eigen::VectorXcd& state)
{
	m_scalarProductTimes.push_back(m_base_time + delta);
	const std::complex<double> sp(m_initial.dot(state));
	m_scalarProducts.push_back(sp);
	//std::cout << m_base_time + delta << "\t" << std::abs(sp) << "\n";
}
