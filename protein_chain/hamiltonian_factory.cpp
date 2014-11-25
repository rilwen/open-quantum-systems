#include "hamiltonian_factory.h"
#include "utils.h"
#include <stdexcept>

template <class M> static void templChain(const std::vector<double>& epsilons, double J, M& m)
{
	m.fill(0.0);
	//if (J > 0) throw std::domain_error("J cannot be positive");
	if (m.rows() != epsilons.size() || m.cols() != epsilons.size()) throw std::domain_error("Matrix size wrong");
	for (unsigned int i = 0; i < epsilons.size(); ++i) {
		m(i,i) = epsilons[i];
		if (i>0) {
			m(i-1,i) = m(i,i-1) = J;
		}
	}
}

void HamiltonianFactory::buildChain(const std::vector<double>& epsilons, double J, double kappa, Eigen::MatrixXcd& m)
{
	templChain(epsilons, J, m);
	setImag(m(epsilons.size() - 1, epsilons.size() - 1), kappa);
}

void HamiltonianFactory::buildChain(const std::vector<double>& epsilons, double J, Eigen::MatrixXd& m)
{
	templChain(epsilons, J, m);
}

template <class M> static void templDipoles(const std::vector<double>& epsilons, double J, M& m)
{
	m.fill(0.0);
	//if (J > 0) throw std::domain_error("J cannot be positive");
	if (m.rows() != epsilons.size() || m.cols() != epsilons.size()) throw std::domain_error("Matrix size wrong");
	for (unsigned int i = 0; i < epsilons.size(); ++i) {
		m(i,i) = epsilons[i];
		for (unsigned int j = 0; j < epsilons.size(); ++j) {
			if (j != i) {
				m(i, j) = J / pow(std::abs(static_cast<double>(i) - j), 3);
			}
		}
	}
}

void HamiltonianFactory::buildDipoles(const std::vector<double>& epsilons, double J, double kappa, Eigen::MatrixXcd& m)
{
	templDipoles(epsilons, J, m);
	setImag(m(epsilons.size() - 1, epsilons.size() - 1), kappa);
}

void HamiltonianFactory::buildDipoles(const std::vector<double>& epsilons, double J, Eigen::MatrixXd& m)
{
	templDipoles(epsilons, J, m);
}

template <class M> static void templPosDip(double epsilons, const std::vector<double>& xi, M& m)
{
	m.fill(0.0);

	for (unsigned int i = 0; i < xi.size(); ++i) {
		m(i,i) = epsilons;
		for (unsigned int j = 0; j < xi.size(); ++j) {
			if (j != i) {
				m(i, j) = -1 / pow(std::abs(static_cast<double>(i) - j + xi[i] - xi[j]), 3);
			}
		}
	}
}

void HamiltonianFactory::buildPosDip(double epsilons, const std::vector<double>& xi, double kappa, Eigen::MatrixXcd& m)
{
	templPosDip(epsilons, xi, m);
	setImag(m(xi.size() - 1, xi.size() - 1), kappa);
}

void HamiltonianFactory::buildPosDip(double epsilons, const std::vector<double>& xi, Eigen::MatrixXd& m)
{
	templPosDip(epsilons, xi, m);
}

template <class M> static void templWheel(const std::vector<double>& epsilons, double J, M& m)
{
        m.fill(0.0);

        for (unsigned int i = 0; i < epsilons.size(); ++i) {
                m(i,i) = epsilons[i];
                if (i>0) {
                        m(i-1,i) = m(i,i-1) = J;
                        m(i-1, epsilons.size()-1) = m(epsilons.size()-1, i-1) = J; // coupling with the hub
                }
        }
        m(epsilons.size()-2,0) = m(0,epsilons.size()-2) = J;

}

void HamiltonianFactory::buildWheel(const std::vector<double>& epsilons, double J, double kappa, Eigen::MatrixXcd& m)
{
        templWheel(epsilons, J, m);
		setImag(m(epsilons.size() - 1, epsilons.size() - 1), kappa);
}

void HamiltonianFactory::buildWheel(const std::vector<double>& epsilons, double J, Eigen::MatrixXd& m)
{
        templWheel(epsilons, J, m);
}

void HamiltonianFactory::buildH(const std::vector<double>& epsilons, double J, Eigen::MatrixXcd& m)
{
        m.fill(0.0);

		for (unsigned int i = 0; i < epsilons.size(); ++i) {
                m(i,i) = epsilons[i];
                if (i>0 && i<6) {
                        m(i-1,i) = m(i,i-1) = J;
                }
        }
		m(2,3) = m(3,2) = 0;
		m(6,1) = m(1,6) = m(6,4) = m(4,6) = J;
}

void HamiltonianFactory::buildHnobar(const std::vector<double>& epsilons, double J, Eigen::MatrixXcd& m)
{
        m.fill(0.0);

		for (unsigned int i = 0; i < epsilons.size(); ++i) {
                m(i,i) = epsilons[i];
                if (i>0 && i<6) {
                        m(i-1,i) = m(i,i-1) = J;
                }
        }
		m(2,3) = m(3,2) = 0;
		m(1,4) = m(4,1) = J;
}