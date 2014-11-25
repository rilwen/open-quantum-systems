#include "partial_function_decomposition_bose.h"
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

static std::string build_pfd_filename(const unsigned int order)
{
	std::stringstream ss;
	ss << "partial_fraction_decomposition_result_" << order << ".dat";
	return ss.str();
}

static bool load_from_file(const unsigned int order, std::vector<std::complex<double> >& xi)
{	
	try {
		std::ifstream inf(build_pfd_filename(order).c_str(), std::ios::in);
		if (inf.good()) {
			for (unsigned int i = 0; i < order; ++i) {
				inf >> xi[i];
			}
		} else {
			return false;
		}
		return true;
	} catch (std::exception&) {
		return false;
	}
}

static void save_to_file(const unsigned int order, const std::vector<std::complex<double> >& xi)
{
	std::ofstream outf(build_pfd_filename(order).c_str(), std::ios::out);
	for (unsigned int i = 0; i < order; ++i) {
		outf << xi[i] << "\n";
	}
	outf.close();
}

PartialFunctionDecompositionBose::PartialFunctionDecompositionBose(const unsigned int order)
	: m_order(order), m_xi(order)
{
	if (!load_from_file(order, m_xi)) {
		Eigen::MatrixXcd A(order, order);
		A.setZero();
		if (order) {
			for (unsigned int m = 0; m < order - 1; ++m) {
				A(m, m + 1) = 2*(m+1)*(2*m + 3);
			}
			const double a = 2*order*(2*order + 1);
			for (unsigned int n = 0; n < order; ++n) {
				A(order - 1, n) -= a;
			}
			Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver;
			solver.compute(A);
			for (unsigned int m = 0; m < order; ++m) {
				m_xi[m] = A.eigenvalues()[m];
			}
		}
		try {
			save_to_file(order, m_xi);
		} catch (std::exception& e) {
			std::cerr << "Error saving PFD to file: " << e.what() << "\n";
		}
	}
}

std::complex<double> PartialFunctionDecompositionBose::evaluate(const std::complex<double>& z) const
{
	if (z != 0.0) {
		std::complex<double> result = 0.5 + 1.0 / z;
		for (unsigned int j = 0; j < m_order; ++j) {
			result += 2.0 * z / (z * z - 4.0 * m_xi[j]);
		}
		return result;
	} else {
		return std::numeric_limits<double>::infinity();
	}
} 

double PartialFunctionDecompositionBose::full_evaluate(double x) const
{
	return 1.0 / ( 1 - exp(-x) );
}

std::complex<double> PartialFunctionDecompositionBose::full_evaluate(const std::complex<double>& x) const
{
	return 1.0 / ( 1.0 - exp(-x) );
}

