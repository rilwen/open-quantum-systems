#include <cmath>
#include <gtest/gtest.h>
#include <boost/make_shared.hpp>
#include <protein_chain/evolver_sparse.h>
#include <protein_chain/evolver_taylor_expansion.h>

using namespace Eigen;

class EvolverSparseTest: public testing::Test
{
};

TEST_F(EvolverSparseTest, EigenSparse)
{
	const size_t first_index = 20;
	const size_t second_index = 120;
	const size_t dim = 200;
	SparseMatrix<std::complex<double> > m(dim,dim);
	m.coeffRef(first_index, first_index) = 1.0;
	m.coeffRef(second_index, second_index) = 0.5;
	m.coeffRef(first_index, second_index) = std::complex<double>(0, -0.2);
	m.coeffRef(second_index, first_index) = std::complex<double>(0, 0.2);
	MatrixXcd dense_matrix(2, 2);
	dense_matrix(0, 0) = m.coeff(first_index, first_index);
	dense_matrix(1, 1) = m.coeff(second_index, second_index);
	dense_matrix(0, 1) = m.coeff(first_index, second_index);
	dense_matrix(1, 0) = m.coeff(second_index, first_index);
	m *= std::complex<double>(0,-1);
	EigenMultiplicator multiplicator(boost::make_shared<EigenMultiplicator::matrix_type>(m));
	EvolverSparse<EigenMultiplicator> sparse_evolver(0.0001, multiplicator, 2);
	
	EvolverTaylorExpansion<4> taylor_evolver(dense_matrix, 0.001);
	Eigen::VectorXcd initialDense(2);
	initialDense[0] = std::complex<double>(0.34, -0.71);
	initialDense[1] = 1;
	Eigen::VectorXcd initialSparse(dim);
	initialSparse.setZero();
	initialSparse[first_index] = initialDense[0];
	initialSparse[second_index] = initialDense[1];
	std::vector<double> times(3);
	times[0] = 1;
	times[1] = 2;
	times[2] = 3;
	std::vector<Eigen::VectorXcd> statesSparse(3);
	std::vector<Eigen::VectorXcd> statesTaylor(3);
	sparse_evolver.evolve(initialSparse, times, statesSparse);
	taylor_evolver.evolve(initialDense, times, statesTaylor);
	const double tol = 2E-4;
	for (size_t i = 0; i < 3; ++i) {
		EXPECT_NEAR(0.0, std::abs(statesTaylor[i][0] - statesSparse[i][first_index]), tol);
		EXPECT_NEAR(0.0, std::abs(statesTaylor[i][1] - statesSparse[i][second_index]), tol);
		for (size_t j = 0; j < dim; ++j)
			if (j != first_index && j != second_index)
				EXPECT_EQ(0.0, std::abs(statesSparse[i][j]));
	}
}
