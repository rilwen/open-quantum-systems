#include <gtest/gtest.h>
#include "exact_hamiltonian.h"

using namespace Eigen;

class ExactHamiltonianTest: public testing::Test
{
};

TEST_F(ExactHamiltonianTest,NearestNeighbourSingleMode)
{
	const size_t nbr_sites = 2;
	const size_t nbr_excited_levels = 1;
	const double J = -0.5;
	const std::complex<double> g1(0,0.25);
	const std::complex<double> g2(0.25);
	VectorXcd omega(1);
	omega[0] = 1;
	MatrixXcd g(2, 1);
	g(0,0) = g1;
	g(1,0) = g2;
	NearestNeighbour nn(nbr_sites, omega, g, J);
	MatrixXcd ham;
	std::vector<unsigned int> bath_basis;
	nn.generateBathBasis(1, 1+ nbr_excited_levels, bath_basis);
	ASSERT_EQ(1+nbr_excited_levels, bath_basis.size());
	ASSERT_EQ(0u, bath_basis[0]);
	ASSERT_EQ(1u, bath_basis[1]);
	nn.generateMatrix(ham, nbr_excited_levels);
	ASSERT_TRUE(ham == ham.adjoint());
	ASSERT_EQ(J, ham(0, 2));
	ASSERT_EQ(J, ham(1, 3));
	ASSERT_EQ(0.0, ham(0, 0));
	ASSERT_EQ(omega[0], ham(1, 1));
	ASSERT_EQ(0.0, ham(2, 2));
	ASSERT_EQ(omega[0], ham(3, 3));
	ASSERT_EQ(g1, ham(1,0));
	ASSERT_EQ(conj(g1), ham(0,1));
	ASSERT_EQ(g2, ham(3,2));
	ASSERT_EQ(conj(g2), ham(2,3));
}

TEST_F(ExactHamiltonianTest,NearestNeighbourNoBath)
{
	const size_t nbr_sites = 2;
	const size_t nbr_excited_levels = 1;
	const double J = -0.5;
	const double g1 = 0.25;
	const double g2 = -0.25;
	VectorXcd omega(1);
	omega[0] = 1;
	MatrixXcd g(2, 1);
	g(0,0) = g1;
	g(1,0) = g2;
	NearestNeighbour nn(nbr_sites, omega, g, J);
	MatrixXcd actual;
	nn.generateMatrix(actual, 0);
	ASSERT_TRUE(actual == actual.adjoint());
	ASSERT_EQ(2, actual.rows());
	ASSERT_EQ(2, actual.cols());
	ASSERT_EQ(0.0, actual(0,0));
	ASSERT_EQ(0.0, actual(1,1));
	ASSERT_EQ(std::complex<double>(J,0), actual(0,1));
	ASSERT_EQ(std::complex<double>(J,0), actual(1,0));
}

TEST_F(ExactHamiltonianTest,NearestNeighbourTwoModes)
{
	const size_t nbr_sites = 2;
	const size_t nbr_modes = 2;
	const size_t nbr_excited_levels = 1;
	const size_t nbr_bath_levels = nbr_excited_levels + 1;
	const double J = -1;
	VectorXcd omega(nbr_modes);
	omega[0] = 0.75;
	omega[1] = 1.5;
	MatrixXcd g(nbr_sites, nbr_modes);
	g(0,0) = std::complex<double>(0, 0.25);
	g(0,1) = std::complex<double>(0, -0.25);
	g(1,0) = 0.25;
	g(1,1) = -0.25;
	NearestNeighbour nn(nbr_sites, omega, g, J);
	MatrixXcd ham;
	std::vector<unsigned int> bath_basis;
	nn.generateBathBasis(nbr_modes, nbr_bath_levels, bath_basis);
	ASSERT_EQ(nbr_modes * nbr_bath_levels * nbr_bath_levels, bath_basis.size());
	ASSERT_EQ(0u, bath_basis[0]);
	ASSERT_EQ(0u, bath_basis[1]);
	ASSERT_EQ(1u, bath_basis[2]);
	ASSERT_EQ(0u, bath_basis[3]);
	ASSERT_EQ(0u, bath_basis[4]);
	ASSERT_EQ(1u, bath_basis[5]);
	ASSERT_EQ(1u, bath_basis[6]);
	ASSERT_EQ(1u, bath_basis[7]);
	nn.generateMatrix(ham, nbr_excited_levels);
	ASSERT_TRUE(ham == ham.adjoint());
	ASSERT_EQ(8, ham.rows());
	ASSERT_EQ(8, ham.cols());
	for (size_t i1 = 0; i1 < 8; ++i1) {
		const size_t m = i1 / 4; // system basis index
		const size_t b1 = i1 % 4; // bath basis index
		for (size_t i2 = 0; i2 < 8; ++i2) {
			const size_t n = i2 / 4; // system basis index
			const size_t b2 = i2 % 4; // bath basis index
			if (m == n) {
				if (b1 == b2) { // h omega[k] a_k^dagger a_k
					std::complex<double> sum = 0;
					for (size_t k = 0; k < nbr_modes; ++k)
						sum += omega[k] * static_cast<double>(bath_basis[b1 * nbr_modes + k]);
					ASSERT_EQ(sum.real(), ham(i1, i2));
					ASSERT_EQ(sum.imag(), 0.0);
				} else { // bath-system coupling
					const unsigned int o1_0 = bath_basis[b1 * nbr_modes];
					const unsigned int o1_1 = bath_basis[b1 * nbr_modes + 1];
					const unsigned int o2_0 = bath_basis[b2 * nbr_modes];
					const unsigned int o2_1 = bath_basis[b2 * nbr_modes + 1];
					if (o1_0 == o2_0) {
						if (o1_1 < o2_1) // annihilation of mode 1
							ASSERT_EQ(conj(g(m, 1)), ham(i1, i2));
						else { // creation of mode 1
							assert(o1_1 > o2_1);
							ASSERT_EQ(g(m, 1), ham(i1, i2));
						}
					} else if (o1_1 == o2_1) {
						if (o1_0 < o2_0) // annihilation of mode 0
							ASSERT_EQ(conj(g(m, 0)), ham(i1, i2));
						else { // creation of mode 0
							assert(o1_0 > o2_0);
							ASSERT_EQ(g(m, 0), ham(i1, i2));
						}
					} else {
						ASSERT_EQ(0.0, ham(i1,i2));
					}
				}
			} else {
				if (b1 == b2) {
					if (abs(static_cast<int>(m - n)) == 1)
						ASSERT_EQ(J, ham(i1, i2)); // V_mn |m><n|
				} else {
					ASSERT_EQ(0.0, ham(i1, i2));
				}
			}
		}
	}
}

TEST_F(ExactHamiltonianTest,Sparse)
{
	const size_t nbr_sites = 2;
	const size_t nbr_modes = 2;
	const size_t nbr_excited_levels = 1;
	const size_t nbr_bath_levels = nbr_excited_levels + 1;
	const double J = -1;
	VectorXcd omega(nbr_modes);
	omega[0] = 0.75;
	omega[1] = 1.5;
	MatrixXcd g(nbr_sites, nbr_modes);
	g(0,0) = std::complex<double>(0, 0.25);
	g(0,1) = std::complex<double>(0, -0.25);
	g(1,0) = 0.25;
	g(1,1) = -0.25;
	NearestNeighbour nn(nbr_sites, omega, g, J);
	MatrixXcd dense;
	nn.generateMatrix(dense, nbr_excited_levels);
	const size_t dim = nn.dimension(nbr_excited_levels);
	SparseMatrix<std::complex<double> > sparse(dim, dim);
	nn.generateMatrix(sparse, nbr_excited_levels);
	ASSERT_EQ(dense.rows(), sparse.rows());
	ASSERT_EQ(dense.cols(), sparse.cols());
	for (int r = 0; r < dense.rows(); ++r)
		for (int c = 0; c < dense.cols(); ++c)
			ASSERT_EQ(dense(r, c), sparse.coeff(r, c));
}
