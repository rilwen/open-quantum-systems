#include <gtest/gtest.h>
#include <protein_chain/reductor.h>

class ReductorTest: public testing::Test
{
};

TEST_F(ReductorTest,Trace)
{
	const size_t dimA = 10;
	const size_t dimB = 5;
	const size_t dim = dimA * dimB;
	Eigen::VectorXcd state(dim);	
	state.setRandom();
	state.normalize();
	Eigen::MatrixXcd rho(dimA, dimA);
	Reductor reductor(dimA, dimB);
	ASSERT_EQ(dimA, reductor.dimA());
	ASSERT_EQ(dimB, reductor.dimB());
	reductor.reduce(state, rho);
	const std::complex<double> tr = rho.trace();
	ASSERT_NEAR(tr.real(), 1.0, 1E-10);
	ASSERT_NEAR(tr.imag(), 0.0, 1E-10);
}
