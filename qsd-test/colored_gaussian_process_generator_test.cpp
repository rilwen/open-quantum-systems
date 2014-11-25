#include <gtest/gtest.h>
#include <qsd/alpha_covariance_calculator.h>
#include <qsd/colored_gaussian_process_generator.h>
#include <qsd/colored_gaussian_process_path.h>
#include <protein_chain/correlation_function_lorentzians.h>
#include <boost/make_shared.hpp>
#include <math/stats/covariance.h>
#include <math/average.h>

using namespace rql::math::stats;

class ColoredGaussianProcessGeneratorTest: public testing::Test
{
};

TEST_F(ColoredGaussianProcessGeneratorTest,Test)
{
	const size_t nbr_steps = 4;
	const size_t dim = 2;
	
	const double dt = 0.5;
	std::vector<boost::shared_ptr<const CorrelationFunction> > alpha(dim, boost::make_shared<CorrelationFunctionLorentzians>(1, 0.25, 0.5));
	alpha[1] = boost::make_shared<CorrelationFunctionLorentzians>(0.1, 0.5, 0.25);
	ColoredGaussianProcessGenerator generator(alpha, nbr_steps, dt);
	AlphaCovarianceCalculator acc(alpha, nbr_steps, dt);

	const double cov_tol = 1E-14;
	std::vector<Eigen::MatrixXd> cov_matrices(dim);
	// Check covariance matrices
	for (size_t m = 0; m < dim; ++m) {
		SingleAlphaCovarianceCalculator sacc(acc.single_alpha_calculator(m));
		Eigen::MatrixXd& cov_real = cov_matrices[m];
		Eigen::VectorXcd tmp(nbr_steps);
		ColoredGaussianProcessGenerator::calculate_covariance_matrix(sacc, cov_real);
		ASSERT_TRUE(cov_real == cov_real.adjoint()) << "Symmetric " << m;
		for (size_t i = 0; i < nbr_steps; ++i) {
			const size_t i_idx = 2*i;
			const double aa_ii = cov_real(i_idx, i_idx);
			const double bb_ii = cov_real(i_idx+1, i_idx+1);
			const double ab_ii = cov_real(i_idx, i_idx+1);
			ASSERT_EQ(aa_ii, bb_ii);
			ASSERT_EQ(ab_ii, 0.0);
			ASSERT_NEAR(aa_ii, sacc.covariance_RR_II(i, i), cov_tol);
			for (size_t j = 0; j < i; ++j) {
				const size_t j_idx = 2*j;
				const double aa_ij = cov_real(i_idx, j_idx);
				const double bb_ij = cov_real(i_idx+1, j_idx+1);
				ASSERT_EQ(aa_ij, bb_ij);
				ASSERT_NEAR(aa_ij, sacc.covariance_RR_II(i, j), cov_tol);
				const double ab_ij = cov_real(i_idx, j_idx + 1);
				const double ba_ij = cov_real(i_idx + 1, j_idx);
				ASSERT_EQ(ab_ij, - ba_ij);
				ASSERT_NEAR(ab_ij, sacc.covariance_RI(i, j), cov_tol);
			}
		}
	}

	// Check simulation
	const size_t nbr_paths = 25000;
	ColoredGaussianProcessGenerator::path_type path(generator.buildPath());
	ColoredGaussianProcessGenerator::Workspace wksp(generator.workspace());
	const size_t cov_dim = 2 * nbr_steps * dim;
	Covariance ccalc(cov_dim);
	std::vector<rql::math::Average<double> > means(cov_dim, rql::math::Average<double>());
	std::vector<double> data(cov_dim);
	for (size_t n = 0; n < nbr_paths; ++n) {
		generator.simulate(wksp, path);
		for (size_t m = 0; m < dim; ++m) {
			for (size_t i = 0; i < nbr_steps; ++i) {
				const size_t idx = m * 2 * nbr_steps + 2 * i;
				data[idx] = path(i, m).real();
				data[idx + 1] = path(i, m).imag();
			}
		}
		ccalc.update(data.begin(), data.end());
		for (size_t i = 0; i < cov_dim; ++i)
			means[i].update(data[i]);
	}	

	Eigen::MatrixXd covariance(ccalc.covariance());		
	const double tol = 1E-2;
	for (size_t m1 = 0; m1 < dim; ++m1)
		for (size_t i1 = 0; i1 < nbr_steps; ++i1) {
			const size_t idx1 = m1*2*nbr_steps + 2*i1;
			for (size_t m2 = 0; m2 < dim; ++m2)
				for (size_t i2 = 0; i2 < nbr_steps; ++i2) {
					const size_t idx2 = m2*2*nbr_steps + 2*i2;
					if (m1 != m2) {
						for (size_t k1 = 0; k1 < 2; ++k1)
							for (size_t k2 = 0; k2 < 2; ++k2)
								ASSERT_NEAR(covariance(idx1 + k1, idx2 + k2), 0.0, tol);
					} else {
						for (size_t k1 = 0; k1 < 2; ++k1)
							for (size_t k2 = 0; k2 < 2; ++k2)
								ASSERT_NEAR(covariance(idx1 + k1, idx2 + k2), cov_matrices[m1](2*i1+k1, 2*i2+k2), tol);
					}
				}
		}

	for (size_t i = 0; i < cov_dim; ++i)
		ASSERT_NEAR(0.0, means[i].value(), 1E-10) << i; // low tolerance for antithetic sampling
}

TEST_F(ColoredGaussianProcessGeneratorTest,TestCommonAlpha)
{
	const size_t nbr_steps = 4;
	const size_t dim = 2;	
	const double dt = 0.5;
	std::vector<boost::shared_ptr<const CorrelationFunction> > alpha(dim, boost::make_shared<CorrelationFunctionLorentzians>(1, 0.25, 0.5));
	ColoredGaussianProcessGenerator generator1(alpha, nbr_steps, dt);
	ColoredGaussianProcessGenerator generator2(alpha[0], dim, nbr_steps, dt);
	const size_t nbr_paths = 20;
	ColoredGaussianProcessGenerator::path_type path1(generator1.buildPath());
	ColoredGaussianProcessGenerator::Workspace wksp1(generator1.workspace());
	ColoredGaussianProcessGenerator::path_type path2(generator2.buildPath());
	ColoredGaussianProcessGenerator::Workspace wksp2(generator2.workspace());
	for (size_t n = 0; n < nbr_paths; ++n) {
		generator1.simulate(wksp1, path1);
		generator2.simulate(wksp2, path2);
		for (size_t m = 0; m < dim; ++m) {
			for (size_t i = 0; i < nbr_steps; ++i) {
				ASSERT_EQ(path1(i, m), path2(i, m));
			}
		}
	}
}

