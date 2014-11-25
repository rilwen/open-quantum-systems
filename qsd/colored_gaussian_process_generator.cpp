#include "colored_gaussian_process_generator.h"
#include "colored_gaussian_process_path.h"
#include <protein_chain/CorrelationFunction.h>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <iostream>
#include <map>
#include <protein_chain/random_vector_gaussian.h>
#include <protein_chain/random_vector_correlated.h>
#include <boost/make_shared.hpp>
#include "alpha_covariance_calculator.h"
#include <Eigen/LU>

ColoredGaussianProcessGenerator::ColoredGaussianProcessGenerator(const std::vector<boost::shared_ptr<const CorrelationFunction> >& alpha, const size_t nbrPts, double dt)
	: m_alpha(alpha), m_nbr_pts(nbrPts), m_dim(alpha.size()), m_dt(dt), m_transforms(alpha.size())
{
	initialize(false);
}

ColoredGaussianProcessGenerator::ColoredGaussianProcessGenerator(const boost::shared_ptr<const CorrelationFunction>& alpha, size_t dim, const size_t nbrPts, double dt)
	: m_alpha(dim, alpha), m_nbr_pts(nbrPts), m_dim(dim), m_dt(dt), m_transforms(dim)
{
	initialize(true);
}

//! Converts NxN matrices Cov(a_k, a_l), Cov(a_k, b_l) into one 2N x 2N covariance matrix Cov(y_i, y_j)
class CovarianceInformationProvider
{
public:
	CovarianceInformationProvider(const SingleAlphaCovarianceCalculator sacc)
		: m_sacc(sacc)
	{
	}

	double operator()(size_t i, size_t j) const
	{
		if (i % 2 == 0) {
			// first variable is a_{i/2}
			if (j % 2 == 0) {
				// second variable is a_{j/2}
				return m_sacc.covariance_RR_II(i/2, j/2);
			} else {
				// second variable is b_{j/2}
				return m_sacc.covariance_RI(i/2, j/2);
			}
		} else {
			// first variable is b_{i/2}
			if (j % 2 == 0) {
				// second variable is a_{j/2}
				return m_sacc.covariance_RI(j/2, i/2);
			} else {
				// second variable is b_{j/2}
				return m_sacc.covariance_RR_II(i/2, j/2);
			}
		}
	}
private:
	CovarianceInformationProvider& operator=(const CovarianceInformationProvider&); // not implemented; declaration added to make the compiler shut up and not give a warning
	// unsafe!
	const SingleAlphaCovarianceCalculator m_sacc;
};

void ColoredGaussianProcessGenerator::initialize(bool common_alpha)
{
	if (!m_nbr_pts)
		throw std::domain_error("At least 1 step needed");
	if (m_dim == 0)
		throw std::domain_error("At least 1 site needed");
	
	AlphaCovarianceCalculator acc(common_alpha ? AlphaCovarianceCalculator(m_alpha[0], m_dim, m_nbr_pts, m_dt) : AlphaCovarianceCalculator(m_alpha, m_nbr_pts, m_dt));
	Eigen::VectorXcd avg_alpha(m_nbr_pts);	
	const size_t nbr_vars = 2*acc.nbr_time_pts();
	if (!common_alpha) {
		// share transforms for identical correlation functions
		typedef std::map<boost::shared_ptr<const CorrelationFunction>, boost::shared_ptr<const Eigen::MatrixXd> > transform_map_type;
		transform_map_type transform_map;
		for (size_t m = 0; m < m_dim; ++m) {
			transform_map_type::const_iterator it = transform_map.find(m_alpha[m]);
			if (it == transform_map.end()) {
				CovarianceInformationProvider cov_info_prov(acc.single_alpha_calculator(m));
				m_transforms[m] = //calculate_transform(acc, m);
					RandomVectorCorrelated::build_transform(cov_info_prov, nbr_vars);
				transform_map[m_alpha[m]] = m_transforms[m];
			} else {
				m_transforms[m] = it->second;
			}
			assert(m_transforms[m]);
		}
	} else {
		CovarianceInformationProvider cov_info_prov(acc.single_alpha_calculator(0));
		m_transforms[0] = //calculate_transform(acc, 0);
			RandomVectorCorrelated::build_transform(cov_info_prov, nbr_vars);
		assert(m_transforms[0]);
		for (size_t m = 1; m < m_dim; ++m)
			m_transforms[m] = m_transforms[0];		
	}	
}

void ColoredGaussianProcessGenerator::calculate_covariance_matrix(const SingleAlphaCovarianceCalculator& sacc, Eigen::MatrixXd& covariance)
{
	CovarianceInformationProvider cov_info_prov(sacc);
	const size_t nbr_vars = 2*sacc.nbr_time_pts();
	covariance.resize(nbr_vars, nbr_vars);		
	for (size_t i = 0; i < nbr_vars; ++i) {
		covariance(i, i) = cov_info_prov(i, i);
		for (size_t j = 0; j < i; ++j) {
			const double covariance_ij = cov_info_prov(i, j);
			covariance(i, j) = covariance(j, i) = covariance_ij;
		}
	}

	if (covariance != covariance.adjoint())
		throw std::runtime_error("Error constructing covariance matrix");
}

ColoredGaussianProcessGenerator::Workspace ColoredGaussianProcessGenerator::workspace() const
{
	return Workspace(*this, m_nbr_pts, m_dim);
}

ColoredGaussianProcessGenerator::path_type ColoredGaussianProcessGenerator::buildPath() const
{
	return path_type(m_nbr_pts, m_dim);
}

typedef Eigen::Block<Eigen::MatrixXd, 2, 2> block_type;

class BlockTransform
{
public:
	BlockTransform(Eigen::MatrixXd& matrix);
	block_type operator()(size_t r, size_t c) { return m_matrix.block<2,2>(2*r, 2*c); }
private:
	size_t offset(size_t r, size_t c);
	Eigen::MatrixXd& m_matrix;
};

BlockTransform::BlockTransform(Eigen::MatrixXd& matrix)
	: m_matrix(matrix)
{
}

boost::shared_ptr<const Eigen::MatrixXd> ColoredGaussianProcessGenerator::calculate_transform(const AlphaCovarianceCalculator& acc, const size_t m) const
{
	const size_t nbr_vars = 2*acc.nbr_time_pts();
	const boost::shared_ptr<Eigen::MatrixXd> transform = boost::shared_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd(nbr_vars, nbr_vars));
	BlockTransform tr(*transform);
	transform->setZero();
	if (!nbr_vars) {
		return transform;
	}

	double norm_prev;
	for (size_t k = 0; k < acc.nbr_time_pts(); ++k) {
		double norm_diff2 = 0;
		double norm_curr = 0;
		for (size_t j = 0; j <= k; ++j) {
			Eigen::Matrix2d tmp;
			tmp(0, 0) = tmp(1, 1) = acc.covariance_RR_II(m, k, j);
			tmp(0, 1) = acc.covariance_RI(m, k, j);
			tmp(1, 0) = - tmp(0, 1);
			//std::cout << "covariance(" << k << ", " << j << ") == " << tmp << std::endl;
			const double cov_norm = tmp.norm();
			for (size_t j2 = 0; j2 < j; ++j2) {
				tmp.noalias() -= tr(k, j2) * tr(j, j2).transpose();
			}
			if (j < k) {
				// solve tmp = A_kj * A_jj^T for A_kj
				Eigen::Matrix2d AjjT = tr(j, j).transpose();
				tr(k, j) = tmp * AjjT.inverse();
			} else {
				// solve tmp = A_kj*A_kj^T
				// j == k  ==> tmp^T = tmp
				// hence A_kj = A_kj^T
				//std::cout << "tmp == " << tmp << std::endl;
				const double trace = tmp.trace();
				const double det = tmp.determinant();				
				if (det < 0) {
					std::stringstream ss;
					ss << "Error: negative determinant: " << det;
					throw std::runtime_error(ss.str().c_str());
				}
				const double s = sqrt(det);			
				const double sqr_t = trace + 2*s;
				if (sqr_t <= 0) {
					std::cout << "k == " << k << std::endl;
					std::cout << "j == " << j << std::endl;
					std::cout << "CovRR == " << acc.covariance_RR_II(m, k, j) << "\n";
					std::cout << "CovRI == " << acc.covariance_RI(m, k, j) << "\n";
					std::cout << "tmp == " << tmp << std::endl;
					/*std::stringstream ss;
					ss << "Error: trace + 2 * sqrt(det) <= 0: " << sqr_t;
					throw std::runtime_error(ss.str().c_str());*/
					tr(k, k).setZero();
				} else {
					const double t = sqrt(sqr_t);
					tr(k, k) = tmp;
					tr(k, k)(0, 0) += s;
					tr(k, k)(1, 1) += s;
					tr(k, k) /= t;
				}
			}
			std::cout << "tr(" << k << "," << j << ") == " << tr(k, j) << std::endl;
			tmp.noalias() -= tr(k, j) * tr(j, j).transpose();
			if (tmp.norm() > 1E-8*cov_norm) {
				std::stringstream ss;
				ss << "Error: did not solve properly: " << tmp.norm();
				throw std::runtime_error(ss.str().c_str());
			}
			norm_curr += pow(tr(k, j).norm(), 2);
			if (j > 0) {
				assert(j <= k);
				tmp = tr(k, j) - tr(k - 1, j - 1);
				norm_diff2 += pow(tmp.norm(), 2);
			}
		}
		norm_curr = sqrt(norm_curr);
		if (k > 10 && norm_diff2 < 1E-12*norm_curr*norm_prev) {
			std::cout << "Copying forward at k == " << k << std::endl;
			for (size_t k2 = k + 1; k2 < acc.nbr_time_pts(); ++k2) {
				const size_t offs = k2 - k;
				for (size_t j = 0; j <= k; ++j) {
					tr(k2, j + offs) = tr(k, j);
				}
			}
			break;
		}
		std::cout << "norms: " << norm_curr << " " << norm_prev << " " << norm_diff2 << "\n";
		norm_prev = norm_curr;		
		std::cout << "Finished k == " << k << std::endl;
	}

	return transform;
}

// Nested class

ColoredGaussianProcessGenerator::Workspace::Workspace(const ColoredGaussianProcessGenerator& owner, const size_t nbrPts, const size_t dim)
	: m_owner(owner), m_random_vectors(dim), m_ab(2*nbrPts)
{
	static const bool USE_ANTITHETIC_SAMPLING = true;
	const boost::shared_ptr<RandomVector> iid = boost::make_shared<RandomVectorGaussian>(2*nbrPts, 1.0, dim, USE_ANTITHETIC_SAMPLING);
	for (size_t m = 0; m < dim; ++m)
		m_random_vectors[m] = boost::make_shared<RandomVectorCorrelated>(iid, owner.m_transforms[m]);		
}
