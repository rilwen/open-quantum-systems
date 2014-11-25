#ifndef __COLORED_GAUSSIAN_PROCESS_GENERATOR_H
#define __COLORED_GAUSSIAN_PROCESS_GENERATOR_H

#include <complex>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <protein_chain/random_vector.h>
#include "core.h"

class AlphaCovarianceCalculator;
class SingleAlphaCovarianceCalculator;
template <class T> class ColoredGaussianProcessPath;
class CorrelationFunction;

class ColoredGaussianProcessGenerator
{
public:
	typedef ColoredGaussianProcessPath<std::complex<double> > path_type;
	class Workspace
	{
	private:
		Workspace(const ColoredGaussianProcessGenerator& owner, size_t nbrPts, size_t dim);
		friend class ColoredGaussianProcessGenerator;
		const ColoredGaussianProcessGenerator& m_owner;
		std::vector<boost::shared_ptr<RandomVector> > m_random_vectors;
		Eigen::VectorXd m_ab;
	};
	QSD_API_DEBUG ColoredGaussianProcessGenerator(const std::vector<boost::shared_ptr<const CorrelationFunction> >& alpha, size_t nbrPts, double dt);
	QSD_API_DEBUG ColoredGaussianProcessGenerator(const boost::shared_ptr<const CorrelationFunction>& alpha, size_t dim, size_t nbrPts, double dt);
	QSD_API_DEBUG Workspace workspace() const;
	//! path(timeIdx, siteIdx) is the path value for timeIdx-th time point and siteIdx-th dimension
	template <class M> void simulate(Workspace& wksp, M& pathMatrix) const;
	QSD_API_DEBUG path_type buildPath() const;
	QSD_API_DEBUG static void calculate_covariance_matrix(const SingleAlphaCovarianceCalculator& sacc, Eigen::MatrixXd& covariance);
private:
	void initialize(bool common_alpha);
	boost::shared_ptr<const Eigen::MatrixXd> calculate_transform(const AlphaCovarianceCalculator& acc, const size_t m) const;

	std::vector<boost::shared_ptr<const CorrelationFunction> > m_alpha;
	size_t m_nbr_pts;
	size_t m_dim;
	double m_dt;
	std::vector<boost::shared_ptr<const Eigen::MatrixXd> > m_transforms;
};

template <class M> void ColoredGaussianProcessGenerator::simulate(ColoredGaussianProcessGenerator::Workspace& wksp, M& path) const
{
	for (size_t m = 0; m < m_dim; ++m) {
		wksp.m_random_vectors[m]->draw(wksp.m_ab); 
		for (size_t i = 0; i < m_nbr_pts; ++i)
			path(i, m) = std::complex<double>(wksp.m_ab[2*i], wksp.m_ab[2*i+1]);
	}
}

#endif // __COLORED_GAUSSIAN_PROCESS_GENERATOR_H
