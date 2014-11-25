#include "random_vector_correlated.h"
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <boost/make_shared.hpp>
#include <math/stats/covariance_decomposition_calculator_iterative.h>

RandomVectorCorrelated::RandomVectorCorrelated(boost::shared_ptr<RandomVector> iidGenerator, const Eigen::MatrixXd& covariance)
	: m_iid_generator(iidGenerator), m_iid_draw(covariance.rows())
{
	assert( iidGenerator->dim() == static_cast<size_t>(covariance.rows()) );	
	m_transform = build_transform(covariance);
}

RandomVectorCorrelated::RandomVectorCorrelated(boost::shared_ptr<RandomVector> iidGenerator, const boost::shared_ptr<const Eigen::MatrixXd>& transform)
	: m_iid_generator(iidGenerator), m_transform(transform)
{
	if (!transform)
		throw std::domain_error("No transform provided");
	m_iid_draw.resize(transform->rows());
}

boost::shared_ptr<const Eigen::MatrixXd> RandomVectorCorrelated::build_transform(const Eigen::MatrixXd& covariance)
{
	assert( covariance == covariance.adjoint() );	
	boost::shared_ptr<Eigen::MatrixXd> transform = boost::make_shared<Eigen::MatrixXd>(covariance.rows(), covariance.cols());
	rql::math::stats::CovarianceDecompositionCalculatorIterative transform_calculator;
	transform_calculator.decompose(covariance, *transform);
	//check_transform(covariance, covariance.rows(), *transform, 1E-8);
	return transform;
}

void RandomVectorCorrelated::draw(std::vector<double>& vec)
{
	compute();
	assert(vec.size() >= dim());
	std::copy(m_draw.data(), m_draw.data() + dim(), vec.begin());
}

void RandomVectorCorrelated::draw(Eigen::VectorXd& vec)
{
	compute();
	vec = m_draw;
}

void RandomVectorCorrelated::compute()
{
	m_iid_generator->draw(m_iid_draw);
	m_draw.noalias() = (*m_transform) * m_iid_draw;
}

void RandomVectorCorrelated::set_seed(unsigned int seed)
{
	m_iid_generator->set_seed(seed);
}
