#ifndef __RANDOM_VECTOR_CORRELATED_H
#define __RANDOM_VECTOR_CORRELATED_H

#include "random_vector.h"
#include "core.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <Eigen/Core>
#include <math/stats/covariance_decomposition_calculator_iterative.h>

class RandomVectorCorrelated: public RandomVector
{
public:
	//! @param[in] iidGenerator Generates iid numbers with sigma = 1
	PROTEIN_CHAIN_API RandomVectorCorrelated(boost::shared_ptr<RandomVector> iidGenerator, const Eigen::MatrixXd& covariance);
	PROTEIN_CHAIN_API RandomVectorCorrelated(boost::shared_ptr<RandomVector> iidGenerator, const boost::shared_ptr<const Eigen::MatrixXd>& transform);
	PROTEIN_CHAIN_API_DEBUG unsigned int dim() const { return m_iid_generator->dim(); }
	PROTEIN_CHAIN_API_DEBUG void draw(std::vector<double>& vec);
	PROTEIN_CHAIN_API_DEBUG void draw(Eigen::VectorXd& vec);
	PROTEIN_CHAIN_API_DEBUG const boost::shared_ptr<const Eigen::MatrixXd>& transform() const { return m_transform; }
	PROTEIN_CHAIN_API_DEBUG static boost::shared_ptr<const Eigen::MatrixXd> build_transform(const Eigen::MatrixXd& covariance);
	template <class M> static boost::shared_ptr<const Eigen::MatrixXd> build_transform(const M& covariance, size_t dimension);
	PROTEIN_CHAIN_API_DEBUG void set_seed(unsigned int seed);
private:
	void compute();
	template <class M> static void check_transform(const M& covariance, size_t dimension, const Eigen::MatrixXd& transform, double tolerance);
	boost::shared_ptr<RandomVector> m_iid_generator;
	Eigen::VectorXd m_iid_draw;
	Eigen::VectorXd m_draw;
	boost::shared_ptr<const Eigen::MatrixXd> m_transform;
};

template <class M> boost::shared_ptr<const Eigen::MatrixXd> RandomVectorCorrelated::build_transform(const M& covariance, const size_t dimension)
{
	boost::shared_ptr<Eigen::MatrixXd> transform = boost::make_shared<Eigen::MatrixXd>(dimension, dimension);
	rql::math::stats::CovarianceDecompositionCalculatorIterative transform_calculator;
	transform_calculator.decompose(covariance, dimension, *transform);
	//check_transform(covariance, dimension, *transform, 1E-8);
	return transform;
}

template <class M> void RandomVectorCorrelated::check_transform(const M& covariance, const size_t dimension, const Eigen::MatrixXd& transform, const double tolerance)
{
	for (size_t i = 0; i < dimension; ++i) {
		for (size_t j = 0; j < dimension; ++j) {
			const double resulting_covariance = transform.row(i).dot(transform.row(j));
			const double expected_covariance = covariance(i, j);
			const double effective_tolerance = expected_covariance != 0 ? tolerance * std::abs(expected_covariance) : tolerance;
			if (std::abs(resulting_covariance - expected_covariance) > effective_tolerance)
				throw std::runtime_error("Transform calculation inaccurate");
		}
	}
}

#endif // __RANDOM_VECTOR_CORRELATED_H
