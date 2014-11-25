#ifndef __EVOLVER_SPARSE_H
#define __EVOLVER_SPARSE_H

#include "evolver_stepping.h"
#include "core.h"
#include <boost/shared_ptr.hpp>
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>

template <class M> class EvolverSparse: public EvolverStepping
{
public:
	//! Multiplicator multiplies by -i*H
	EvolverSparse(double timeStep, M multiplicator, size_t dim);
	EvolverSparse(double timeStep);
	void step(double step_size, Eigen::VectorXcd& state);
private:
	M m_multiplicator;
	Eigen::VectorXcd m_f;
};

class EigenMultiplicator
{
public:
	typedef Eigen::SparseMatrix<std::complex<double> > matrix_type;
	typedef boost::shared_ptr<const matrix_type> matrix_pointer;
	typedef Eigen::VectorXcd vector_type;
	PROTEIN_CHAIN_API EigenMultiplicator();
	PROTEIN_CHAIN_API_DEBUG EigenMultiplicator(matrix_pointer hamTimesMinusI);
	PROTEIN_CHAIN_API EigenMultiplicator(const Eigen::SparseMatrix<std::complex<double> >& hamTimesMinusI);
	PROTEIN_CHAIN_API void operator()(const vector_type& state, vector_type& f) const;
private:
	matrix_pointer m_matrix;
};

template <class M> EvolverSparse<M>::EvolverSparse(double timeStep, M multiplicator, size_t dim)
	: EvolverStepping(timeStep), m_multiplicator(multiplicator), m_f(dim)
{
}

template <class M> EvolverSparse<M>::EvolverSparse(double timeStep)
	: EvolverStepping(timeStep)
{
}

template <class M> void EvolverSparse<M>::step(double step_size, Eigen::VectorXcd& state)
{
	m_multiplicator(state, m_f);
	state += m_f*step_size;
}

#endif // __EVOLVER_SPARSE_H
