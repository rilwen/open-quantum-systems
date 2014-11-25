#include "evolver_sparse.h"
#include <boost/make_shared.hpp>

EigenMultiplicator::EigenMultiplicator(matrix_pointer hamTimesMinusI)
	: m_matrix(hamTimesMinusI)
{
}

EigenMultiplicator::EigenMultiplicator(const Eigen::SparseMatrix<std::complex<double> >& hamTimesMinusI)
	: m_matrix(boost::make_shared<matrix_type>(hamTimesMinusI))
{
}

EigenMultiplicator::EigenMultiplicator()
	: m_matrix(new matrix_type())
{
}

void EigenMultiplicator::operator()(const vector_type& state, vector_type& f) const
{
	f = (*m_matrix) * state;
}

template class EvolverSparse<EigenMultiplicator>;
