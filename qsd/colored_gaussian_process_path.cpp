#include "colored_gaussian_process_path.h"
#include <complex>
#include <Eigen/Core>



template <class T> ColoredGaussianProcessPath<T>::ColoredGaussianProcessPath(size_t len, size_t dim)
	: m_size(len), m_dim(dim), m_data(len, dim)
{
	m_data.setZero();
}

template class ColoredGaussianProcessPath<double>;
template class ColoredGaussianProcessPath<std::complex<double> >;
