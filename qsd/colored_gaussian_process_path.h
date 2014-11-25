#ifndef __COLORED_GAUSSIAN_PROCESS_PATH_H
#define __COLORED_GAUSSIAN_PROCESS_PATH_H

#include <cassert>
#include <Eigen/Core>

template <class T> struct ColoredGaussianProcessPathTraits
{
};

template <> struct ColoredGaussianProcessPathTraits<double>
{
	typedef Eigen::MatrixXd storage_type;
};

template <> struct ColoredGaussianProcessPathTraits<std::complex<double> >
{
	typedef Eigen::MatrixXcd storage_type;
};

template <class T> class ColoredGaussianProcessPath
{
public:
	ColoredGaussianProcessPath(size_t len, size_t dim);
	size_t size() const { return m_size; }
	size_t dim() const { return m_dim; }
	const T& operator()(size_t stepIdx, size_t dimIdx) const { return m_data(stepIdx, dimIdx); }
	T& operator()(size_t stepIdx, size_t dimIdx) { return m_data(stepIdx, dimIdx); }
private:
	size_t m_size;
	size_t m_dim;
	typename ColoredGaussianProcessPathTraits<T>::storage_type m_data;
};


#endif // __COLORED_GAUSSIAN_PROCESS_PATH_H
