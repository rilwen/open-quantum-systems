#ifndef __UTILS_H
#define __UTILS_H

#include <complex>
#include <sstream>
#include <iosfwd>
#include <vector>

template <class Out, class In> Out strConvert(const In& i)
{
	std::stringstream ss;
	ss << i;
	Out o;
	ss >> o;
	return o;
}

template <class M> void dumpMatrix(const M& matrix, std::ostream& stream)
{
	const unsigned int r = matrix.rows();
	const unsigned int c = matrix.cols();
	for (unsigned int i = 0; i < r; ++i) {
		for (unsigned int j = 0; j < c; ++j) {
			stream << matrix(i,j);
			if (j + 1 != c) {
				stream << " ";
			}
		}
		stream << "\n";
	}
}

template <class T> void dumpStdVector(const std::vector<T>& vector, std::ostream& stream)
{
	const unsigned int l = vector.size();
	for (unsigned int i = 0; i < l; ++i) {
		if (i>0) stream << " ";
		stream << vector[i];
	}
	stream << "\n";
}

inline std::complex<double>& setReal(std::complex<double>& z, double r)
{
#ifdef _WIN32
	z.real(r);
#else
	z.real(r);
#endif
	return z;
}

inline std::complex<double>& setImag(std::complex<double>& z, double i)
{
#ifdef _WIN32
	z.imag(i);
#else
	z.imag(i);
#endif
	return z;
}

template <class T>
void appendParam(std::stringstream& ss, const char* name, const T& value)
{
	ss << "_" << name << value;
}

#endif // __UTILS_H
