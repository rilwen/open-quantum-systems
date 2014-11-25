#ifndef __EIGENSTATE_PROPERTIES_H
#define __EIGENSTATE_PROPERTIES_H

#include <cmath>

struct EigenstateProperties
{
	template <class V> static double pr(const V& eigenstate);
	template <class V> static double mu(const V& eigenstate);
	template <class V> static double l1norm(const V& eigenstate);
};

template <class V> double EigenstateProperties::pr(const V& eigenstate)
{
	const unsigned int l = eigenstate.size();
	double sum = 0;
	for (unsigned int i = 0; i < l; ++i) {
		const double x = eigenstate[i];
		sum += x*x*x*x;
	}
	return 1/sum;
}

template <class V> double EigenstateProperties::mu(const V& eigenstate)
{
	const unsigned int l = eigenstate.size();
	double sum = 0;
	for (unsigned int i = 0; i < l; ++i) {
		const double x = eigenstate[i];
		sum += x;
	}
	return sum;
}

template <class V> double EigenstateProperties::l1norm(const V& eigenstate)
{
	const unsigned int l = eigenstate.size();
	double sum = 0;
	for (unsigned int i = 0; i < l; ++i) {
		const double x = eigenstate[i];
		sum += std::abs(x);
	}
	return sum;
}

#endif // __EIGENSTATE_PROPERTIES_H
