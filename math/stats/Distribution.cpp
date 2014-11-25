#include "Distribution.h"
#include "../MathCore.h"
#include <cstdio>

namespace rql { namespace math { namespace stats {

	bool prob_is_ok(double prob)
	{
		if (prob < - Distribution<double>::prob_epsilon()) {
			fprintf(stderr, "Probability %g lower than -DISTRIBUTION_PROB_EPS\n", prob);
			return false;
		}
		if (prob > 1+Distribution<double>::prob_epsilon()) {
			fprintf(stderr, "Probability %g higher than 1+DISTRIBUTION_PROB_EPS\n", prob);
			return false;
		}
		return true;
	}

	bool probs_are_ok(const std::vector<double>& probabilities)
	{
		double sum = 0;
		for (std::vector<double>::const_iterator i = probabilities.begin(); i != probabilities.end(); ++i) {
			assert(prob_is_ok(*i));
			sum += *i;
		}
		if (std::abs(sum - 1) >= Distribution<double>::prob_epsilon()) {
			fprintf(stderr, "Sum of probabilities is %g while it should be 1\n", sum);
			fprintf(stderr, "Whole distribution:\n");
			if (probabilities.size() < 20)
				for (std::vector<double>::const_iterator i = probabilities.begin(); i != probabilities.end(); ++i) {
					fprintf(stderr, "%g\n", *i);
				}
			return false;
		}
		else return true;
	}

}}}