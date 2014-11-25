#ifndef __RANDOM_VECTOR_H
#define __RANDOM_VECTOR_H

#include <vector>
#include <Eigen/Core>
#include "core.h"

class RandomVector
{
public:
	PROTEIN_CHAIN_API virtual ~RandomVector();
	PROTEIN_CHAIN_API_DEBUG virtual unsigned int dim() const = 0;
	PROTEIN_CHAIN_API_DEBUG virtual void draw(std::vector<double>& vec) = 0;
	PROTEIN_CHAIN_API_DEBUG virtual void draw(Eigen::VectorXd& vec) = 0;
	PROTEIN_CHAIN_API_DEBUG virtual void set_seed(unsigned int seed) = 0;
};

#endif // __RANDOM_VECTOR_H
