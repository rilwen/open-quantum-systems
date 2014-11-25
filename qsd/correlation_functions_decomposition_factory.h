#ifndef __QSD_CORRELATION_FUNCTIONS_DECOMPOSITION_FACTORY_H
#define __QSD_CORRELATION_FUNCTIONS_DECOMPOSITION_FACTORY_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include "core.h"

class CorrelationFunctionsDecomposition;
class CorrelationFunctionLorentzians;
class CorrelationFunctionDecomposable;

namespace CorrelationFunctionsDecompositionFactory
{
	QSD_API_DEBUG CorrelationFunctionsDecomposition decomposeLorentzians(const std::vector<boost::shared_ptr<const CorrelationFunctionLorentzians> >& alphas);
	QSD_API_DEBUG CorrelationFunctionsDecomposition decomposeVirtual(const std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >& alphas);
}

#endif // __QSD_CORRELATION_FUNCTIONS_DECOMPOSITION_FACTORY_H
