#include "Interpolator.h"
#include "InterpolatorImpl.h"

namespace rql {
	namespace interp {

		Interpolator::Interpolator(std::shared_ptr<InterpolatorImpl> impl)
			: Pimpl<InterpolatorImpl>(impl)
		{
		}		
	}
}