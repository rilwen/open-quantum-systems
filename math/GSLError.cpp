#include "GSLError.h"

namespace rql {
	GSLError::GSLError(const char* msg, int errorCode)
		: std::runtime_error(msg), m_error_code(errorCode)
	{
	}
}
