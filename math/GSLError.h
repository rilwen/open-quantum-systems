#ifndef __MATH_GSL_ERROR_H
#define __MATH_GSL_ERROR_H

#include <stdexcept>

namespace rql {
	class GSLError: public std::runtime_error
	{
	public:
		GSLError(const char* msg, int errorCode);
		int errorCode() const { return m_error_code; }
	private:
		int m_error_code;
	};
}

#endif // __MATH_GSL_ERROR_H
