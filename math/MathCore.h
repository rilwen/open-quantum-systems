#ifndef __MATH_CORE_H
#define __MATH_CORE_H

#ifdef _MSC_VER
#ifdef BUILDING_RQL_MATH
#define RQL_MATH_API __declspec(dllexport)
#ifndef NDEBUG
#define RQL_MATH_API_DEBUG __declspec(dllexport)
#else // NDEBUG
#define RQL_MATH_API_DEBUG
#endif // NDEBUG
#else // BUILDING_RQL_MATH
#define RQL_MATH_API __declspec(dllimport)
#ifndef NDEBUG
#define RQL_MATH_API_DEBUG __declspec(dllimport)
#else // NDEBUG
#define RQL_MATH_API_DEBUG
#endif // NDEBUG
#endif // BUILDING_RQL_MATH
#else // _MSC_VER
#define RQL_MATH_API
#define RQL_MATH_API_DEBUG
#endif // _MSC_VER

#endif // __MATH_CORE_H