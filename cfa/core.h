#ifndef __CFA_CORE_H
#define __CFA_CORE_H

#ifdef _MSC_VER
#ifdef BUILDING_CFA
#define CFA_API __declspec(dllexport)
#ifndef NDEBUG
#define CFA_API_DEBUG __declspec(dllexport)
#else // NDEBUG
#define CFA_API_DEBUG
#endif // NDEBUG
#else // BUILDING_CFA
#define CFA_API __declspec(dllimport)
#ifndef NDEBUG
#define CFA_API_DEBUG __declspec(dllimport)
#else // NDEBUG
#define CFA_API_DEBUG
#endif // NDEBUG
#endif // BUILDING_CFA
#else // _MSC_VER
#define CFA_API
#define CFA_API_DEBUG
#endif // _MSC_VER

#endif // __CFA_CORE_H
