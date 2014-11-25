#ifndef __QSD_CORE_H
#define __QSD_CORE_H

#ifdef _MSC_VER
#ifdef BUILDING_QSD
#define QSD_API __declspec(dllexport)
#ifndef NDEBUG
#define QSD_API_DEBUG __declspec(dllexport)
#else // NDEBUG
#define QSD_API_DEBUG
#endif // NDEBUG
#else // BUILDING_QSD
#define QSD_API __declspec(dllimport)
#ifndef NDEBUG
#define QSD_API_DEBUG __declspec(dllimport)
#else // NDEBUG
#define QSD_API_DEBUG
#endif // NDEBUG
#endif // BUILDING_QSD
#else // _MSC_VER
#define QSD_API
#define QSD_API_DEBUG
#endif // _MSC_VER

#endif // __QSD_CORE_H
