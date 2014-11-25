#ifndef __PSEUDOMODE_CORE_H
#define __PSEUDOMODE_CORE_H

#ifdef _MSC_VER
#ifdef BUILDING_PSEUDOMODE
#define PSEUDOMODE_API __declspec(dllexport)
#ifndef NDEBUG
#define PSEUDOMODE_API_DEBUG __declspec(dllexport)
#else // NDEBUG
#define PSEUDOMODE_API_DEBUG
#endif // NDEBUG
#else // BUILDING_PSEUDOMODE
#define PSEUDOMODE_API __declspec(dllimport)
#ifndef NDEBUG
#define PSEUDOMODE_API_DEBUG __declspec(dllimport)
#else // NDEBUG
#define PSEUDOMODE_API_DEBUG
#endif // NDEBUG
#endif // BUILDING_PSEUDOMODE
#else // _MSC_VER
#define PSEUDOMODE_API
#define PSEUDOMODE_API_DEBUG
#endif // _MSC_VER

#endif // __PSEUDOMODE_CORE_H
