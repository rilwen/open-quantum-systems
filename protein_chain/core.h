#ifndef __PROTEIN_CHAIN_CORE_H
#define __PROTEIN_CHAIN_CORE_H

#ifdef _MSC_VER
#ifdef BUILDING_PROTEIN_CHAIN
#define PROTEIN_CHAIN_API __declspec(dllexport)
#ifndef NDEBUG
#define PROTEIN_CHAIN_API_DEBUG __declspec(dllexport)
#else // NDEBUG
#define PROTEIN_CHAIN_API_DEBUG
#endif // NDEBUG
#else // BUILDING_PROTEIN_CHAIN
#define PROTEIN_CHAIN_API __declspec(dllimport)
#ifndef NDEBUG
#define PROTEIN_CHAIN_API_DEBUG __declspec(dllimport)
#else // NDEBUG
#define PROTEIN_CHAIN_API_DEBUG
#endif // NDEBUG
#endif // BUILDING_PROTEIN_CHAIN
#else // _MSC_VER
#define PROTEIN_CHAIN_API
#define PROTEIN_CHAIN_API_DEBUG
#endif // _MSC_VER

#endif // __PROTEIN_CHAIN_CORE_H
