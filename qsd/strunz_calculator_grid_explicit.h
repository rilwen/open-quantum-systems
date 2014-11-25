#ifndef __STRUNZ_CALCULATOR_GRID_EXPLICIT_H
#define __STRUNZ_CALCULATOR_GRID_EXPLICIT_H

#include <cassert>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <Eigen/Core>
#include <math/MathUtils.h>
#include "core.h"

class CorrelationFunction;

class StrunzCalculatorGridExplicit
{
private:
	typedef Eigen::Map<Eigen::MatrixXcd,Eigen::Aligned> map_type;
	typedef Eigen::Map<const Eigen::MatrixXcd,Eigen::Aligned> const_map_type;
public:	
	class Workspace
	{
	public:
	  Workspace(Workspace&&) = default;
	private:
		friend class StrunzCalculatorGridExplicit;
		//		friend class StrunzSimulatorGrid;
		QSD_API_DEBUG Workspace(const StrunzCalculatorGridExplicit& owner, size_t nbrTimePts, size_t nbrSites);		
		size_t operatorIdx(size_t timeIdx, size_t siteIdx) const { return timeIdx * m_nbr_sites + siteIdx; }		
		std::complex<double>* statePtr(size_t operatorIdx) { return m_data.get() + operatorIdx * m_matrix_size; }
		const std::complex<double>* statePtr(size_t operatorIdx) const { return m_data.get() + operatorIdx * m_matrix_size; }
		map_type state(size_t operatorIdx) { return map_type(statePtr(operatorIdx), m_nbr_sites, m_nbr_sites); }
		const_map_type stateConst(size_t operatorIdx) const { return const_map_type(statePtr(operatorIdx), m_nbr_sites, m_nbr_sites); }
		const_map_type stateConst(size_t timeIdx, size_t siteIdx) const { return stateConst(operatorIdx(timeIdx, siteIdx)); }
		void reset();

		const StrunzCalculatorGridExplicit& m_owner;		
		Eigen::MatrixXcd m_HeffTimesMinusI;
		const size_t m_nbr_time_pts;
		const size_t m_nbr_sites;
		const size_t m_nbr_operators;
		const size_t m_matrix_size;		
		rql::math::eigen_aligned_data<std::complex<double> > m_data;
		size_t m_time_idx;
		Eigen::MatrixXcd m_f;
		const std::complex<double>* m_d1_data_begin;
	};

	QSD_API_DEBUG StrunzCalculatorGridExplicit(const std::vector<boost::shared_ptr<const CorrelationFunction> >& alphas, const Eigen::MatrixXcd& Hel, size_t nbrTimePoints, double dt);
	Workspace workspace() const { return Workspace(*this, m_nbr_time_pts, m_nbr_sites); }
	size_t nbrSites() const { return m_nbr_sites; }
	size_t nbrTimePts() const { return m_nbr_time_pts; }
	size_t nbrOperators() const { return nbrSites() * m_nbr_time_pts; }
	double dt() const { return m_dt; }
	//! Calculates the effective hamiltonian for t=timeIdx*dt()
	QSD_API_DEBUG void step(Workspace& wksp) const;
	//! Returns current effective hamiltonian from the workspace
	const Eigen::MatrixXcd& effectiveHamiltonianTimesMinusI(const Workspace& wksp) const { return wksp.m_HeffTimesMinusI; }
	//! Reset workspace
	void reset(Workspace& wksp) const { wksp.reset(); }
	void calculateEffectiveHamiltonianTimesMinusI(const Workspace& wksp, Eigen::MatrixXcd& heff) const;	
	const Eigen::MatrixXcd& Hel() const { return m_Hel; }
private:	
	StrunzCalculatorGridExplicit& operator=(const StrunzCalculatorGridExplicit&); // not implemented
private:
	std::vector<boost::shared_ptr<const CorrelationFunction> > m_alphas;	
	Eigen::MatrixXcd m_Hel;
	const size_t m_nbr_time_pts;
	const double m_dt;	
	const size_t m_nbr_sites;
	Eigen::MatrixXcd m_alpha_cache; // stores precalculated values of alpha correlation functions
};


#endif // __STRUNZ_CALCULATOR_GRID_EXPLICIT_H
