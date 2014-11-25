#include "strunz_calculator_grid_explicit.h"
#include <cassert>
#include <stdexcept>
#include <protein_chain/CorrelationFunction.h>

StrunzCalculatorGridExplicit::StrunzCalculatorGridExplicit(const std::vector<boost::shared_ptr<const CorrelationFunction> >& alphas, const Eigen::MatrixXcd& Hel, size_t nbrTimePoints, double dt)
	: m_alphas(alphas), m_Hel(Hel), m_nbr_time_pts(nbrTimePoints), m_dt(dt), m_nbr_sites(Hel.rows()), m_alpha_cache(m_nbr_sites, m_nbr_time_pts + 1)
{
	for (size_t t = 0; t <= m_nbr_time_pts; ++t)
		for (size_t m = 0; m < m_nbr_sites; ++m)
			m_alpha_cache(m, t) = (*m_alphas[m])(t*m_dt);
}


void StrunzCalculatorGridExplicit::calculateEffectiveHamiltonianTimesMinusI(const Workspace& wksp, Eigen::MatrixXcd& Heff_times_minus_I) const
{
	assert( &wksp.m_owner == this );
	if (wksp.m_time_idx >= m_nbr_time_pts)
		throw std::domain_error("Time index too large");
	Heff_times_minus_I = m_Hel;
	Heff_times_minus_I *= std::complex<double>(0,-1);
	if (wksp.m_time_idx > 0) {
		for (size_t si = 0; si <= wksp.m_time_idx; ++si) {
			const double weight = ((si == 0 || si == wksp.m_time_idx) ? 0.5 : 1) * m_dt;
			const size_t tau_idx = wksp.m_time_idx - si;
			for (size_t m = 0; m < m_nbr_sites; ++m) {
				const std::complex<double> wm = m_alpha_cache(m, tau_idx) * weight;
				const_map_type d2op = wksp.stateConst(si, m);
				for (size_t n = 0; n < m_nbr_sites; ++n)
					Heff_times_minus_I(m, n) += wm*d2op(m, n);
			}
		}
	}
}

void StrunzCalculatorGridExplicit::step(Workspace& wksp) const
{
	assert( &wksp.m_owner == this );
	calculateEffectiveHamiltonianTimesMinusI(wksp, wksp.m_HeffTimesMinusI);
	const size_t op_bnd_idx = wksp.operatorIdx(wksp.m_time_idx + 1, 0);
	const Eigen::MatrixXcd& heff_times_minus_i = wksp.m_HeffTimesMinusI;
	Eigen::MatrixXcd& f = wksp.m_f;
	for (size_t i = 0; i < op_bnd_idx; ++i) {
		const_map_type d2i = wksp.stateConst(i);
		f.noalias() = heff_times_minus_i * d2i;
		f.noalias() -= d2i * heff_times_minus_i;
		wksp.state(i) += f*m_dt;
	}
	++wksp.m_time_idx;
}

// Nested classes

StrunzCalculatorGridExplicit::Workspace::Workspace(const StrunzCalculatorGridExplicit& owner, size_t nbrTimePts, size_t nbrSites)
	: m_owner(owner), m_HeffTimesMinusI(nbrSites, nbrSites), m_nbr_time_pts(nbrTimePts), m_nbr_sites(nbrSites), m_nbr_operators(nbrTimePts * nbrSites), m_matrix_size(nbrSites*nbrSites), m_data((m_nbr_operators)*m_matrix_size), m_f(nbrSites,nbrSites)
	, m_d1_data_begin(m_data.get() + m_nbr_operators*m_matrix_size)
{
	reset();
}

void StrunzCalculatorGridExplicit::Workspace::reset()
{
	m_time_idx = 0;
	std::fill(m_data.get(), m_data.get() + m_data.size(), std::complex<double>(0.0,0.0));
	for (size_t i = 0; i < m_nbr_time_pts; ++i)
		for (size_t m = 0; m < m_nbr_sites; ++m)
			state(operatorIdx(i, m))(m, m) = -1;
}
