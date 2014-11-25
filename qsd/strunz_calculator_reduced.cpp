#include "strunz_calculator_reduced.h"
#include <cassert>
#include <sstream>
#include <stdexcept>
#include <iostream>

StrunzCalculatorReduced::StrunzCalculatorReduced(const CorrelationFunctionsDecomposition& corr_funcs_decomp, const Eigen::MatrixXcd& Hel)
	: m_decomp(corr_funcs_decomp), m_nbr_sites(Hel.rows()), m_HelTimesMinusI(Hel * std::complex<double>(0,-1))
{
	if (Hel.rows() != Hel.cols())
		throw std::domain_error("Hel matrix must be square");
	if (corr_funcs_decomp.nbr_sites() != m_nbr_sites)
		throw std::domain_error("Wrong number of alphas");
	m_nbr_operators = corr_funcs_decomp.total_nbr_peaks();
}

void StrunzCalculatorReduced::add_interaction_hamiltonian_times_minus_i(const Data& state, Eigen::MatrixXcd& h) const
{
	assert(h.rows() == static_cast<int>(m_nbr_sites));
	assert(h.cols() == static_cast<int>(m_nbr_sites));
	for (size_t m = 0; m < m_nbr_sites; ++m) {
		const size_t curr_peak_nbr = m_decomp.nbr_peaks(m);
		for (size_t k = 0; k < curr_peak_nbr; ++k)
			h.row(m) += state(m, k).row(m);
	}
}

void StrunzCalculatorReduced::calculate_interaction_hamiltonian(const Data& state, Eigen::MatrixXcd& h) const
{
	h.setZero(m_nbr_sites, m_nbr_sites);
	add_interaction_hamiltonian_times_minus_i(state, h);
	h *= std::complex<double>(0.0, 1.0);
}

void StrunzCalculatorReduced::calculate_effective_hamiltonian_times_minus_i(const StrunzCalculatorReduced::Data& state, Eigen::MatrixXcd& h) const
{
	h = m_HelTimesMinusI;
	add_interaction_hamiltonian_times_minus_i(state, h);
	/*for (size_t m = 0; m < m_nbr_sites; ++m) {
		const size_t curr_peak_nbr = m_decomp.nbr_peaks(m);
		for (size_t k = 0; k < curr_peak_nbr; ++k)
			h.row(m) += state(m, k).row(m);
	}*/
}

void StrunzCalculatorReduced::calculate_time_derivative(const StrunzCalculatorReduced::Data& state, const Eigen::MatrixXcd& h_eff_times_minus_i, StrunzCalculatorReduced::Data& derivative) const
{
	for (size_t m = 0; m < m_nbr_sites; ++m) {
		const size_t curr_peak_nbr = m_decomp.nbr_peaks(m);
		for (size_t k = 0; k < curr_peak_nbr; ++k) {
			const_map_type st_op(state(m, k));
			map_type dr_op(derivative(m, k));
			dr_op.noalias() = h_eff_times_minus_i*st_op;
			dr_op.noalias() -= st_op*h_eff_times_minus_i;
			/*if (std::abs(dr_op.trace()) > 1E-2) {
				std::stringstream ss;
				ss << "Trace too large: " << dr_op.trace() << " with norms: " << h_eff_times_minus_i.norm() << " and " << st_op.norm();
				throw std::runtime_error(ss.str().c_str());
			}*/
			dr_op += st_op*m_decomp.exponent(m, k);
			dr_op(m,m) -= m_decomp.scale(m, k);
		}
	}
}

std::complex<double> StrunzCalculatorReduced::operator_trace_theoretical(size_t site_idx, size_t peak_idx, double time) const
{
	const std::complex<double> w = m_decomp.exponent(site_idx, peak_idx);
	if (std::abs(w) <= 1E-8) // avoids numerical issues
		return - m_decomp.scale(site_idx, peak_idx) * time * exp(w*time);
	else
		return (1.0 - exp(w*time))*m_decomp.scale(site_idx, peak_idx) / w;
}

double StrunzCalculatorReduced::total_absolute_trace_deviation_from_theoretical(const Data& state, double time) const
{
	double sum = 0;
	for (size_t m = 0; m < m_nbr_sites; ++m) {
		const size_t curr_peak_nbr = m_decomp.nbr_peaks(m);
		for (size_t k = 0; k < curr_peak_nbr; ++k) {
			sum += std::abs( operator_trace_numerical(state, m, k) - operator_trace_theoretical(m, k, time) );
		}
	}
	return sum;
}

void StrunzCalculatorReduced::normalize_traces(Data& state, double time) const
{
	for (size_t m = 0; m < m_nbr_sites; ++m) {
		const size_t curr_peak_nbr = m_decomp.nbr_peaks(m);
		for (size_t k = 0; k < curr_peak_nbr; ++k) {
			const std::complex<double> actual_trace = operator_trace_numerical(state, m, k);
			if (actual_trace != 0.0) {
				const std::complex<double> theoretical_trace = operator_trace_theoretical(m, k, time);
				state(m, k) *= (theoretical_trace / actual_trace);
			}
		}
	}
}

void StrunzCalculatorReduced::calculate_operator(const Data& state, const size_t site_idx, Eigen::MatrixXcd& m) const
{
	m.setZero(m_nbr_sites, m_nbr_sites);
	const size_t curr_peak_nbr = m_decomp.nbr_peaks(site_idx);
	for (size_t k = 0; k < curr_peak_nbr; ++k) {
		m += state(site_idx, k);
	}
	
	/*
	// check commutator condition
	double comm_cond_dev = 0;
	for (size_t p = 0; p < m_nbr_sites; ++p) {
		if (p != site_idx) {
			comm_cond_dev += fabs( m(site_idx, p) ) + fabs( m(p, site_idx) ); 
		}
	}
	std::cout << "Commutator deviation: " << comm_cond_dev << std::endl;
	*/
}

// Nested classes

StrunzCalculatorReduced::Data::Data(const StrunzCalculatorReduced& owner)
	: m_matrix_dim(owner.m_nbr_sites)
{
	std::vector<size_t> pc(owner.m_decomp.peak_counts());
	m_ptrs = rql::math::Jagged2DArray<std::complex<double>*>(pc.begin(), pc.end());
	const size_t matrix_size = m_matrix_dim*m_matrix_dim;
	m_data = boost::shared_ptr<store_type>(new store_type(matrix_size*owner.m_nbr_operators));
	std::complex<double>* ptr = m_data->get();
	for (size_t m = 0; m < owner.m_nbr_sites; ++m) {
		const size_t curr_nbr_peaks = owner.m_decomp.nbr_peaks(m);
		for (size_t k = 0; k < curr_nbr_peaks; ++k) {
			m_ptrs(m, k) = ptr;
			ptr += matrix_size;
		}
	}
	assert(ptr == m_data->end());
	clear();
}

void StrunzCalculatorReduced::Data::clear()
{
	std::fill(m_data->begin(), m_data->end(), std::complex<double>(0,0));
}

StrunzCalculatorReduced::Data& StrunzCalculatorReduced::Data::copy_data_from(const Data& other)
{
	assert(m_data->size() == other.m_data->size());
	std::copy(other.m_data->begin(), other.m_data->end(), m_data->begin());
	return *this;
}

StrunzCalculatorReduced::Data& StrunzCalculatorReduced::Data::operator*=(double x)
{
	const store_type::iterator end = m_data->end();
	for (store_type::iterator it = m_data->begin(); it != end; ++it)
		(*it) *= x;
	return *this;
}

StrunzCalculatorReduced::Data& StrunzCalculatorReduced::Data::operator/=(double x)
{
	const store_type::iterator end = m_data->end();
	for (store_type::iterator it = m_data->begin(); it != end; ++it)
		(*it) /= x;
	return *this;
}

StrunzCalculatorReduced::Data& StrunzCalculatorReduced::Data::operator+=(const Data& other)
{
	assert(m_data->size() == other.m_data->size());
	const store_type::iterator end = m_data->end();
	store_type::const_iterator oit = other.m_data->begin();
	for (store_type::iterator it = m_data->begin(); it != end; ++it) {
		(*it) += *oit;
		++oit;
	}
	return *this;
}

StrunzCalculatorReduced::Data& StrunzCalculatorReduced::Data::add(const Data& a, const double x)
{
	assert(m_data->size() == a.m_data->size());
	const store_type::iterator end = m_data->end();
	store_type::const_iterator ait = a.m_data->begin();
	for (store_type::iterator it = m_data->begin(); it != end; ++it) {
		(*it) += (*ait)*x;
		++ait;
	}
	return *this;
}
