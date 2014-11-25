#ifndef __STRUNZ_CALCULATOR_REDUCED_H
#define __STRUNZ_CALCULATOR_REDUCED_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <math/MathUtils.h>
#include <math/Jagged2DArray.h>
#include "correlation_functions_decomposition.h"
#include "core.h"

class CorrelationFunctionLorentzians;

//! Reduced equations for spectral density functions decomposed into sums of Lorentzians
class StrunzCalculatorReduced
{
private:
	typedef Eigen::Map<Eigen::MatrixXcd,Eigen::Aligned> map_type;
	typedef Eigen::Map<const Eigen::MatrixXcd,Eigen::Aligned> const_map_type;
public:	
	class Data
	{
	public:				
		map_type operator()(const size_t site_idx, const size_t peak_idx)
		{
			return map_type(m_ptrs(site_idx, peak_idx), m_matrix_dim, m_matrix_dim);
		}
		const_map_type operator()(const size_t site_idx, const size_t peak_idx) const
		{
			return const_map_type(m_ptrs(site_idx, peak_idx), m_matrix_dim, m_matrix_dim);
		}
		QSD_API_DEBUG void clear();
		QSD_API_DEBUG Data& copy_data_from(const Data& other);
		//! add a*x
		QSD_API_DEBUG Data& add(const Data& a, double x);
		QSD_API_DEBUG Data& operator*=(double x);
		QSD_API_DEBUG Data& operator/=(double x);
		QSD_API_DEBUG Data& operator+=(const Data& other);
	private:	
		QSD_API_DEBUG Data(const StrunzCalculatorReduced& owner);
		friend class StrunzCalculatorReduced;
		typedef rql::math::eigen_aligned_data<std::complex<double> > store_type;
		rql::math::Jagged2DArray<std::complex<double>*> m_ptrs;
		boost::shared_ptr<store_type> m_data;
		size_t m_matrix_dim;
	};
public:	
	QSD_API_DEBUG StrunzCalculatorReduced(const CorrelationFunctionsDecomposition& corr_funcs_decomp, const Eigen::MatrixXcd& Hel);
	size_t nbrSites() const { return m_nbr_sites; }	
	Data data() const { return Data(*this); }
	QSD_API_DEBUG void calculate_effective_hamiltonian_times_minus_i(const Data& state, Eigen::MatrixXcd& h) const;
	void calculate_interaction_hamiltonian(const Data& state, Eigen::MatrixXcd& h) const;
	QSD_API_DEBUG void calculate_time_derivative(const Data& state, const Eigen::MatrixXcd& h_eff_times_minus_i, Data& derivative) const;
	QSD_API_DEBUG void calculate_effective_hamiltonian_times_minus_i_and_time_derivative(const Data& state, Eigen::MatrixXcd& h_eff_times_minus_i, Data& derivative) const
	{
		calculate_effective_hamiltonian_times_minus_i(state, h_eff_times_minus_i);
		calculate_time_derivative(state, h_eff_times_minus_i, derivative);
	}
	void calculate_operator(const Data& state, size_t site_idx, Eigen::MatrixXcd& m) const;
	std::complex<double> operator_trace_theoretical(size_t site_idx, size_t peak_idx, double time) const;
	std::complex<double> operator_trace_numerical(const Data& state, size_t site_idx, size_t peak_idx) const { return state(site_idx, peak_idx).trace(); }
	//! sum_m sum_k |TraceTheoretical(D_mk) - TraceNumerical(D_mk)|
	QSD_API_DEBUG double total_absolute_trace_deviation_from_theoretical(const Data& state, double time) const;
	//! Rescale operators so that they have correct traces.
	QSD_API_DEBUG void normalize_traces(Data& state, double time) const;
	const CorrelationFunctionsDecomposition& correlation_functions_decomposition() const { return m_decomp; }
private:
	void add_interaction_hamiltonian_times_minus_i(const Data& state, Eigen::MatrixXcd& h) const;
private:
	std::vector<boost::shared_ptr<const CorrelationFunctionLorentzians> > m_alphas;
	size_t m_nbr_sites;
	Eigen::MatrixXcd m_HelTimesMinusI;
	size_t m_nbr_operators;
	CorrelationFunctionsDecomposition m_decomp;
};

#endif // __STRUNZ_CALCULATOR_REDUCED_H
