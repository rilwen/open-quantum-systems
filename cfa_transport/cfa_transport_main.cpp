#include <Eigen/Core>
#include <boost/make_shared.hpp>
#include <protein_chain/absorption_spectrum_calculator_time_dependent.h>
#include <protein_chain/correlation_function_lorentzians.h>
#include <protein_chain/command_line_arguments_reader.h>
#include <protein_chain/hamiltonian_properties.h>
#include <protein_chain/hamiltonian_factory.h>
#include <protein_chain/rho_analyzer.h>
#include <protein_chain/exact_hamiltonian.h>
#include <cfa/classical_field_approximation.h>
#include <cfa/matrix_field_approximation.h>
#include <cfa/matrix_field_approximation_reduced.h>
#include <cfa/matrix_field_approximation_higher_order.h>
#include <cfa/reduced_operator_approximation_lorentzians.h>
#include <cfa/reduced_operator_approximation_lorentzians_higher_order.h>
#include <cfa/cfa_processor_density_matrix.h>
#include <cfa/cfa_processor_pair.h>
#include <cfa/cfa_processor_mean_energy.h>
#include <cfa/cfa_simulator.h>
#include <cfa/cfa_utils.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>

static const int RING = 0;
static const int CHAIN = 1;

typedef cfa::MatrixFieldApproximation method_type;
//typedef cfa::MatrixFieldApproximationHigherOrder method_type;
//typedef cfa::ReducedOperatorApproximationLorentzians method_type;
//typedef cfa::ReducedOperatorApproximationLorentziansHigherOrder method_type;
static const int hamiltonian_type = CHAIN;

boost::shared_ptr<const CorrelationFunctionLorentzians> buildFMOAlpha(const double omega0, const double Gamma, const double scale)
{
	/*
	const size_t nbrpeaks = 7;
	Eigen::VectorXd omega0_vec(nbrpeaks);
	Eigen::VectorXd hwhm(nbrpeaks);
	Eigen::VectorXd height(nbrpeaks);
	if (hamiltonian_type != RING) {
		std::cerr << "Set hamiltonian type to RING!" << std::endl;
		return 0;
	}
	std::cout << "RING_PRL: Ignoring omega0, Gamma and scale..." << std::endl;
	omega0_vec[0] = 0.16; hwhm[0] = 0.06; height[0] = 0.5;
	omega0_vec[1] = 0.43; hwhm[1] = 0.05; height[1] = 0.43;
	omega0_vec[2] = 0.86; hwhm[2] = 0.05; height[2] = 0.7;
	omega0_vec[3] = 1.21; hwhm[3] = 0.05; height[3] = 0.7;
	omega0_vec[4] = 1.38; hwhm[4] = 0.06; height[4] = 1.9;
	omega0_vec[5] = 1.56; hwhm[5] = 0.07; height[5] = 1.6;
	omega0_vec[6] = 0.05; hwhm[6] = 0.04; height[6] = 0.12;
	*/

	//height /= 10;

	//Eigen::VectorXd omega0(2);
	//Eigen::VectorXd hwhm(2);
	//Eigen::VectorXd height(2);
	//omega0[0] = 1; hwhm[0] = 0.1; height[0] = CorrelationFunctionLorentzians::height(0.64, 0.1);
	//omega0[1] = 0.3; hwhm[1] = 0.05; height[1] = CorrelationFunctionLorentzians::height(0.025, 0.05);

	const size_t nbrPeaks = 1;
	Eigen::VectorXd omega0_vec(nbrPeaks);
	Eigen::VectorXd hwhm(nbrPeaks);
	Eigen::VectorXd height(nbrPeaks);
	
	omega0_vec[0] = omega0;
	hwhm[0] = Gamma;
	height[0] = CorrelationFunctionLorentzians::height(scale, Gamma);
	
	/*omega0[0] = 1;
	hwhm[0] = 0.25;
	height[0] = 0.5;*/
	 
	/*omega0[0] = 1;
	hwhm[0] = 0.25;
	height[0] = CorrelationFunctionLorentzians::height(0.64, hwhm[0]);*/

	/*omega0[0] = 0.5;
	omega0[1] = 1.5;
	hwhm.fill(0.1);
	height.fill(CorrelationFunctionLorentzians::height(1, hwhm[0]));*/

	/*omega0[0] = 1;
	hwhm.fill(0.5);
	height.fill(CorrelationFunctionLorentzians::height(1, hwhm[0]));*/

/*	Eigen::VectorXd omega0(2);
	Eigen::VectorXd hwhm(2);
	Eigen::VectorXd height(2);
	omega0[0] = 0.5;
	omega0[1] = 1;
	hwhm.fill(1E-8);
	height[0] = CorrelationFunctionLorentzians::height(0.25, hwhm[0]);*/

	return boost::make_shared<CorrelationFunctionLorentzians>(omega0_vec, hwhm, height);
}

void prepare_single_mode(Eigen::VectorXd& omega, Eigen::MatrixXcd& g, const size_t nbr_sites, double omega0, double g0)
{
	omega.resize(1);
	omega[0] = omega0;
	g.resize(1, nbr_sites);
	g.fill(-g0);
	g(0,0) = g0;
}

void prepare_two_modes(Eigen::VectorXd& omega, Eigen::MatrixXcd& g, const size_t nbr_sites, double omega0, double omega1, double g0)
{
	omega.resize(2);
	omega[0] = omega0;
	omega[1] = omega1;
	g.resize(2, nbr_sites);
	g.fill(-g0);
	g.col(0).fill(g0);
}

void prepare_one_mode_independent_baths(Eigen::VectorXd& omega, Eigen::MatrixXcd& g, const size_t nbr_sites, double omega0, double g0)
{
	static const unsigned int nbr_modes_per_site = 1;
	omega.resize(nbr_modes_per_site);
	omega[0] = omega0;
	ExactHamiltonian::independentBaths(g0, nbr_modes_per_site, nbr_sites, g, omega);
	g.transposeInPlace();
}

void prepare_two_modes_independent_baths(Eigen::VectorXd& omega, Eigen::MatrixXcd& g, const size_t nbr_sites, double omega0, double omega1, double g0)
{
	static const unsigned int nbr_modes_per_site = 2;
	omega.resize(nbr_modes_per_site);
	omega[0] = omega0;
	omega[1] = omega1;
	ExactHamiltonian::independentBaths(g0, nbr_modes_per_site, nbr_sites, g, omega);
	g.transposeInPlace();
}

void prepare_three_modes_independent_baths(Eigen::VectorXd& omega, Eigen::MatrixXcd& g, const size_t nbr_sites, double omega0, double omega1, double omega2, double g0)
{
	static const unsigned int nbr_modes_per_site = 3;
	omega.resize(nbr_modes_per_site);
	omega[0] = omega0;
	omega[1] = omega1;
	omega[2] = omega2;
	ExactHamiltonian::independentBaths(g0, nbr_modes_per_site, nbr_sites, g, omega);
	g.transposeInPlace();
}

void min_max_omega_for_discretization(const boost::shared_ptr<const CorrelationFunctionLorentzians> alpha, double& w0, double& w1, double f)
{
	const size_t n = alpha->omegas().size();
	assert(n);
	double max_hwhm = 0;
	w0 = w1 = alpha->omegas()[0];
	for (size_t i = 0; i < n; ++i) {
		max_hwhm = std::max(max_hwhm, alpha->gammas()[i]);
		const double w = alpha->omegas()[i];
		w0 = std::min(w0, w);
		w1 = std::max(w1, w);
	}
	w0 -= f*max_hwhm;
	//w0 = std::max(w0, 0.0);	
	w1 += f*max_hwhm;
}

int main(int argc, char* argv[])
{
	const time_t t0 = clock();

	size_t nbr_sites, init;
	size_t nbr_bath_modes_per_site;
	double dt, Tmax;
	double C, f;
	double kappa;
	double omega0 = 0;
	double hwhm = 0;
	double scale = 0;
	try {		
		CommandLineArgumentsReader reader(argc, argv);
		reader.pop(nbr_sites);
		reader.pop(init);
		//reader.pop(nbr_bath_modes_per_site);
		reader.pop(C);
		reader.pop(dt);
		reader.pop(Tmax);
		/*reader.pop(f);
		if (reader.hasNext()) {
			reader.pop(kappa);
		} else {
			kappa = 0;
		}*/
		reader.pop(omega0);
		//reader.pop(hwhm);
		reader.pop(scale);
		//reader.pop(omega1);
		//reader.pop(g0);
	} catch (...) {
		//std::cerr << "Usage: " << argv[0] << " <nbr_sites> <nbr_bath_modes_per_site> <dt> <t_max> <omega0> <omega1> <g0>" << std::endl;
		//std::cerr << "Usage: " << argv[0] << " <nbr_sites> <init_site> <nbr_bath_modes_per_site> <C> <dt> <t_max> <f> [kappa]" << std::endl;
		std::cerr << "Usage: " << argv[0] << " <nbr_sites> <init_site> <C> <dt> <t_max> <omega> <scale>" << std::endl;
		return -1;
	}

	double J;
	switch (hamiltonian_type) {
		case RING:
			J = HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(nbr_sites, C);
			break;
		case CHAIN:
			J = C/2;
			break;
		default:
			std::cerr << "Unknown hamiltonian type: " << hamiltonian_type << std::endl;
			return -1;
	};

	std::cout << "J == " << J << std::endl;
	
	dt /= std::abs(C);
	Tmax /= std::abs(C);

	std::cout << "Tmax: " << Tmax << std::endl;

	try {		
		const size_t nbr_steps = static_cast<size_t>(floor((Tmax + 0.5*dt)/dt));
		const size_t nbr_bath_modes = nbr_sites * nbr_bath_modes_per_site;
		double w0;
		double w1;		
		Eigen::MatrixXcd H0(nbr_sites, nbr_sites);		
		
		switch (hamiltonian_type) {
			case RING:
				HamiltonianFactory::buildRing(std::vector<double>(nbr_sites, 0.0), J, H0);
				break;
			case CHAIN:
				HamiltonianFactory::buildChain(std::vector<double>(nbr_sites, 0.0), J, 0.0, H0);
				break;
			default:
				std::cerr << "Unknown hamiltonian type: " << hamiltonian_type << std::endl;
				return -1;
		};
		
		Eigen::VectorXd omega;
		Eigen::MatrixXcd g;
		prepare_one_mode_independent_baths(omega, g, nbr_sites, omega0, sqrt(scale));
		std::cout << "g: " << g << std::endl;
		/*if (omega1 == 0) {
			prepare_one_mode_independent_baths(omega, g, nbr_sites, omega0, g0);
		} else {
			prepare_two_modes_independent_baths(omega, g, nbr_sites, omega0, omega1, g0);
		}
		std::cout << omega << std::endl;
		std::cout << g << std::endl;*/
		//const boost::shared_ptr<const CorrelationFunctionLorentzians> alpha = buildFMOAlpha(omega0, hwhm, scale);
		//min_max_omega_for_discretization(alpha, w0, w1, f);
		//std::cerr << "w0 == " << w0 << "\n";
		//std::cerr << "w1 == " << w1 << "\n";
		//cfa::discretise_correlation_function(alpha, w0, w1, nbr_sites, nbr_bath_modes_per_site, omega, g);
		//std::cerr << omega.minCoeff() << " " << omega.maxCoeff() << "\n";

		cfa::Simulator<method_type> simulator(method_type(H0, g, omega)); // for general MFA methods
		//cfa::Simulator<method_type> simulator(method_type(H0, alpha)); // for Lorentzian method
		//cfa::Simulator<method_type> simulator(method_type(H0, cfa::discretise_correlation_function(alpha, w0, w1, nbr_bath_modes_per_site))); // for Lorentzian method -- discrettize to peaks

		Eigen::VectorXcd initial_state(nbr_sites);
		initial_state.setZero();
		initial_state[init] = 1;
		//initial_state.fill(1/sqrt(static_cast<double>(nbr_sites)));		
		cfa::ProcessorDensityMatrix<method_type> proc_rho(initial_state, nbr_steps);
		cfa::ProcessorMeanEnergy<method_type> proc_mean_energy(initial_state, nbr_steps);
		auto joined_procs = cfa::join_processors(proc_rho, proc_mean_energy);
		simulator.simulate(initial_state, dt, nbr_steps, joined_procs);
		
		std::string outf_name;
		{
			std::stringstream ss;
			ss << "roa" << method_type::name() << "_omega" << omega0 << "_Gamma" << hwhm << "_scale" << scale << "_N" << nbr_sites << "_dt" << dt << "_tMax" << Tmax << "_C" << C << "_J" << J << ".txt";	
			//ss << "roa_RingPRL" << method_type::name() << "_N" << nbr_sites << "_dt" << dt << "_tMax" << Tmax << "_C" << C << "_J" << J << ".txt";
			outf_name = ss.str();
		}
		std::ofstream outf(outf_name, std::fstream::out);
		std::cout << "Writing to file: " << outf_name << std::endl;

		// for H0 * rho product
		Eigen::MatrixXcd rho_times_H0(nbr_sites, nbr_sites);
		
		for (size_t i = 0; i <= nbr_steps; ++i) {
			const double t = i*dt;
			const double t_in_PRL_units = t* std::abs(C);
			Eigen::MatrixXcd rho(proc_rho.rho(i));			
			const std::complex<double> trace = rho.trace();
			rho *= proc_rho.rho_trace(i) / trace; // normalisation
			
			// Average system energy
			rho_times_H0.noalias() = H0 * rho;
			const double energy = rho_times_H0.trace().real();
			
			try {
				RhoAnalyzer rho_analyzer(rho, method_type::WRAP_RHO_EIGENVALUES, method_type::RHO_EIGENVALUE_TOLERANCE);			
				//rho_analyzer.correct_rho(rho);
				/*for (size_t site_idx = 0; site_idx < nbr_sites; ++site_idx) {
					if (rho(site_idx,site_idx).real() > 1) {
						rho(site_idx, site_idx) = 0;
					}
				}*/			
				outf << t;//t_in_PRL_units;
				for (size_t site_idx = 0; site_idx < nbr_sites; ++site_idx) {
					outf << "\t" << std::setprecision(16) << rho(site_idx, site_idx).real();
				}			
				outf << "\t" << rho_analyzer.vonNeumannEntropy();
				outf << "\t" << rho_analyzer.coherence();
				outf << "\t" << trace.real();
				outf << "\t" << proc_rho.rho_trace(i).real();
				outf << "\t" << (rho - rho.adjoint()).norm();
				outf << "\t" << energy;
			} catch (std::exception& e) {
				std::stringstream ss;
				ss << "Error at time " << t_in_PRL_units << ": ";
				ss << e.what();
				throw std::runtime_error(ss.str().c_str());
			}
			outf << "\n";
		}

		outf.close();

	} catch (std::exception& e) {
		std::cerr << "Exception: " << e.what() << std::endl;
		return -1;
	}


	const time_t t1 = clock();
	std::cout << "Running time [s]: " << (t1 - t0) / static_cast<double>(CLOCKS_PER_SEC) << std::endl;
}
