#include <ctime>
#include <iomanip>
#include <iostream>
#include <vector>
#include <boost/make_shared.hpp>
#include <Eigen/Core>
#include <math/average.h>
#include <math/Jagged2DArray.h>
#include <math/MathUtils.h>
#include <protein_chain/command_line_arguments_reader.h>
#include <protein_chain/correlation_function_lorentzians.h>
#include <protein_chain/correlation_function_modes.h>
#include <protein_chain/exact_hamiltonian.h>
#include <protein_chain/evolver_taylor_expansion.h>
#include <protein_chain/hamiltonian_properties.h>
#include <protein_chain/hamiltonian_factory.h>
#include <protein_chain/rho_analyzer.h>
#include <qsd/strunz_simulator_reduced_stochastic_nonlinear.h>
#include <sstream>
#include <fstream>

boost::shared_ptr<const CorrelationFunctionLorentzians> buildFMOAlpha(const double scale, const double gamma, const double omega0)
{
	/*Eigen::VectorXd omega0(5);
	Eigen::VectorXd hwhm(5);
	Eigen::VectorXd height(5);
	omega0[0] = 0.005; hwhm[0] = 0.005; height[0] = 1E-2;
	omega0[1] = 0.02; hwhm[1] = 0.008; height[1] = 4E-2;
	omega0[2] = 0.05; hwhm[2] = 0.01; height[2] = 2E-2;
	omega0[3] = 0.1; hwhm[3] = 0.02; height[3] = 1.2E-2;
	omega0[4] = 0.17; hwhm[4] = 0.055; height[4] = 7.9E-2;*/

	/*
	//// PRL
	const size_t npeaks = 7;
	Eigen::VectorXd omega(npeaks);
	Eigen::VectorXd hwhm(npeaks);
	Eigen::VectorXd height(npeaks);
	omega[0] = 0.16; hwhm[0] = 0.06; height[0] = 0.5;
	omega[1] = 0.43; hwhm[1] = 0.05; height[1] = 0.43;
	omega[2] = 0.86; hwhm[2] = 0.05; height[2] = 0.7;
	omega[3] = 1.21; hwhm[3] = 0.05; height[3] = 0.7;
	omega[4] = 1.38; hwhm[4] = 0.06; height[4] = 1.9;
	omega[5] = 1.56; hwhm[5] = 0.07; height[5] = 1.6;
	omega[6] = 0.05; hwhm[6] = 0.04; height[6] = 0.12;
	*/

	Eigen::VectorXd omega(1);
	Eigen::VectorXd hwhm(1);
	Eigen::VectorXd height(1);

	omega[0] = omega0;
	hwhm[0] = gamma;
	height[0] = CorrelationFunctionLorentzians::height(scale, gamma);
	
	//NARROW PEAK
	/*omega[0] = 1;
	hwhm[0] = 0.125;
	height[0] = 1;//CorrelationFunctionLorentzians::height(0.86, 0.05);*/

	//WIDE PEAK
	/*omega[0] = 1;
	hwhm[0] = 0.05;
	height[0] = CorrelationFunctionLorentzians::height(0.64, hwhm[0]);*/

/*	omega[0] = 1;
	hwhm[0] = 0.1;
	height[0] = CorrelationFunctionLorentzians::height(0.3, hwhm[0]);*/

	//Eigen::VectorXd omega(2);
	//Eigen::VectorXd hwhm(2);
	//Eigen::VectorXd height(2);
	//omega[0] = 1; hwhm[0] = 0.1; height[0] = CorrelationFunctionLorentzians::height(0.64, 0.1);
	//omega[1] = 0.3; hwhm[1] = 0.05; height[1] = CorrelationFunctionLorentzians::height(0.025, 0.05);

//	Eigen::VectorXd omega(1);
//	Eigen::VectorXd hwhm(1);
//	Eigen::VectorXd height(1);
//	omega[0] = 1;
////	omega[1] = 0.46;
////	omega[2] = 0;
//	hwhm.fill(1E-8);
//	// scale = g^2
//	height.fill(CorrelationFunctionLorentzians::height(0.25, hwhm[0]));//CorrelationFunctionLorentzians::height(0.64, hwhm[0]);

	return boost::make_shared<CorrelationFunctionLorentzians>(omega, hwhm, height);
}

int main(int argc, char* argv[])
{	
	const time_t t0 = clock();
	CommandLineArgumentsReader clar(argc, argv);
	size_t nbrSites, initSite;
	double C;
	double dt;
	double tMax;
	size_t nbr_paths;
	double scale;
	double gamma;
	double omega0;
	try {
		clar.pop(nbrSites);
		clar.pop(initSite);		
		clar.pop(C);	
		clar.pop(scale);		
		clar.pop(gamma);
		clar.pop(omega0);
		clar.pop(dt);
		clar.pop(tMax);		
		clar.pop(nbr_paths);		
	} catch (...) {
		std::cerr << "Usage: " << argv[0] << " nbrSites initSite C scale gamma omega0 dt tMax nbr_paths" << std::endl;
		return -1;
	}
	const size_t nbrSteps = static_cast<size_t>(ceil(tMax/dt));

	/* RING */
	//const double J = HamiltonianProperties::absorption_band_shift_to_interation_strength_ring_nn(nbrSites, C);
	/* CHAIN */
	const double J = C/2;

	std::cout << "J == " << J << std::endl;

	dt /= std::abs(C);
	tMax /= std::abs(C);
	std::cout << "tMax: " << tMax << std::endl;
	
	std::string outf_name;
	{
		std::stringstream ss;
		//ss << "strunz_RingPRL_N" << nbrSites << "_J" << J << "_dt" << dt << "_tMax" << tMax << "_nbrpaths" << nbr_paths << ".txt"; 
		ss << "strunz_Chain_N" << nbrSites << "_scale" << scale << "_gamma" << gamma << "_omega" << omega0 << "_J" << J << "_dt" << dt << "_tMax" << tMax << "_nbrpaths" << nbr_paths << ".txt";
		outf_name = ss.str();
	}
	std::cout << "Output filename: " << outf_name << std::endl;

	try {
		const time_t t0 = clock();		
		const boost::shared_ptr<const CorrelationFunctionLorentzians> monomerAlpha = buildFMOAlpha(scale, gamma, omega0);
		//const boost::shared_ptr<const CorrelationFunction> monomerAlpha = boost::make_shared<CorrelationFunctionModes>();
		/*const size_t nbrModesPerSite = 2;
		const boost::shared_ptr<const CorrelationFunction> monomerAlpha = buildModeAlpha(nbrModesPerSite);*/
		//const boost::shared_ptr<const ExactHamiltonian> exactHamiltonian = buildSingleModeExactHamiltonian(nbrSites, J, nbrModesPerSite);

		Eigen::VectorXcd initial(nbrSites);
		initial.fill(0.0);
		initial[initSite] = 1.0;
		
		Eigen::MatrixXcd hamiltonian(nbrSites, nbrSites);
		//HamiltonianFactory::buildRing(std::vector<double>(nbrSites, 0), J, hamiltonian);
		HamiltonianFactory::buildChain(std::vector<double>(nbrSites, 0), J, 0.0, hamiltonian);
		
		typedef StrunzSimulatorReducedStochasticNonlinear simulator_type;
		simulator_type simulator(monomerAlpha, nbrSites, hamiltonian, dt, nbrSteps);
		
		std::vector<Eigen::VectorXcd> states;
		const time_t tsim_0 = clock();
		simulator_type::Workspace wksp(simulator.workspace());
		const time_t tsim_1 = clock();
		std::cout << "Simulator workspace: " << (tsim_1 - tsim_0) / static_cast<double>(CLOCKS_PER_SEC) << std::endl;

		std::vector<rql::math::Average<Eigen::MatrixXcd> > density_matrices(nbrSteps + 1, rql::math::Average<Eigen::MatrixXcd>(Eigen::MatrixXcd::Zero(nbrSites, nbrSites)));
		rql::math::Jagged2DArray<rql::math::Average<double> > occupation_numbers(nbrSteps + 1, nbrSites);
		const time_t ttraj_0 = clock();
		for (size_t i = 0; i < nbr_paths; ++i) {
			simulator.simulate(wksp, initial, states);
			for (size_t t = 0; t < nbrSteps + 1; ++t) {
				for (size_t s = 0; s < nbrSites; ++s) {
					occupation_numbers(t, s).update(pow(std::abs(states[t][s]), 2));
				}
				density_matrices[t].update(states[t] * states[t].adjoint());
			}
			std::cout << "Path no " << i << " completed" << std::endl;
		}
		const time_t ttraj_1 = clock();
		std::cout << "Trajectory: " << (ttraj_1 - ttraj_0) / static_cast<double>(CLOCKS_PER_SEC) << std::endl;

		const time_t t1 = clock();
		std::cout << "Total: " << (t1 - t0) / static_cast<double>(CLOCKS_PER_SEC) << std::endl;
		
		std::ofstream outf(outf_name, std::fstream::out);

		// for H0 * rho product
		Eigen::MatrixXcd rho_times_H0(nbrSites, nbrSites);
		
		for (size_t tIdx = 0; tIdx < nbrSteps + 1; ++tIdx) {
			const double t = tIdx * dt;
			const double t_in_PRL_units = t * std::abs(C);
			outf << t;//t_in_PRL_units;
			for (size_t s = 0; s < nbrSites; ++s) {
				outf << "\t" << occupation_numbers(tIdx, s).value();
			}
			const Eigen::MatrixXcd& rho = density_matrices[tIdx].value();
			
			// Average system energy
			rho_times_H0.noalias() = hamiltonian * rho;
			const double energy = rho_times_H0.trace().real();
			
			// calculate coherence
			try {
				RhoAnalyzer rho_analyzer(rho);
				const double von_neumann_entropy = rho_analyzer.vonNeumannEntropy();
				const double coherence = rho_analyzer.coherence();
				outf << std::setprecision(16) << "\t" << von_neumann_entropy << "\t" << coherence;
			} catch (std::exception& e) {
				std::stringstream ss;
				ss << "Bad density matrix at t == " << t_in_PRL_units << ": " << e.what();
				throw std::runtime_error(ss.str().c_str());
			}
			
			// print system energy
			outf << std::setprecision(16) << "\t" << energy;
			
			outf << "\n";
		}
		outf.close();

		return 0;
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return -1;
	}

	const time_t t1 = clock();
	std::cout << "Running time [s]: " << (t1 - t0) / static_cast<double>(CLOCKS_PER_SEC) << std::endl;
}
