#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <pseudomode/pseudomode_simulator.h>
#include <protein_chain/absorption_spectrum_calculator_time_dependent.h>
#include <protein_chain/command_line_arguments_reader.h>
#include <protein_chain/correlation_function_lorentzians.h>
#include <protein_chain/hamiltonian_factory.h>
#include <protein_chain/reductor.h>
#include <protein_chain/rho_analyzer.h>

struct SimulatorAndHamiltonian {
	std::shared_ptr<PseudomodeSimulator> simulator;
	std::shared_ptr<Eigen::MatrixXcd> hamiltonian;
};

// Gamma: scale
// gamma: HWHM
SimulatorAndHamiltonian buildSimulator(size_t nbrSites, double J, double dt, size_t nbrExcitBathLevels, double omega, double Gamma, double gamma)
{
	//NARROW PEAK
	/*PseudomodeSimulator::mode_param_triple parTrip;
	parTrip[PseudomodeSimulator::SMALL_GAMMA_IDX] = 0.125; // hwhm
	parTrip[PseudomodeSimulator::LARGE_GAMMA_IDX] = CorrelationFunctionLorentzians::scale(1, parTrip[PseudomodeSimulator::SMALL_GAMMA_IDX]); // convert height to scale
	parTrip[PseudomodeSimulator::OMEGA_IDX] = 1;	 // omega0 (centre of the peak)
	PseudomodeSimulator::mode_param_array bath_params(nbrSites, PseudomodeSimulator::mode_param_row(1, parTrip));*/

	//WIDE PEAK
	/*PseudomodeSimulator::mode_param_triple parTrip;
	parTrip[PseudomodeSimulator::SMALL_GAMMA_IDX] = 0.5; // hwhm
	parTrip[PseudomodeSimulator::LARGE_GAMMA_IDX] = CorrelationFunctionLorentzians::scale(0.5, parTrip[PseudomodeSimulator::SMALL_GAMMA_IDX]); // convert height to scale
	parTrip[PseudomodeSimulator::OMEGA_IDX] = 1;	 // omega0 (centre of the peak)
	PseudomodeSimulator::mode_param_array bath_params(nbrSites, PseudomodeSimulator::mode_param_row(1, parTrip));*/

	////NMSQ
	//PseudomodeSimulator::mode_param_triple parTrip;
	//parTrip[PseudomodeSimulator::SMALL_GAMMA_IDX] = 0.5; // hwhm
	//parTrip[PseudomodeSimulator::LARGE_GAMMA_IDX] = 0.64; // scale parameter
	//parTrip[PseudomodeSimulator::OMEGA_IDX] = 1;	 // omega0 (centre of the peak)
	//PseudomodeSimulator::mode_param_array bath_params(nbrSites, PseudomodeSimulator::mode_param_row(1, parTrip));

	//test case: weak narrow
	PseudomodeSimulator::mode_param_triple parTrip;
	parTrip[PseudomodeSimulator::SMALL_GAMMA_IDX] = gamma; // hwhm
	parTrip[PseudomodeSimulator::LARGE_GAMMA_IDX] = Gamma; // scale parameter
	parTrip[PseudomodeSimulator::OMEGA_IDX] = omega;	 // omega0 (centre of the peak)
	PseudomodeSimulator::mode_param_array bath_params(nbrSites, PseudomodeSimulator::mode_param_row(1, parTrip));
	
	/*PseudomodeSimulator::mode_param_triple parTrip;
	parTrip[PseudomodeSimulator::SMALL_GAMMA_IDX] = 0.25; // hwhm
	parTrip[PseudomodeSimulator::LARGE_GAMMA_IDX] = CorrelationFunctionLorentzians::scale(0.5, parTrip[PseudomodeSimulator::SMALL_GAMMA_IDX]); // convert height to scale
	parTrip[PseudomodeSimulator::OMEGA_IDX] = 1; // omega0 (centre of the peak)
	PseudomodeSimulator::mode_param_array bath_params(nbrSites, PseudomodeSimulator::mode_param_row(1, parTrip));*/


	// two peaks (example)
	//PseudomodeSimulator::mode_param_row peaks(2);
	//peaks[0][PseudomodeSimulator::SMALL_GAMMA_IDX] = 0.1;
	//peaks[0][PseudomodeSimulator::LARGE_GAMMA_IDX] = 0.64;
	//peaks[0][PseudomodeSimulator::OMEGA_IDX] = 1;
	//peaks[1][PseudomodeSimulator::SMALL_GAMMA_IDX] = 0.05;
	//peaks[1][PseudomodeSimulator::LARGE_GAMMA_IDX] = 0.025;
	//peaks[1][PseudomodeSimulator::OMEGA_IDX] = 0.3;
	//PseudomodeSimulator::mode_param_array bath_params(nbrSites, peaks);
	
	SimulatorAndHamiltonian result;
	result.hamiltonian = std::make_shared<Eigen::MatrixXcd>(nbrSites, nbrSites);
	HamiltonianFactory::buildChain(std::vector<double>(nbrSites, 0.0), J, 0.0, *(result.hamiltonian));
	result.simulator = std::make_shared<PseudomodeSimulator>(*(result.hamiltonian), bath_params, dt, nbrExcitBathLevels);
	return result;
}

void print_scalar_products(const std::vector<std::complex<double> >& scal_prods, const std::vector<double>& scal_prod_times)
{
	if (scal_prods.empty())
		return;
	const double tmax = scal_prod_times.back();
	const double output_time_interval = tmax / 201;
	double prev_time = 0;
	std::cout << 0 << "\t" << scal_prods[0].real() << "\t" << scal_prods[0].imag() << "\n";
	for (size_t i = 0; i < scal_prods.size(); ++i) {
		if (scal_prod_times[i] - prev_time > output_time_interval) {
			std::cout << scal_prod_times[i] << "\t" << scal_prods[i].real() << "\t" << scal_prods[i].imag() << "\n";
			prev_time = scal_prod_times[i];
		}
	}
	if (prev_time < tmax)
		std::cout << tmax << "\t" << scal_prods.back().real() << "\t" << scal_prods.back().imag() << "\n";
}

void print_density_matrix_info(std::vector<Eigen::MatrixXcd>& rhos, const std::vector<double>& scal_prod_times, const size_t nbrSites)
{
	if (rhos.empty())
		return;
	
}

int main(int argc, char* argv[])
{	
	const time_t t0 = clock();
	CommandLineArgumentsReader clar(argc, argv);
	size_t nbrSites;
	double C;
	double dt;
	double tmax;	
	double w0;
	double w1;
	size_t nbrExcitedBathLevels;
	double omega;
	double scale;
	double gamma;
	try {
		clar.pop(nbrSites);
		clar.pop(C);		
		//clar.pop(w0); // irrelevant for transport
		//clar.pop(w1); // irrelevant for transport
		clar.pop(dt);
		clar.pop(tmax);
		clar.pop(nbrExcitedBathLevels); // min. 20, but recommend 40
		clar.pop(omega);
		clar.pop(scale);
		clar.pop(gamma);
	} catch (...) {
		std::cerr << "Usage: " << argv[0] << " nbrSites C dt tmax nbrExcitedBathLevels omega scale gamma" << std::endl;
		return -1;
	}

	/* CHAIN */
	const double J = C/2;

	dt /= std::abs(C);
	tmax /= std::abs(C);

	const size_t nbrSteps = static_cast<size_t>(ceil(tmax/dt));	

	try {
		std::string outf_name;
		{
			std::stringstream ss;
			ss << "pseudoMode_transport2_omega" << omega << "_scale" << scale << "_gamma" << gamma << "_N" << nbrSites << "_dt" << dt << "_tMax" << tmax << "_m" << nbrExcitedBathLevels << "_J" << J << ".txt";
			outf_name = ss.str();
		}		
		std::cout << "Will save result in: " << outf_name << std::endl;
		
		const SimulatorAndHamiltonian sah = buildSimulator(nbrSites, J, dt, nbrExcitedBathLevels, omega, scale, gamma);
		const std::shared_ptr<PseudomodeSimulator> pmsim = sah.simulator;
		std::cout << "Built a simulator" << std::endl;

		Eigen::VectorXcd initialEl(nbrSites);
		//initialEl.fill(1/sqrt(static_cast<double>(nbrSites)));
		initialEl.setZero();
		initialEl[0] = 1;
		Eigen::VectorXcd initial(pmsim->dim());
		pmsim->convertState(initialEl, initial);
		std::cout << "Built the initial state" << std::endl;
		if (std::abs(initial.norm() - 1.0) > 1E-10) {
			throw std::runtime_error("Normalization error");
		}

		std::vector<Eigen::MatrixXcd> rhos;
		std::vector<std::complex<double> > scal_prods;
		pmsim->simulate(initial, nbrSteps, rhos, scal_prods);
		std::cout << "Carried out the simulation" << std::endl;
		assert(scal_prods.size() == nbrSteps + 1);
		assert(scal_prods[0] == 1.0);

		std::vector<double> scal_prod_times(scal_prods.size());
		for (size_t i = 0; i < scal_prods.size(); ++i)
			scal_prod_times[i] = i * dt;		
		
		std::ofstream outf(outf_name, std::fstream::out);
		Eigen::MatrixXcd rho_times_Hel(nbrSites, nbrSites);
		for (size_t i = 0; i < rhos.size(); ++i) {
			outf << scal_prod_times[i];
			Eigen::MatrixXcd& rho = rhos[i];
			const double orig_trace = rho.trace().real();
			rho /= rho.trace();
			for (size_t n = 0; n < nbrSites; ++n) {
				const auto v =rho(n, n); 
				if (v.real() < 0 || v.real() > 1) {
					std::cerr << "Occupation number outside allowed range at t == " << scal_prod_times[i] << std::endl;
					return -1;
				}
				outf << "\t" << v.real() << "\t" << v.imag();
			}
			const double nonherm = (rho - rho.adjoint()).norm();
			RhoAnalyzer rho_analyzer(rho, true);
			outf << "\t" << rho_analyzer.coherence();
			outf << "\t" << rho_analyzer.vonNeumannEntropy();
			
			// calculate average system energy
			rho_times_Hel.noalias() = rho * (*sah.hamiltonian);
			const double energy = rho_times_Hel.trace().real();
			outf << "\t" << energy;
			outf << "\t" << orig_trace;
			outf << "\t" << rho_analyzer.original_eigenvalues().minCoeff();
			outf << "\t" << rho_analyzer.original_eigenvalues().maxCoeff();
			outf << "\t" << nonherm;
			outf << "\n";
		}	
		outf.close();
		const time_t t1 = clock();		
		std::cout << "Running time [s]: " << (t1 - t0) / static_cast<double>(CLOCKS_PER_SEC) << std::endl;
		return 0;

		/*print_scalar_products(scal_prods, scal_prod_times);	
		return 0;*/

		/*
		const size_t nbrOmegas = 256;
		const double dw = (w1 - w0) / (nbrOmegas - 1);
		std::vector<double> sigma;
		AbsorptionSpectrumCalculatorTimeDependent abs_calc(w0, w1, dw);
		abs_calc.calculateAbsorptionSpectrum(scal_prod_times, scal_prods, sigma);

		for (size_t i = 0; i < nbrOmegas; ++i) {
			const double w = w0 + i*dw;
			std::cout << w << "\t" << sigma[i];
			//std::cout << "\t" << sigmaExact[i];
			std::cout << "\n";
		}
		*/
	

		return 0;
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return -1;
	}	
}
