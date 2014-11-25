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
#include <protein_chain/correlation_function_drude.h>
#include <protein_chain/correlation_function_lorentzians.h>
#include <protein_chain/exact_hamiltonian.h>
#include <protein_chain/evolver_taylor_expansion.h>
#include <protein_chain/hamiltonian_properties.h>
#include <protein_chain/hamiltonian_factory.h>
#include <protein_chain/rho_analyzer.h>
#include <qsd/strunz_simulator_reduced_stochastic_linear.h>
#include <qsd/strunz_simulator_reduced_stochastic_nonlinear.h>
#include <sstream>
#include <fstream>
#include <boost/math/special_functions/fpclassify.hpp>

static const size_t nbrSites = 3;

boost::shared_ptr<const CorrelationFunctionLorentzians> Drude11Lor(double lambda, double omega_cutoff, double kBT, size_t pfdOrder)
{
	static const size_t nbr_peaks = 22;
	static const double lambda_Eisfeld = 35;
	static const double omega_cutoff_Eisfeld = 166;
	Eigen::VectorXd omega(nbr_peaks);
	Eigen::VectorXd hwhm(nbr_peaks);
	Eigen::VectorXd scale(nbr_peaks);
	Eigen::VectorXd height(nbr_peaks);
	/*omega[0] = 20177.28213; omega[1] = 1797.669703; omega[2] = 849.5608213; omega[3] = 301.4274205; omega[4] = 177.0297857; omega[5] = 495.6284577; omega[6] = 141.6411068;
	omega[7] = 239.9693839; omega[8] = 40.1414793; omega[9] = 358.2380846; omega[10] = 342.4112912;

	scale[0] = 17.32642035; scale[1] = 11154.7129; scale[2] = 64860.57973; scale[3] = -1067947.525; scale[4] = 623149.5892; scale[5] = 270997.5806; scale[6] = 21412.67442;
	scale[7] = 5.501245116; scale[8] = 4340.402898; scale[9] = -28413.01945; scale[10] = 94562.0744;

	hwhm[0] = 38.1001007; hwhm[1] = 944.3095657; hwhm[2] = 736.153262; hwhm[3] = 616.695154; hwhm[4] = 661.5906125; hwhm[5] = 493.879183; hwhm[6] = 179.321305;
	hwhm[7] = 7.70118E-05; hwhm[8] = 112.3387044; hwhm[9] = 190.8494941; hwhm[10] = 241.5086996;*/

	scale[0] = 17.32642665;	scale[1] = 11807.35028;	scale[2] = 204.5171125;	scale[3] = 19373.10798;	scale[4] = 200.0083777;	scale[5] = 4757.388623;	scale[6] = 15753.17797;
	scale[7] = 5.489001814;	scale[8] = 5652.634563;	scale[9] = 233.9980317;	scale[10] = 1802.796124;

	omega[0] = 20301.64124; omega[1] = 1528.093874; omega[2] = 734.0829801; omega[3] = 225.7808956; omega[4] = 19.85946792; omega[5] = 416.0557036; omega[6] = 80.36144649;
	omega[7] = 134.5547866; omega[8] = 40.88679341; omega[9] = 304.1627838; omega[10] = 366.3911778;

	hwhm[0] = 38.3806493; hwhm[1] = 2329.334014; hwhm[2] = 1081.754556; hwhm[3] = 910.2629255; hwhm[4] = 1402.675083; hwhm[5] = 1207.678113; hwhm[6] = 280.4383199;
	hwhm[7] = 7.73958E-05; hwhm[8] = 165.0275369; hwhm[9] = 229.7312668; hwhm[10] = 392.8450466;

	for (size_t m = 0; m < nbr_peaks/2; ++m) {
		scale[m + nbr_peaks/2] = -scale[m];
		omega[m + nbr_peaks/2] = -omega[m];
		hwhm[m + nbr_peaks/2] = hwhm[m];
	}

	for (size_t m = 0; m < nbr_peaks; ++m) 
		height[m] = CorrelationFunctionLorentzians::height(scale[m], hwhm[m]);

	const double scale_omega = omega_cutoff / omega_cutoff_Eisfeld;
	const double scale_lambda = lambda / lambda_Eisfeld;
	omega *= scale_omega;
	hwhm *= scale_omega;
	height *= scale_lambda;
	scale *= (scale_lambda * scale_omega);
	return boost::shared_ptr<const CorrelationFunctionLorentzians>(new CorrelationFunctionLorentzians(omega, hwhm, height, kBT, pfdOrder));
}

void build_subunit_hamiltonian(Eigen::MatrixXcd& H, const double energy_shift)
{
	H.resize(nbrSites, nbrSites);
	/*// Eisfeld's Hamiltonian from NJP 2011
	H(0,0) = 410; 
	H(1,1) = 530;
	H(2,2) = 210;
	H(3,3) = 320;
	H(4,4) = 480;
	H(5,5) = 630;
	H(6,6) = 440;
	H(0,1) = -87.7;
	H(0,2) = 5.5;
	H(0,3) = -5.9;
	H(0,4) = 6.7;
	H(0,5) = -13.7;
	H(0,6) = -9.9;
	H(1,2) = 30.8;
	H(1,3) = 8.2;
	H(1,4) = 0.7;
	H(1,5) = 11.8;
	H(1,6) = 4.3;
	H(2,3) = -53.5;
	H(2,4) = -2.2;
	H(2,5) = -9.6;
	H(2,6) = 6.0;
	H(3,4) = -70.7;
	H(3,5) = -17.0;
	H(3,6) = -63.3;
	H(4,5) = 81.1;
	H(4,6) = -1.3;
	H(5,6) = 39.7;*/
	
	H(0,0) = 12653; // Edward's Hamiltonian
	H(0,1) = -45.59;
	H(0,2) = 3.42;
	H(1,2) = 273;
	H(1,1) = 12325;
	H(2,2) = 12025;

/*	H(0,0) = 100; // Ishazaki&Flemming's dimer Hamiltonian
	H(0,1) = 100;
	H(1,1) = 0;*/
	

	for (size_t m = 0; m < nbrSites; ++m) {
		for (size_t n = 0; n < m; ++n) {
			H(m, n) = H(n, m);
		}
	}
	for (size_t m = 0; m < nbrSites; ++m) {
		H.diagonal()[m] += energy_shift;
	}
}

boost::shared_ptr<const CorrelationFunctionDecomposable> build_correlation_function(double lambda, double omega_c, double kBT, size_t pfdOrder, bool use_lorentzian)
{
	if (use_lorentzian) {
		return CorrelationFunctionDrude::approximate_by_lorentzian(lambda, omega_c, kBT);
	} else {
		return boost::shared_ptr<const CorrelationFunctionDecomposable>(new CorrelationFunctionDrude(lambda, omega_c, kBT, pfdOrder));
	}
}

void build_subunit_alphas_Edward(std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >& alphas, double kBT, size_t pfdOrder, bool use_lorentzian)
{
	alphas.resize(nbrSites);
	// Edward's bath
	alphas[0] = Drude11Lor(3, 106, kBT, pfdOrder);
		//build_correlation_function(3E-1, 106, kBT, pfdOrder, use_lorentzian); // boost::shared_ptr<const CorrelationFunctionDecomposable>(new CorrelationFunctionDrude(30, 106, kBT, pfdOrder));
	alphas[1] = alphas[2] = Drude11Lor(3, 106, kBT, pfdOrder);
		//build_correlation_function(3E-1, 106, kBT, pfdOrder, use_lorentzian);// boost::shared_ptr<const CorrelationFunctionDecomposable>(new CorrelationFunctionDrude(300, 106, kBT, pfdOrder));	
}

void build_subunit_alphas_Ishizaki(std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> >& alphas, double kBT, double lambda, size_t pfdOrder, bool use_lorentzian)
{
	alphas.resize(nbrSites);
	alphas[0] = alphas[1] = alphas[2] = build_correlation_function(lambda, 53.08, kBT, pfdOrder, use_lorentzian); // boost::shared_ptr<const CorrelationFunctionDecomposable>(new CorrelationFunctionDrude(30, 106, kBT, pfdOrder));	
}

void save_corr_func_to_file(const boost::shared_ptr<const CorrelationFunctionDecomposable> alpha, const double dt, const size_t nbrSteps)
{
	std::ofstream outf("corr_func.txt", std::fstream::out);
	for (size_t i = 0; i < nbrSteps; ++i) {
		const double t = i*dt;
		outf << t;
		const std::complex<double> a = (*alpha)(t);
		outf << "\t" << a.real();
		outf << "\t" << a.imag();
		outf << "\n";
	}
	outf << "*********************************\n";
	outf << "*********************************\n";
	outf << "*********************************\n";
	for (size_t i = 0; i < alpha->nbr_exponents(); ++i) {
		outf << alpha->scale(i).real() << "\t" << alpha->scale(i).imag() << "\t" << alpha->exponent(i).real() << "\t" << alpha->exponent(i).imag();
		outf << "\n";
	}
	outf.close();
}

int main_off(int argc, char* argv[])
{	
	CommandLineArgumentsReader clar(argc, argv);
	double dt;
	double tMax;
	size_t nbr_paths;
	double kBT;
	size_t pfdOrder;
	bool use_lorentzian = false;
	double energy_shift = 0;
	try {
		clar.pop(dt);
		clar.pop(tMax);
		clar.pop(nbr_paths);
		clar.pop(kBT);
		clar.pop(pfdOrder);
		if (clar.hasNext()) {
			clar.pop(energy_shift);
		}
		if (clar.hasNext()) {
			int use_lorentzian_code;
			clar.pop(use_lorentzian_code);
			use_lorentzian = use_lorentzian_code != 0;
		}		
	} catch (...) {
		std::cerr << "Usage: " << argv[0] << " dt tMax nbr_paths kBT pfdOrder [energy_shift] [use_lorentzian]" << std::endl;
		return -1;
	}
	const size_t nbrSteps = static_cast<size_t>(ceil(tMax/dt));
	std::cout << "nbrSteps == " << nbrSteps << std::endl;

	try {
		const time_t t0 = clock();		

		Eigen::MatrixXcd hamiltonian(nbrSites, nbrSites);
		build_subunit_hamiltonian(hamiltonian, energy_shift);
		std::cout << "H0 == " << hamiltonian << std::endl;

		std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> > alphas(nbrSites); //(**)
		//build_subunit_alphas_Edward(alphas, kBT, pfdOrder, use_lorentzian); //(**)

		build_subunit_alphas_Ishizaki(alphas, kBT, 3, pfdOrder, use_lorentzian); //(**)

		// Eisfeld's NJP bath (**)
		// T = 77K ==> kBT = 53.51485374
		//const boost::shared_ptr<const CorrelationFunctionDecomposable> alpha = boost::shared_ptr<const CorrelationFunctionDecomposable>(new CorrelationFunctionDrude(35, 166, kBT, pfdOrder));
		//const boost::shared_ptr<const CorrelationFunctionDecomposable> alpha = Drude11Lor(35, 166, kBT, pfdOrder);
		//save_corr_func_to_file(alpha, dt, nbrSteps);
		
		const time_t tcorr = clock();
		std::cout << "Build Drude correlation functions: " << (tcorr - t0) / static_cast<double>(CLOCKS_PER_SEC) << std::endl;

		/*for (size_t i = 0; i < nbrSteps; ++i) {
			for (size_t k = 0; k < nbrSites; ++k) {
				std::cout << (*alphas[k])(i*dt) << "\t";
			}
			std::cout << "\n";
		}*/

		Eigen::VectorXcd initial(nbrSites);
		initial.fill(0.0);
		initial[0] = 1.0;
		
		
		const time_t tsimulator_0 = clock();
		typedef StrunzSimulatorReducedStochasticNonlinear simulator_type;
		simulator_type simulator(alphas, hamiltonian, dt, nbrSteps); //different baths (**)
		//simulator_type simulator(alpha, nbrSites, hamiltonian, dt, nbrSteps); // same baths (**)
		//simulator_type simulator(alphas[1], nbrSites, hamiltonian, dt, nbrSteps);
		const time_t tsimulator_1 = clock();
		std::cout << "Simulator: " << (tsimulator_1 - tsimulator_0) / static_cast<double>(CLOCKS_PER_SEC) << std::endl;
		
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
					const double p = pow(std::abs(states[t][s]), 2);
					if (!boost::math::isfinite(p)) {
						std::cerr << "Bad occupation probability at path " << i << ", time " << (t*dt) << ", site " << s << "\n";
						return -1;
					}
					occupation_numbers(t, s).update(p);
				}
				density_matrices[t].update(states[t] * states[t].adjoint());
			}
			std::cout << "Path no " << i << " completed" << std::endl;
		}
		const time_t ttraj_1 = clock();
		std::cout << "Trajectory: " << (ttraj_1 - ttraj_0) / static_cast<double>(CLOCKS_PER_SEC) << std::endl;

		

		const time_t t1 = clock();
		std::cout << "Total: " << (t1 - t0) / static_cast<double>(CLOCKS_PER_SEC) << std::endl;

		std::stringstream ss;
		ss << "strunz_Hsub3_11Lor_";
		if (use_lorentzian) {
			ss << "L";
		}
		ss << "_dt" << (dt) << "_tMax" << (tMax) << "_nbrpaths" << nbr_paths << "_kBT" << (kBT) << "_PFDorder" << pfdOrder << "_enShft" << energy_shift << ".txt"; 
		std::ofstream outf(ss.str(), std::fstream::out);

		for (size_t tIdx = 0; tIdx < nbrSteps + 1; ++tIdx) {
			const double t = tIdx * dt * 5.31;
			outf << t;
			for (size_t s = 0; s < nbrSites; ++s) {
				outf << "\t" << occupation_numbers(tIdx, s).value();
			}
			const Eigen::MatrixXcd& rho = density_matrices[tIdx].value();
			// calculate coherence
			try {
				RhoAnalyzer rho_analyzer(rho);
				const double von_neumann_entropy = rho_analyzer.vonNeumannEntropy();
				const double coherence = rho_analyzer.coherence();
				outf << std::setprecision(16) << "\t" << von_neumann_entropy << "\t" << coherence;
			} catch (std::exception& e) {
				std::stringstream ss;
				ss << "Bad density matrix at t == " << t << ": " << e.what();
				throw std::runtime_error(ss.str().c_str());
			}
			outf << "\n";
		}
		outf.close();

		return 0;
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return -1;
	}
}
