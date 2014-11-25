#ifndef __EXACT_HAMILTONIAN_HPP
#define __EXACT_HAMILTONIAN_HPP

#include <Eigen/Core>
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>
#include <vector>
#include "core.h"

class ExactHamiltonian
{
public:
	PROTEIN_CHAIN_API ExactHamiltonian(size_t nbrSites);
	PROTEIN_CHAIN_API virtual ~ExactHamiltonian() {}
	PROTEIN_CHAIN_API virtual size_t dimension(size_t nbrExcitedBathLevels) const = 0;
	PROTEIN_CHAIN_API virtual void generateMatrix(Eigen::MatrixXcd& matrix, size_t nbrExcitedBathLevels) const = 0;
	PROTEIN_CHAIN_API virtual void generateMatrix(Eigen::SparseMatrix<std::complex<double> >& matrix, size_t nbrExcitedBathLevels) const = 0;
	PROTEIN_CHAIN_API virtual void convertState(const Eigen::VectorXcd& excitonState, Eigen::VectorXcd& fullState, size_t nbrExcitedBathLevels) const = 0;
	PROTEIN_CHAIN_API size_t nbrSites() const { return m_nbr_sites; }
	
	//! g: matrix of coefficients multiplying creation operators
	PROTEIN_CHAIN_API static void independentBaths(const Eigen::VectorXcd& siteGs, size_t nbrSites, Eigen::MatrixXcd& g, Eigen::VectorXd& omegas);
	
	//! g: matrix of coefficients multiplying creation operators
	PROTEIN_CHAIN_API static void independentBaths(const std::complex<double>& g0, size_t nbrModesPerSite, size_t nbrSites, Eigen::MatrixXcd& g, Eigen::VectorXd& omegas);
private:
	size_t m_nbr_sites;
};

//! Just one mode
class DiagonalizedChain: public ExactHamiltonian
{
public:
	//! @param[in] E Vector of energies
	//! @param[in] omega Oscillator frequency
	//! @param[in] g Coupling matrix for creation operators
	PROTEIN_CHAIN_API DiagonalizedChain(const Eigen::VectorXd& E, double omega, const Eigen::MatrixXcd& g);
	PROTEIN_CHAIN_API void generateMatrix(Eigen::MatrixXcd& matrix, size_t nbrExcitedBathLevels) const;
	PROTEIN_CHAIN_API void generateMatrix(Eigen::SparseMatrix<std::complex<double> >& matrix, size_t nbrExcitedBathLevels) const;
	PROTEIN_CHAIN_API size_t dimension(size_t nbrExcitedBathLevels) const { return nbrSites() * (nbrExcitedBathLevels + 1); }
	PROTEIN_CHAIN_API void convertState(const Eigen::VectorXcd& excitonState, Eigen::VectorXcd& fullState, size_t nbrExcitedBathLevels) const;
private:
	template <class M> void generateMatrixImpl(M& matrix, size_t nbrExcitedBathLevels) const;
	Eigen::VectorXd m_E;
	double m_omega;
	Eigen::MatrixXcd m_g;  
};

//! Hamiltonian with system-bath coupling "sum_m |m><m| sum_k (g_km b_k^dagger + g_km^* b_k)"
class DiagCoupling: public ExactHamiltonian
{
public:
	//! @param[in] nbrSites Number of chain sites
	//! @param[in] omega Oscillator frequencies (can be imaginary, but real part must be > 0)
	//! @param[in] g Bath-chain coupling matrix for creation operators (g.cols() == omega.size(), g.rows() == nbrSites)
	PROTEIN_CHAIN_API DiagCoupling(size_t nbrSites, const Eigen::VectorXcd& omega, const Eigen::MatrixXcd& g);
	PROTEIN_CHAIN_API void generateMatrix(Eigen::MatrixXcd& matrix, size_t nbrExcitedBathLevels) const;
	PROTEIN_CHAIN_API void generateMatrix(Eigen::SparseMatrix<std::complex<double> >& matrix, size_t nbrExcitedBathLevels) const;
	PROTEIN_CHAIN_API size_t dimension(size_t nbrExcitedBathLevels) const;
	PROTEIN_CHAIN_API void convertState(const Eigen::VectorXcd& excitonState, Eigen::VectorXcd& fullState, size_t nbrExcitedBathLevels) const;
	//! @param[in] nbrModes Number of bath modes
	//! @param[in] nbrLevels Number of excitation levels (incl. vacuum) for each bath mode
	//! @param[out] basis Generated basis, on return has size bathBasisDimension(nbrModes, nbrLevels)*nbrModes and contains ordered basis vectors
	PROTEIN_CHAIN_API static void generateBathBasis(size_t nbrModes, size_t nbrLevels, std::vector<unsigned int>& basis);
	//! @param[in] nbrModes Number of bath modes
	//! @param[in] nbrLevels Number of excitation levels (incl. vacuum) for each bath mode
	PROTEIN_CHAIN_API static size_t bathBasisDimension(size_t nbrModes, size_t nbrLevels);
private:
	virtual std::complex<double> h0(size_t m, size_t n) const = 0;
	template <class M> void generateMatrixImpl(M& matrix, size_t nbrExcitedBathLevels) const;
	Eigen::VectorXcd m_omega;
	Eigen::MatrixXcd m_g;	
};

//! Predefined electronic Hamiltonian
class PredefinedElectronic: public DiagCoupling
{
public:
	//! @param[in] Hel Electronic hamiltonian (Hel.rows() == Hel.cols() == number of sites)
	//! @param[in] omega Oscillator frequencies (can be imaginary, but real part must be > 0)
	//! @param[in] g Bath-chain coupling matrix for creation operators (g.cols() == omega.size(), g.rows() == number of sites)
	PROTEIN_CHAIN_API PredefinedElectronic(const Eigen::MatrixXcd& Hel, const Eigen::VectorXcd& omega, const Eigen::MatrixXcd& g);
private:
	std::complex<double> h0(size_t m, size_t n) const { return m_Hel(m, n); }
	Eigen::MatrixXcd m_Hel;
};

//! Chain with nearest-neighbour interaction
class NearestNeighbour: public DiagCoupling
{
public:
	//! @param[in] nbrSites Number of chain sites
	//! @param[in] omega Oscillator frequencies (can be imaginary, but real part must be > 0)
	//! @param[in] g Bath-chain coupling matrix for creation operators (g.cols() == omega.size(), g.rows() == nbrSites)
	//! @param[in] J Nearest-neighbour interaction constant (defaults to -1).
	PROTEIN_CHAIN_API NearestNeighbour(size_t nbrSites, const Eigen::VectorXcd& omega, const Eigen::MatrixXcd& g, double J = -1);
private:
	std::complex<double> h0(size_t m, size_t n) const;
	double m_J;
};

#endif // __EXACT_HAMILTONIAN_HPP
