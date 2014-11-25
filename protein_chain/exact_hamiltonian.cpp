#include "exact_hamiltonian.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

template <class M, class V> struct MatrixSetter
{
	inline static void setZero(M& matrix, size_t dim);
	inline static void set(M& matrix, size_t i, size_t j, const V& value);
};

template <> struct MatrixSetter<Eigen::MatrixXcd, std::complex<double> >
{
	inline static void zero(Eigen::MatrixXcd& matrix, size_t dim)
	{
		matrix.setZero(dim, dim);
	}

	inline static void set(Eigen::MatrixXcd& matrix, size_t i, size_t j, const std::complex<double>& value)
	{
		matrix(i, j) = value;
	}
};

template <class V> struct MatrixSetter<Eigen::SparseMatrix<V>, V>
{
	inline static void zero(Eigen::SparseMatrix<V>& matrix, size_t dim)
	{
		matrix.setZero();		
	}

	inline static void set(Eigen::SparseMatrix<V>& matrix, size_t i, size_t j, const V& value)
	{
		if (value != 0.0)
			matrix.coeffRef(i, j) = value;
	}
};

ExactHamiltonian::ExactHamiltonian(size_t nbrSites)
	: m_nbr_sites(nbrSites)
{
	assert(m_nbr_sites);
}

static void makeIndependent(Eigen::VectorXd& omegas, const size_t nbrSites)
{
	const size_t nbrModesPerSite = omegas.size();
	Eigen::VectorXd newOmegas(nbrSites * nbrModesPerSite);
	for (size_t i = 0; i < nbrSites; ++i) {
		for (size_t j = 0; j < nbrModesPerSite; ++j) {
			newOmegas[i * nbrModesPerSite + j] = omegas[j];
		}
	}
	omegas.swap(newOmegas);
}

void ExactHamiltonian::independentBaths(const Eigen::VectorXcd& siteGs, const size_t nbrSites, Eigen::MatrixXcd& g, Eigen::VectorXd& omegas)
{
	const int nbr_modes_per_site = siteGs.size();
	g.setZero(nbrSites, nbr_modes_per_site * nbrSites);
	for (size_t siteIdx = 0; siteIdx < nbrSites; ++siteIdx) {
		for (int modeIdx = 0; modeIdx < nbr_modes_per_site; ++modeIdx) {
			g(siteIdx, siteIdx * nbr_modes_per_site + modeIdx) = siteGs[modeIdx];
		}
	}
	makeIndependent(omegas, nbrSites);
}

void ExactHamiltonian::independentBaths(const std::complex<double>& g0, const size_t nbrModesPerSite, size_t nbrSites, Eigen::MatrixXcd& g, Eigen::VectorXd& omegas)
{
	g.setZero(nbrSites, nbrModesPerSite * nbrSites);
	for (size_t siteIdx = 0; siteIdx < nbrSites; ++siteIdx) {
		g.block(siteIdx, siteIdx * nbrModesPerSite, 1, nbrModesPerSite).fill(g0);
	}
	makeIndependent(omegas, nbrSites);
}

DiagonalizedChain::DiagonalizedChain(const Eigen::VectorXd& E, double omega, const Eigen::MatrixXcd& g)
	: ExactHamiltonian(E.size()), m_E(E), m_omega(omega), m_g(g)
{
	assert(g.rows() == static_cast<int>(nbrSites()));
	assert(g.cols() == static_cast<int>(nbrSites()));  
}

template <class M> void DiagonalizedChain::generateMatrixImpl(M& matrix, size_t nbrExcitedBathLevels) const
{
	const size_t dim = dimension(nbrExcitedBathLevels);
	const size_t nbl = nbrExcitedBathLevels + 1;
	MatrixSetter<M,std::complex<double> >::zero(matrix, dim);
	assert(matrix.rows() == dim);
	assert(matrix.cols() == dim);
	for (size_t i1 = 0; i1 < nbrSites(); ++i1) {
		const double Ei1 = m_E[i1];
		for (size_t j1 = 0; j1 < nbl; ++j1) {
			const double sqrt_j1 = sqrt(static_cast<double>(j1));
			const double sqrt_j1p1 = sqrt(static_cast<double>(j1 + 1));
			const size_t idx1 = i1*nbl+ j1;
			MatrixSetter<M,std::complex<double> >::set(matrix, idx1, idx1, Ei1 + m_omega*j1);
			for (size_t i2 = 0; i2 < nbrSites(); ++i2) {
				const std::complex<double> g12 = m_g(i1, i2);
				if (j1 > 0) {
					const size_t j2 = j1 - 1; // creation operator
					const size_t idx2 = i2*nbl + j2;
					MatrixSetter<M,std::complex<double> >::set(matrix, idx1, idx2, g12 * sqrt_j1);
				}
				if (j1 < nbl - 1) {
					const size_t j2 = j1 + 1; // annihilation operator
					const size_t idx2 = i2*nbl + j2;
					MatrixSetter<M,std::complex<double> >::set(matrix, idx1, idx2, conj(g12) * sqrt_j1p1);
				}	
			}
		}
	}
}

void DiagonalizedChain::generateMatrix(Eigen::MatrixXcd& m, size_t nbrExcitedBathLevels) const
{
	generateMatrixImpl(m, nbrExcitedBathLevels);
}

void DiagonalizedChain::generateMatrix(Eigen::SparseMatrix<std::complex<double> >& m, size_t nbrExcitedBathLevels) const
{
	generateMatrixImpl(m, nbrExcitedBathLevels);
}

void DiagonalizedChain::convertState(const Eigen::VectorXcd& excitonState, Eigen::VectorXcd& fullState, size_t nbrExcitedBathLevels) const
{
	assert( excitonState.size() == static_cast<int>(nbrSites()) );
	fullState.setZero(dimension(nbrExcitedBathLevels));
	for (int i = 0; i < excitonState.size(); ++i)
		fullState[i * (nbrExcitedBathLevels + 1)] = excitonState[i];
}

DiagCoupling::DiagCoupling(size_t nbrSites, const Eigen::VectorXcd& omega, const Eigen::MatrixXcd& g)
	: ExactHamiltonian(nbrSites), m_omega(omega), m_g(g)
{
	assert( m_omega.size() == m_g.cols() );
	assert( static_cast<size_t>(m_g.rows()) == nbrSites );
}

template <class M> void DiagCoupling::generateMatrixImpl(M& matrix, size_t nbrExcitedBathLevels) const
{
	std::vector<unsigned int> basis;	
	const size_t nbr_modes = m_omega.size();
	generateBathBasis(nbr_modes, nbrExcitedBathLevels + 1, basis);
	const size_t bath_basis_dim = basis.size() / m_omega.size();
	const size_t total_basis_dim = bath_basis_dim * nbrSites();
	MatrixSetter<M,std::complex<double> >::zero(matrix, total_basis_dim);
	assert(matrix.rows() == total_basis_dim);
	assert(matrix.cols() == total_basis_dim);
	for (size_t m = 0; m < nbrSites(); ++m) {
		for (size_t n = 0; n < nbrSites(); ++n) {
			const std::complex<double> Vmn = h0(m, n);
			for (size_t i1 = 0; i1 < bath_basis_dim; ++i1) {
				const size_t b_start = i1 * nbr_modes;
				std::complex<double> sum = Vmn;
				if (m == n) {
					for (size_t k = 0; k < nbr_modes; ++k) {
						sum += m_omega[k] * static_cast<double>(basis[b_start + k]);
					}
				}
				const size_t j1 = m * bath_basis_dim + i1;
				size_t j2 = n * bath_basis_dim + i1;
				MatrixSetter<M,std::complex<double> >::set(matrix, j1, j2, sum);
				if (m == n) {
					// bath-system coupling
					for (size_t k = 0; k < nbr_modes; ++k) {						
						const unsigned int exc_level = basis[b_start + k];
						if (exc_level > 0) {
							const double f = sqrt(static_cast<double>(exc_level));
							const size_t i2 = i1 - static_cast<size_t>(pow(static_cast<double>(nbrExcitedBathLevels + 1), static_cast<int>(k)));
							assert( i2 != i1 );
							assert( basis[i2 * nbr_modes + k] == exc_level - 1 ); // i2 should be the less-excited bath state
							const std::complex<double> mel = f * m_g(m, k);
							j2 = m * bath_basis_dim + i2; // index of the less-excited state
							MatrixSetter<M,std::complex<double> >::set(matrix, j1, j2, mel); // creation operator
							MatrixSetter<M,std::complex<double> >::set(matrix, j2, j1, conj(mel)); // annihilation operator
						}
					}
				}
			}
		}
	}
}

void DiagCoupling::generateMatrix(Eigen::MatrixXcd& matrix, size_t nbrExcitedBathLevels) const
{
	generateMatrixImpl(matrix, nbrExcitedBathLevels);
}

void DiagCoupling::generateMatrix(Eigen::SparseMatrix<std::complex<double> >& m, size_t nbrExcitedBathLevels) const
{
	generateMatrixImpl(m, nbrExcitedBathLevels);
}

void DiagCoupling::generateBathBasis(const size_t nbrModes, const size_t nbrLevels, std::vector<unsigned int>& basis)
{
	const size_t dim = bathBasisDimension(nbrModes, nbrLevels);
	basis.resize(dim * nbrModes);
	std::vector<unsigned int>::iterator bit = basis.begin();
	for (size_t i = 0; i < dim; ++i) {
		size_t i2 = i;
		for (size_t j = 0; j < nbrModes; ++j) {
			assert( bit != basis.end() );
			*bit = i2 % nbrLevels;
			i2 /= nbrLevels;
			++bit;
		}
	}
}

size_t DiagCoupling::bathBasisDimension(size_t nbrModes, size_t nbrLevels)
{
	return static_cast<size_t>(pow(static_cast<double>(nbrLevels), static_cast<int>(nbrModes)));
}

size_t DiagCoupling::dimension(size_t nbrExcitedBathLevels) const
{
	return nbrSites() * bathBasisDimension(m_omega.size(), nbrExcitedBathLevels + 1);
}

void DiagCoupling::convertState(const Eigen::VectorXcd& excitonState, Eigen::VectorXcd& fullState, size_t nbrExcitedBathLevels) const
{
	assert( excitonState.size() == static_cast<int>(nbrSites()) );
	const size_t bathDim = bathBasisDimension(m_omega.size(), nbrExcitedBathLevels + 1);
	fullState.setZero(dimension(nbrExcitedBathLevels));
	for (int i = 0; i < excitonState.size(); ++i)
		fullState[i * bathDim] = excitonState[i];
}

PredefinedElectronic::PredefinedElectronic(const Eigen::MatrixXcd& Hel, const Eigen::VectorXcd& omega, const Eigen::MatrixXcd& g)
	: DiagCoupling(Hel.rows(), omega, g), m_Hel(Hel)
{
	if (Hel.rows() != Hel.cols())
		throw std::domain_error("Electronic hamiltonian must be a square matrix");
}

NearestNeighbour::NearestNeighbour(size_t nbrSites, const Eigen::VectorXcd& omega, const Eigen::MatrixXcd& g, double J)
	: DiagCoupling(nbrSites, omega, g), m_J(J)
{
	assert(g.rows() == static_cast<int>(nbrSites));
	assert(g.cols() == omega.size());
}

std::complex<double> NearestNeighbour::h0(size_t m, size_t n) const
{
	if (m + 1 == n || n + 1 == m)
		return m_J;
	else
		return 0;
}
