#include "diagonalizer_realhermitian_test.h"
#include "diagonalizer_realhermitian.h"
#include <cassert>
#include <vector>

using namespace Eigen;

void DiagonalizerRealHermitianTest::test(DiagonalizerRealHermitian& diagonalizer)
{
	assert( diagonalizer.dim() == dim() );
	MatrixXd m(dim(),dim());
	m(0,0) = -0.760160219372685;
	m(0,1) = -1.798287145164638;
	m(0,2) = -0.329360458723475;
	m(0,3) = -3.209614858617021;
	m(0,4) = 0.686778468359319;
	m(1,1) = -1.183810411330305;
	m(1,2) = -1.213699081438693;
	m(1,3) = -0.311828847760508;
	m(1,4) = -0.637383186114842;
	m(2,2) = -1.251804530839432;
	m(2,3) = 1.796771426918028;
	m(2,4) = -0.114860541786266;
	m(3,3) = 2.500445933634741;
	m(3,4) = 1.574155463616271;
	m(4,4) =  1.142138002541798;
	for (unsigned int i = 0; i < dim(); ++i) {
		for (unsigned int j = 0; j < i; ++j) {
			m(i, j) = m(j, i);
		}
	}
	VectorXd eigvals(dim());
	MatrixXd eigvecs(dim(),dim());
	diagonalizer.diagonalize(m, eigvecs, eigvals);
	const double tol = 1E-14;
	EXPECT_NEAR(-3.974999037407877, eigvals[0], tol);
	EXPECT_NEAR(-2.862643083238481, eigvals[1], tol);
	EXPECT_NEAR(0.136582399565815, eigvals[2], tol);
	EXPECT_NEAR(1.966701997107453, eigvals[3], tol);
	EXPECT_NEAR(5.181166498607206, eigvals[4], tol);
	MatrixXd v(5,5);
	v(0,0) = -0.738920670305702;
	v(0,1) = 0.087881155529801;
	v(0,2) = -0.183957838411956;
	v(0,3) = 0.469707488198307;
	v(0,4) = 0.437958370983997;
	v(1,0) = -0.495205103938613;
	v(1,1) = -0.482868722060628;
	v(1,2) = 0.557035775161652;
	v(1,3) = -0.459545512981308;
	v(1,4) = -0.011780024511297;
	v(2,0) = -0.022948705837562;
	v(2,1) = -0.756234437895661;
	v(2,2) = -0.596341638497358;
	v(2,3) = 0.102843372387863;
	v(2,4) = -0.247755368844335;
	v(3,0) = -0.424462788964514;
	v(3,1) = 0.351811985216194;
	v(3,2) = -0.083631850511528;
	v(3,3) = 0.007636639073194;
	v(3,4) = -0.830064493377708;
	v(4,0) = 0.167549532118020;
	v(4,1) = -0.251897800222942;
	v(4,2) = 0.541528413521545;
	v(4,3) = 0.746697877008326;
	v(4,4) = -0.240133109649572;
	for (unsigned int i = 0; i < dim(); ++i) {
		for (unsigned int j = 0; j <= i; ++j) {
			const double sp = eigvecs.col(i).dot(eigvecs.col(j));
			
		}
	}
	for (unsigned int i = 0; i < dim(); ++i) {
		for (unsigned int j = 0; j < dim(); ++j) {
			const double sp = eigvecs.col(i).dot(v.col(j));
			if (i == j) EXPECT_NEAR(1, std::abs(sp), tol) << i;
			else EXPECT_NEAR(0, sp, tol) << i << " " << j;
		}
	}
}
