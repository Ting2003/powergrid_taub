// ----------------------------------------------------------------//
// Filename : algebra.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// linear algebra library
// ----------------------------------------------------------------//
// - Zigang Xiao - Wed Jan 26 18:15:02 CST 2011
//   * file created

#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

#include "global.h"
#include "vec.h"
#include "umfpack.h"
#include "cholmod.h"
class Vec;
class Algebra{
public:
	static void solve(const Matrix & A, const Vec & b, Vec & x);
	static void solve_CK(Matrix & A, cholmod_dense *&x, 
			cholmod_dense *b, cholmod_common *cm);
	static void LU_decomposition(int n, UF_long * Ap, UF_long * Ai, 
			double * Ax,
		void ** Numeric);
	static void CK_decomp(Matrix &A, cholmod_factor *&L, 
			cholmod_common *cm);
};

#endif
