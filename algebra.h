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
			cholmod_dense *b, cholmod_common *cm, size_t &peak_mem, size_t &CK_mem);
	static void factor_to_triplet(cholmod_factor *L, float *&L_h, size_t &L_h_nz);
	static void solve_CK_for_back_sub(int &my_id, int &L_nz, float *L, float *b, int &base_nz, int &base_n);
	static void LU_decomposition(int n, UF_long * Ap, UF_long * Ai, 
			double * Ax,
		void ** Numeric);
	static void CK_decomp(Matrix &A, cholmod_factor *&L, 
			cholmod_common *cm, size_t &peak_mem, size_t &CK_mem);
};

#endif
