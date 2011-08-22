// ----------------------------------------------------------------//
// Filename : block.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of block class
// ----------------------------------------------------------------//

#ifndef __BLOCK_H__
#define __BLOCK_H__
#include <fstream>
#include "triplet.h"
#include "global.h"
#include "vec.h"
#include "net.h"
#include "util.h"
#include "umfpack.h"
#include "cholmod.h"
#include "mpi_class.h"

using namespace std;

class Block{
	typedef vector<Net *> NetPtrVector;
public:
	Block(size_t count=0);
	~Block();
	void free_block_cholmod(cholmod_common *cm);
	void LU_decomposition();
	void CK_decomp(Matrix & A, cholmod_common *cm);
	void solve_CK(cholmod_common *cm); // solve with cholesky decomp
	void solve_LU();

	// allocate space for the matrix and vectors accoding to count
	void allocate_resource(cholmod_common*cm);
	void allocate_resource_LU();
	void allocate_mpi_resource(cholmod_common *cm);
	void get_col_compressed(Matrix & A);
	void update_x();

	void update_rhs(int &my_id);

	void update_block_geometry(MPI_CLASS &mpi_class);
	
	NetPtrVector boundary_netlist;

	// ****** variable for CK solver ****
	cholmod_factor * L;
	
	// vector b
	cholmod_dense * b_ck, *b_new_ck;
	// pointer to b_ck, b_new_ck, and x_ck;
	double *bp, *bnewp, *xp, *x_old;
	// solution
	cholmod_dense *x_ck;
	
	// ****** variable for CK solver *****

	// ****** variable for LU solver *****
	void *Numeric;
	UF_long * Ap;   // n+1
	UF_long * Ai;   // nz
	double * Ax;    // nz

	//vector b
	Vec b, bnew;
	// solution
	Vec x, xold;
	// ******** variable for LU solver ******

	// equal to matrix size and b size
	size_t count;

	Node ** nodes;

	float lx, ly, ux, uy;
};

#endif
