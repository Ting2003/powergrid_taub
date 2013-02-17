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
//#include "umfpack.h"
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

	void solve_CK_tr(cholmod_common *cm); // solve with cholesky decomp
	// allocate space for the matrix and vectors accoding to count
	void allocate_resource(cholmod_common*cm);
	void allocate_mpi_resource(cholmod_common *cm);

	void update_x();

	void update_rhs(double *bnewp, double *bp, int &my_id);

	void update_block_geometry(MPI_CLASS &mpi_class);
	
	NetPtrVector boundary_netlist;

	cholmod_factor * L;
	
	// vector b
	cholmod_dense * b_ck, *b_new_ck, *bnew_temp;
	// pointer to b_ck, b_new_ck, and x_ck;
	double *bp, *bnewp, *xp, *x_old, *bnewp_temp;
	// solution
	cholmod_dense *x_ck;

	// number of *representative* nodes in this block
	// equal to matrix size and b size
	size_t count;

	Node ** nodes;

	float lx, ly, ux, uy;
};

#endif
