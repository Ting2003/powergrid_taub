// ----------------------------------------------------------------//
// Filename : block.cpp
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// implementation of block class
// ----------------------------------------------------------------//
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added opeartor !=
//   * added this log

#include <cassert>
#include "cholmod.h"
#include "block.h"
#include "node.h"
#include "util.h"
#include "umfpack.h"

Block::Block(size_t _count):
        b_ck(NULL),
	b_new_ck(NULL),
	x_old(NULL),
	x_new(NULL),
	b_init(NULL),
	b_new(NULL),
	//bp(NULL),
	//bnewp(NULL),
	//xp(NULL),
	x_ck(NULL),	
	count(_count),
	nodes(NULL),
	lx(-1.0), ly(-1.0),
	ux(-1.0), uy(-1.0){}

Block::~Block(){
    delete [] nodes;
    delete [] x_old;
    delete [] b_init;
    delete [] b_new;
    delete [] x_new;
}

void Block::free_block_cholmod(cholmod_common *cm){
    cholmod_free_factor(&L, cm);
    cholmod_free_dense(&b_ck, cm);
    cholmod_free_dense(&b_new_ck, cm);
    cholmod_free_dense(&x_ck, cm);
}

void Block::CK_decomp(Matrix & A, cholmod_common *cm){
	Algebra::CK_decomp(A, L, cm);
}

void Block::solve_CK(cholmod_common *cm){
	x_ck = cholmod_solve(CHOLMOD_A, L, b_new_ck, cm);
}

void Block::allocate_resource(cholmod_common *cm){
	if( count == 0 ) return;
	nodes = new Node *[count];

	b_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	x_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	bp = static_cast<double*>(b_ck->x);
	xp = static_cast<double*>(x_ck->x);
	b_new_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	bnewp = static_cast<double*>(b_new_ck->x);
	x_old = new double [count];
	x_new = new float [count];
}

// update rhs of each block with its boundary netlist
void Block::update_rhs(int &my_id){
	size_t size = boundary_netlist.size();
	size_t k=0, l=0;
	//b_new = b;
	for(size_t i=0;i<count;i++){
		bnewp[i] = bp[i];
	}

	// for each net in this block
	for(size_t i=0;i<size;i++){
		Net * net = boundary_netlist[i];
		double G = 1.0/net->value;

		Node * a = net->ab[0]->rep;
		Node * b = net->ab[1]->rep;

		// if a is inside block
		if(a->flag_bd == 0){
			k = a->rid;
			if(!a->isX())
				bnewp[k] += G * b->value;
		}
		else if(b->flag_bd ==0){
			l = b->rid;
			if(!b->isX()) //b_new[l] += G *a->value;
				bnewp[l] += G * a->value;
		}
	} // end of for i
}

/////////////////////////////////////////////////////////////////
// methods for BlockInfo

// update block 4 corners
void Block::update_block_geometry(MPI_CLASS &mpi_class){
	// compute the geometrical information for the blocks
	lx = mpi_class.block_geo[0];
	ly = mpi_class.block_geo[1];
	ux = mpi_class.block_geo[2];
	uy = mpi_class.block_geo[3];
}
