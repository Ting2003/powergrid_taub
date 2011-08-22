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
	Numeric(NULL),
	Ap(NULL),
	Ai(NULL), 
	Ax(NULL),

	/*b_ck(NULL),
	b_new_ck(NULL),
	x_old(NULL),
	x_ck(NULL),
	*/
	count(_count),
	nodes(NULL),
	lx(-1.0), ly(-1.0),
	ux(-1.0), uy(-1.0){}

Block::~Block(){
	umfpack_dl_free_numeric (&Numeric) ;                             delete [] Ap;
	delete [] Ai;
	delete [] Ax; 
	delete [] nodes;
    	//delete [] x_old;
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

void Block::LU_decomposition(){
	Algebra::LU_decomposition(count, Ap, Ai, Ax, &Numeric);
}

void Block::solve_CK(cholmod_common *cm){
	x_ck = cholmod_solve(CHOLMOD_A, L, b_new_ck, cm);
	//cholmod_solve_new(CHOLMOD_A, L, b_new_ck, x_ck, cm);
}

void Block::solve_LU(){ 
	//the matrix needs to be initialized 
	assert( Ap != NULL && Ai != NULL && Ax != NULL);             
	int status; 
	double Control [UMFPACK_CONTROL];                            
	umfpack_dl_defaults (Control) ; 
	double *null = (double *) NULL;                              
	status = umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax, 
		x.get_val(), bnew.get_val(), Numeric, Control, null) ;               
	if( status < 0 ){ 
		umfpack_dl_report_status (Control, status) ;         
		report_exit("umfpack_dl_solve failed\n") ;           
	}
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
}

void Block::allocate_resource_LU(){
	if( count == 0 ) return;
	nodes = new Node *[count];
	
	b.resize (count);
	b.fill_zeros();
	bnew.resize(count);
	bnew.fill_zeros();

	x.resize(count);
	x.fill_zeros();
	xold.resize(count);
	xold.fill_zeros();
}

void Block::get_col_compressed(Matrix & A){
	size_t nz = A.size();
	size_t n_row = count;
	size_t n_col = count;

	// NOTE: DO NOT MODIFY. size must be n_col+1, see UMFPACK manual
	UF_long * Ti = new UF_long[nz];
	UF_long * Tj = new UF_long[nz];
	double * Tx = new double[nz];
	A.to_arrays((size_t*)Ti,(size_t*)Tj,Tx);

	Ap = new UF_long[n_col+1]; 
	Ai = new UF_long[nz];
	Ax = new double [nz];

	int status;
	double Control [UMFPACK_CONTROL];
	umfpack_dl_defaults (Control) ;
	status = umfpack_dl_triplet_to_col(n_row, n_col, nz, Ti, Tj, Tx, 
			Ap, Ai, Ax, (UF_long *) NULL);
	if( status < 0 ){
		umfpack_dl_report_status (Control, status) ;
		report_exit("umfpack_dl_triplet_to_col failed\n") ;
	}

	delete [] Ti;
	delete [] Tj;
	delete [] Tx;
}

// update rhs of each block with its boundary netlist
void Block::update_rhs(int &my_id){
	size_t size = boundary_netlist.size();
	size_t k=0, l=0;

	bnew = b;
	int temp = 1;
		
	// for each net in this block
	for(size_t i=0;i<size;i++){
		Net * net = boundary_netlist[i];
		double G = 1.0/net->value;

		Node * a = net->ab[0]->rep;
		Node * c = net->ab[1]->rep;
	
		// if a is inside block
		if(a->flag_bd == 0){
			k = a->rid;
			if(!a->isX()){
				bnew[k] += G * c->value;
			}
		}
		else if(c->flag_bd ==0){
			l = c->rid;
			if(!c->isX()){
				bnew[l] += G * a->value;
			}
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
