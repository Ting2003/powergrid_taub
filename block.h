// ----------------------------------------------------------------//
// Filename : block.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of block class
// ----------------------------------------------------------------//

#ifndef __BLOCK_H__
#define __BLOCK_H__
#include <fstream>
#include <map>
#include "triplet.h"
#include "global.h"
#include "vec.h"
#include "net.h"
#include "util.h"
#include "transient.h"
//#include "umfpack.h"
#include "cholmod.h"
#include "mpi_class.h"

using namespace std;

class Block{
	typedef vector<Net *> NetPtrVector;
	typedef vector<Node *> NodePtrVector;
	typedef NetPtrVector NetList;
public:
	Block(size_t count=0);
	~Block();
	void free_block_cholmod(cholmod_common *cm);
	void LU_decomposition();
	void CK_decomp(Matrix & A, cholmod_common *cm);
	void solve_CK_DC(int my_id); // solve with cholesky decomp

	void solve_CK_tr(); // solve with cholesky decomp
	// allocate space for the matrix and vectors accoding to count
	void allocate_resource();
	void allocate_mpi_resource(cholmod_common *cm);

	void update_x();

	bool node_in_block(Node *nd);
	int net_in_block(Net *net);

	void sort_nodes();
	void copy_node_voltages_block();
	double modify_voltage(int &my_id);
	void copy_array(double *x_old, double *xp);
	void build_nd_IdMap();
	void stamp_matrix(int &my_id, MPI_CLASS &mpi_class);
	void stamp_matrix_tr(int &my_id, MPI_CLASS &mpi_class, Tran &tran);
	void stamp_resistor(int &my_id, Net * net);
	void stamp_resistor_tr(int &my_id, Net * net);

	void make_A_symmetric(double *bp, int &my_id);
	void stamp_VDD(int &my_id, Net * net);
	void stamp_VDD_tr(int &my_id, Net * net);
	void stamp_current(int &my_id, Net * net, MPI_CLASS &mpi_class);
	void stamp_inductance_dc(Net * net, int &my_id);
	void stamp_capacitance_tr(Net *net, Tran &tran, int &my_id);
	void stamp_inductance_tr(Net * net, Tran &tran, int &my_id);

	void make_A_symmetric_tr(int &my_id, Tran &tran);
	void stamp_current_tr(int &my_id, double &time);
	void stamp_current_tr_net(Net * net, double &time, int &my_id);

	void stamp_bd_net(int my_id, Net *net);
	void clear_A();
	void CK_decomp();
	void copy_vec(double *bnewp, double *bp);
	void reset_array(double *bp);
	void modify_rhs_tr_0(double * b, double *x, int &my_id);
	void modify_rhs_c_tr_0(Net *net, double * rhs, double *x, int &my_id);

	void modify_rhs_l_tr_0(Net *net, double *rhs, double *x, int &my_id);

	void modify_rhs_c_tr(Net *net, double * rhs, double *x);

	void modify_rhs_l_tr(Net *net, double *rhs, double *x);

	void modify_rhs_tr(double * b, double *x);

	void stamp_current_tr_1(double &time);
	void current_tr(Net *net, double &time);
	void stamp_current_tr_net_1(Net * net, double &time);
	void update_rhs(double *bnewp, double *bp, int &my_id);
	// NetPtrVector boundary_netlist;

	cholmod_common c, *cm;
	Matrix A;
	cholmod_factor * L;
	
	// vector b
	cholmod_dense * b_ck, *b_new_ck, *bnew_temp;
	// pointer to b_ck, b_new_ck, and x_ck;
	double *bp, *bnewp, *xp, *x_old, *bnewp_temp;
	// solution
	cholmod_dense *x_ck;

	double *Lx;
	int *Li, *Lp, *Lnz;

	// number of *representative* nodes in this block
	// equal to matrix size and b size
	size_t count;

	// Node ** nodes;
	NodePtrVector replist;
	map<Node *, size_t> nd_IdMap;
	Node *nd_GND;
	NetList net_set[NUM_NET_TYPE];// should be the same as size of NET_TYPE
	NetPtrVector bd_netlist;

	double lx, ly, ux, uy;
	friend class Circuit;
};

bool compare_node_ptr(const Node *a, const Node *b);	
#endif
