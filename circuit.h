// ----------------------------------------------------------------//
// Filename : circuit.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of Circuit class
// used to construct the circuit network
// ----------------------------------------------------------------//
// - Ting Yu - Tue Feb 8 5:45 pm 2011
//   * added the ostream<< func in .h file
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added solve() and related function
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#ifndef __CIRCUIT_H__
#define __CIRCUIT_H__

#include <string>
#include <vector>
#include <tr1/unordered_map>
#include <map>
#include <list>
#include <cmath>
#include "cholmod.h"
#include "global.h"
#include "node.h"
#include "net.h"
#include "vec.h"
#include "triplet.h"
#include "block.h"
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
using namespace std;
using namespace std::tr1;

typedef vector<double> DoubleVector;
typedef vector<Net *> NetPtrVector;
typedef list<Net *> NetPtrList;
typedef vector<Node *> NodePtrVector;
typedef NetPtrVector NetList;

// functor of translating Node * to void *
namespace std{ 
	namespace tr1{
		template<> struct hash< Node * >{
			size_t operator()( const Node * x ) const {
				return hash< const char* >()( (char*) x );
			}
		};
	}
}

class Circuit{
public:
	Circuit(string name="");
	~Circuit();
	void check_sys() const;
	friend class Block;
	// can be written as inline to speed up
	Node * get_node(string name);
	Net * get_net(string name);
	string get_name() const;

	static size_t get_total_num_layer();

	// add a node into nodelist
	bool add_node(Node * nd);

	// add a net into netset
	bool add_net(Net * net);

	bool has_node(string name) const;
	bool has_net(string name) const;

	// sort nodes according to predefined order
	void sort_nodes();

	// solve for node voltage
	void solve(int &my_id, int &num_procs, MPI_Comm &comm);
	
	void set_blocklist(Node * nd);

	static void set_parameters(double, double, double, size_t, int);
	static void get_parameters(double&, double&, double&, size_t&, int&);

	friend ostream & operator << (ostream & os, const Circuit & ckt);
	friend class Parser;

	// C style output
	void print(int &color);
	cholmod_common c, *cm;

	// mpi related functions
	void mpi_setup(int &my_id, int &num_procs, MPI_Comm &comm);
	void MPI_Assign_Task(int &num_tasks, int &num_procs);
	void block_mpi_setup(int &num_procs);
	void mpi_create_matrix(int &my_id, int &num_procs, MPI_Comm &comm); 

	// mpi related variables
	// member
	// ******* processor 0  variable ********
	Matrix *A_g; // local is Matrix *A
	vector<long> Ti_g; // receiv by Ti
	vector<long> Tj_g; // receiv by Tj
	vector<double> Tx_g; // receiv by Tx

	// start_task and end_task
	int *start_task;
	int *end_task;
	// tasks_n stores # of blocks for each processor
	int *tasks_n;
	// b_new_info is the global solution array
	float *b_new_info;

	// # of nz for all the blocks, receive by 
	int *nz_A;
	// # of n for all the blocks, receive by L_n
	int *L_n_d;

	// base for matrix nz for each processor
	int *base_nz_A;
	// base for rhs n for each processor
	int *base_n_d;

	// # for matrix nz send for each processor
	int *send_nz_A;
	// # for rhs n send for each processor
	int *send_n;
	
	// ****** other processor *******
	Matrix *A;
	vector<long> Ti;
	vector<long> Tj;
	vector<double> Tx;

	BlockInfo block_info_mpi;

	// solution array for each processor
	float *b_x_d;

	// # of matrix nz in each block within a processor
	int *nz_dd_A;
	// # of rhs n in each block within a processor
	int *L_n_dd;
	
	// nz and L_n are # of nz in L and # of n in rhs
	// nz here is the non-zero of triplet including coord
	int nz; 
	int L_n;
	// block_size is the number of tasks for each cpu
	int block_size;

	int num_blocks;

	size_t total_n, total_nz_A;

private:
	// member functions
	void solve_LU();
	void solve_LU_core();

	bool solve_IT(int &my_id, int&num_procs, MPI_Comm &comm);
	void solve_block_LU();
	void decomp_matrix(int &my_id, Matrix *A);
	bool solve_pcg();
	//bool solve_block_pcg();

	
	// initialize things before solve_iteration
	void solve_init();

	// updates nodes value in each iteration
	double solve_iteration(int &my_id, int&num_procs, MPI_Comm &comm);

	void block_init(Matrix *A);
	void update_block_geometry();

	// methods of stamping the matrix
	void stamp_by_set(Matrix & A, double* b);
	void stamp_resistor(Matrix & A, Net * net);
	void stamp_current(double* b, Net * net);
	void stamp_VDD(Matrix & A, double* b, Net * net);
	
	void make_A_symmetric(Matrix &A, double *bp);
	void make_A_symmetric_block();

	void stamp_block_matrix(Matrix *A);
	void stamp_boundary_matrix();
	void stamp_boundary_net(Net * net);
	void stamp_block_resistor(Net *net, Matrix * A);
	void stamp_block_current(Net * net, Matrix * A);
	void stamp_block_VDD(Net * net, Matrix * A);

	void update_block_rhs(Block & block, int dir);

	//  ******* method for PCG method  ********
	// solve circuit with preconditioned pcg method
	void copy_node_voltages_block(bool from=true);

	// after solving, copy node voltage from replist to nodes
	void get_voltages_from_LU_sol(double* x);
	void get_voltages_from_block_LU_sol();
	void get_vol_mergelist();

	Vec compute_precondition(const Matrix & ML, const Matrix & D, 
			const Matrix & MU, Vec &r);
	void init_precondition(const Matrix &A, Matrix &ML, Matrix &D,
			Matrix &MU);
	// searching all the blocks for global residue
	Vec get_block_res(const Vec& b, const Vec& xk1);
	// searching all the blocks for global zk1
	Vec compute_block_precondition( Vec &r);

	void set_len_per_block();
	void find_block_size ();
	void block_boundary_insert_net(Net * net);
	void find_block_base();

	void partition_circuit();
	double modify_voltage(int &my_id, Block &block, size_t &base, double* x_old);

	void node_voltage_init();
	void solve_one_block(size_t block_id);

	void select_omega();

	void set_type(CIRCUIT_TYPE type){circuit_type = type;};

	void get_samples();

	bool check_diverge() const;

	void merge_along_dir(Node *, DIRECTION dir);
	Node * merge_along_dir_one_pass(Node *, DIRECTION dir, bool remove);
	void merge_node(Node * node);

	// ************** member variables *******************
	NodePtrVector nodelist;		// a set of nodes
	NodePtrVector replist;		// a set of representative nodes
	NodePtrVector mergelist;	// nodes for merging
	NetList net_set[NUM_NET_TYPE];// should be the same as size of NET_TYPE
	// defines the net direction in layers
	static vector<LAYER_DIR> layer_dir;
	vector<int> layers;
	
	// mapping from name to Node object pointer
	unordered_map<string, Node*> map_node;

	// mapping from Net pointer to their index in netlist
	unordered_map<Net*, size_t> net_id;

	// circuit name
	string name;

	// blocks
	BlockInfo block_info;
	size_t x_min, y_min, x_max, y_max;

	// control variables
	static double EPSILON;
	static double OMEGA;
	static double OVERLAP_RATIO;
	static size_t MAX_BLOCK_NODES;
	static int MODE; // 0 = IT, 1 = LU

	CIRCUIT_TYPE circuit_type;

	NodePtrVector sample;

	double VDD;
};

inline size_t Circuit::get_total_num_layer(){return layer_dir.size();}

// adds a node into nodelist
inline bool Circuit::add_node(Node * node){
	nodelist.push_back(node);
	map_node[node->name] = node;
	return true;
}

// adds a net into netset
inline bool Circuit::add_net(Net * net){
	if( net->type == RESISTOR )
		net_id[net] = net_set[net->type].size();
	net_set[net->type].push_back(net);
	return true;
}

// fina a node by name
inline bool Circuit::has_node(string name) const{
	if( map_node.find(name) != map_node.end() ) return true;
	return false;
}

// get a node by name
inline Node * Circuit::get_node(string name){
	unordered_map<string, Node*>::const_iterator it = map_node.find(name);
	if( it != map_node.end() ) return it->second;
	else return NULL;
}

inline void Circuit::merge_node(Node * node){
	for(DIRECTION dir = WEST; dir <= NORTH; dir=DIRECTION(dir+1)){
		// test whether this line has been processed
		if( node->end[dir] != node ) continue;

		// probe for one step, if the line is only one step, don't do it.
		Node * next = node->get_nbr_node(dir);
		if( next == NULL || !next->is_mergeable() ) continue;
		merge_along_dir(node, dir);
	}
}

/*
// find a net by name
inline bool Circuit::has_net(string name) const{
	if( map_net.find(name) != map_net.end() ) return true;
	return false;
}


// get a net by name
inline Net * Circuit::get_net(string name){return map_net[name];}
*/

bool compare_node_ptr(const Node *a, const Node *b);
ostream & operator << (ostream & os, const NodePtrVector & nodelist);
ostream & operator << (ostream & os, const NetList & nets);
//ostream & operator << (ostream & os, const vector<Block > & block_info);

#endif
