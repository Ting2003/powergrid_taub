// ----------------------------------------------------------------//
// Filename : circuit.cpp
// Author : Zigang Xiao <zxiao2@illinois.edu>
//          Ting Yu <tingyu1@illinois.edu>
//
// implementation file of circuit.h
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 30 18:35:56 CST 2011
//   * Add UMFPACK support
// - Zigang Xiao - Tue Jan 25 17:19:21 CST 2011
//   * added framework of PCG
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added solve() and related function
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <utility>
#include <cassert>
#include <vector>
#include "cholmod.h"
#include "umfpack.h"
#include "circuit.h"
#include "util.h"
#include "algebra.h"
#include "node.h"
#include "mpi.h"

using namespace std;

double Circuit::EPSILON = 1e-5;
size_t Circuit::MAX_BLOCK_NODES =5500;
double Circuit::OMEGA = 1.2;
double Circuit::OVERLAP_RATIO = 0.2;
int    Circuit::MODE = 0;
//const int MAX_ITERATION = 1000;
const int SAMPLE_INTERVAL = 5;
const size_t SAMPLE_NUM_NODE = 10;
const double MERGE_RATIO = 0.3;

//////////////////////////////////////////////////////////////////////////
// Constructor and utility functions goes here

vector<LAYER_DIR> Circuit::layer_dir(MAX_LAYER);

// constructor of Circuit class, name is optional
Circuit::Circuit(string _name):name(_name),
	x_min(INFTY),y_min(INFTY),x_max(0),y_max(0),
	circuit_type(UNKNOWN), VDD(0.0){
	// add ground node
	Node * gnd = new Node(string("0"), Point(-1,-1,-1));
	gnd->rep = gnd;
	this->add_node(gnd);

	for(int i=0;i<MAX_LAYER;i++)
		layer_dir[i]=NA;
	peak_mem = 0;
	CK_mem = 0;
}

// Trick: do not release memory to increase runtime
Circuit::~Circuit(){
	for(size_t i=0;i<nodelist.size();i++) delete nodelist[i];
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		for(it=ns.begin();it!=ns.end();++it)
			delete *it;
	}
}

void Circuit::check_sys() const{
	clog<<"**** CHECKING SYSTEM ENVIRONMENT ****"<<endl;
	clog<<"* int size     = "<< sizeof(int)<<endl;
	clog<<"* long size    = "<< sizeof(long)<<endl;
	clog<<"* size_t size  = "<< sizeof(size_t)<<endl;
	clog<<"* UF_long size = "<< sizeof(UF_long)<<endl;
	clog<<"* Max nodelist = "<<(size_t)nodelist.max_size()<<endl;
	clog<<"****            END              ****"<<endl<<endl;
}

// functor to be used in STL sort function
// order: y > x > z > flag 
// input: two node a, b
// return true if a < b, false o/w
// note that ground note are put to last
bool compare_node_ptr(const Node * a, const Node * b){
	if( a->is_ground() ) return false;
	if (b->is_ground() ) return true;

	if( a->pt.y == b->pt.y ){
		if( a->pt.x == b->pt.x ){
			if( a->pt.z == b->pt.z ){
				return a->isX();
			}
			else{
				return (a->pt.z > b->pt.z);// top down
			}
		}
		else
			return ( a->pt.x < b->pt.x );
	}
	else
		return (a->pt.y < b->pt.y);
}

// sort the nodes according to their coordinate 
void Circuit::sort_nodes(){
	//clog<<" before sorting. "<<endl;
	
	sort(nodelist.begin(), nodelist.end(), compare_node_ptr);
	// update node id mapping, 
	// NOTE: ground node will be the last
}

string Circuit::get_name() const{return this->name;}

ostream & operator << (ostream & os, const NodePtrVector & nodelist){
	for(size_t i=0;i<nodelist.size();i++)
		os<<*nodelist[i]<<endl;
	return os;
}

ostream & operator << (ostream & os, const NetList & nets){
	NetList::const_iterator it;
	for(it=nets.begin();it!=nets.end();++it)
		if( (*it) != NULL ) os<<**it<<endl;
	return os;
}

ostream & operator << (ostream & os, const Circuit & ckt){
	os<<"Circuit ["<<ckt.name<<"] info:"<<endl;

	os<<"==== Nodes ===="<<endl;
	os<<ckt.nodelist;

	os<<"==== Reps  ===="<<endl;
	os<<ckt.replist;

	os<<"==== Nets  ===="<<endl;
	os<<ckt.net_set[RESISTOR];

	return os;
}

void Circuit::print(){
	// uncomment this if want to output to a file
	//freopen("output.txt","w",stdout);

	// don't output ground node
	for(size_t i=0;i<nodelist.size()-1;i++){
		printf("%s  %.5e\n", nodelist[i]->name.c_str(), 
				nodelist[i]->value);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Computation Functions

// initialization before solving the circuit
// 1. sort the nodes
// 2. set node representatives
// 3. find node in which block, update count
// 4. get representative lists
void Circuit::solve_init(){
	sort_nodes();
	size_t size = nodelist.size() - 1;
	Node * p = NULL;
	for(size_t i=0, nr=0;i<size;i++){
		p=nodelist[i];
		//clog<<i<<" "<<*p<<endl;

		// test if it can be merged
		if( p->is_mergeable() ){
			mergelist.push_back(p);
			continue;
		}
		Net * net = p->nbr[TOP];
		merge_node(p);

		// find the VDD value
		if( p->isX() ) VDD = p->get_value();
		// test short circuit
		if( !p->isX() && // X must be representative 
		    net != NULL &&
		    fzero(net->value) ){
			// TODO: ensure ab[1] is not p itself
			assert( net->ab[1] != p );
			p->rep = net->ab[1]->rep;
		} // else the representative is itself
		// push the representatives into list
		if( p->rep == p ) {
			replist.push_back(p);
			p->rid = nr++;
		}
	}// end of for i
	size_t n_merge = mergelist.size();
	size_t n_nodes = nodelist.size();
	size_t n_reps  = replist.size();
	double ratio = n_merge / (double) (n_merge + n_reps);
	/*clog<<"mergeable  "<<n_merge<<endl;
	clog<<"replist    "<<n_reps <<endl;
	clog<<"nodelist   "<<n_nodes<<endl;
	clog<<"ratio =    "<<ratio  <<endl;*/
	net_id.clear();
}

// partition the circuit to X_BLOCKS * Y_BLOCKS blocks
// according to the node size. Assuming they are distributed
// uniformly at random
void Circuit::partition_circuit(){
	size_t num_nodes = replist.size();
	size_t num_blocks =  num_nodes / MAX_BLOCK_NODES;
	if( num_nodes % MAX_BLOCK_NODES > 0 ) ++num_blocks;
	size_t len_x = x_max-x_min;
	size_t len_y = y_max-y_min;
	size_t X_BLOCKS, Y_BLOCKS;

	// Extreme case: only one node in x/y axis
	if( num_blocks == 1 ){
		X_BLOCKS = Y_BLOCKS = 1;
	}
	else if( len_x == 0 ){
		X_BLOCKS = 1;
		Y_BLOCKS = num_blocks;
	}
	else if (len_y == 0 ){
		Y_BLOCKS = 1;
		X_BLOCKS = num_blocks;
	}
	else{// compute length of x and y according to their ratio
		double ratio = double(len_x) / double(len_y);
		X_BLOCKS = sqrt(num_blocks/ratio);
		Y_BLOCKS = X_BLOCKS*ratio;
		if( X_BLOCKS * Y_BLOCKS < num_blocks ) Y_BLOCKS+=1; 
		num_blocks = X_BLOCKS * Y_BLOCKS;
	}
	//clog<<"num_blocks: "<<X_BLOCKS<<" / "<<Y_BLOCKS <<endl;
	block_info.X_BLOCKS = X_BLOCKS;
	block_info.Y_BLOCKS = Y_BLOCKS;
	block_info.resize(X_BLOCKS * Y_BLOCKS);
}

// build up block info
// 1. Find block divide point
// 2. Set each node in replist into block
// 3. Compute block size
// 4. Insert boundary netlist into map
void Circuit::block_init(){
	block_info.set_len_per_block(x_min, x_max, y_min, y_max, OVERLAP_RATIO);
	block_info.update_block_geometry();
	find_block_size();
	copy_node_voltages_block();

	//stamp_boundary_matrix();
	stamp_block_matrix();
}

void Circuit::block_boundary_insert_net(Net * net){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
	if(nd[0]->is_ground() || nd[1]->is_ground()) return;

	vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};

	// idea: for each block where node_a is in
	// find whether node_b is also in the block
	// if yes, then this net is not a boundary net
	// for this block
	vector<size_t>::const_iterator it;
	for(size_t j=0;j<2;j++){
		vector<size_t> *p = ls[j], *q = ls[1-j];
		// for each block k in ls[j]
		for(size_t k=0;k<p->size();k++){
			// find whether block_id in another list
			size_t block_id = (*p)[k];
			it = find( (*q).begin(), (*q).end(), block_id);

			// block_id not in another list, net is boundary
			if( it == (*q).end() ){
				Block & blk = block_info[block_id];
				blk.boundary_netlist.push_back(net);
			}
		}// end of for k
	}// end of for j
}

// stamp the nets by sets, block version
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_block_matrix(){
	size_t num_blocks = block_info.size();
	Matrix A[num_blocks];
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			for(it=ns.begin();it!=ns.end();++it){
				Net * net = *it;
				if( net == NULL ) continue;
				if(net->ab[0]->is_ground() || 
				   net->ab[1]->is_ground()) 
					continue;
				assert( fzero(net->value) == false );
				stamp_block_resistor(*it, A);
			}
			break;
		case CURRENT:
			for(it=ns.begin();it!=ns.end();++it)
				stamp_block_current((*it), A);
			break;
		case VOLTAGE:
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_block_VDD((*it), A);
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	make_A_symmetric_block();
	//size_t count = 0;
	// after stamping, convert A to column compressed form
	for(size_t i=0;i<num_blocks;i++){
		if(block_info[i].count>0){
			A[i].set_row(block_info[i].count);
			block_info[i].CK_decomp(A[i], cm, peak_mem, CK_mem);
			block_info[i].solve_CK_setup();
			//count += block_info[i].L_h_nz;
		}
	}
	//clog<<"peak memory for cholmod: "<<peak_mem / 1e9<<" e+06"<<endl;
	//clog<<"CK_mem is: "<<CK_mem / 1e6 <<" G"<<endl;
}

// 1. mark rep nodes into corresponding blocks
// 2. find block size of replist nodes
// 3. allocate rhs size of each block
// 4. find local block index for each node
void Circuit::find_block_size(){
	// start assign replist node into blocks
	size_t n = replist.size();
	Node * p = NULL;
	
	// find which blocks the node belongs to,
	// and its index in the block
	for(size_t i=0;i<n;i++){
		p=replist[i];
		set_blocklist(p);
	}
	// for each block, allocate resource
	for(size_t i=0;i<block_info.size();i++){
		Block & block = block_info[i];
		block.allocate_resource(cm);
		//for(int k=0;k<4;k++) block.A_nbr[k].set_row(block.count);
	}
}

bool Circuit::solve(int &my_id, int&num_procs){
	if( MODE == 0 ){
		int flag = 0;
		flag = solve_IT_setup(my_id, num_procs);
		// if solve with block-iterative method
		// else do not continue
		if(flag !=0)	return true;
		else return false;
	}
	else{
		solve_LU();
		return false;
	}
}

// solve Circuit
int Circuit::solve_IT_setup(int &my_id, int&num_procs){
	// did not find any `X' node
	if( circuit_type == UNKNOWN )
		circuit_type = C4;
	solve_init();
	/*if( replist.size() <= 2*MAX_BLOCK_NODES ){
	  clog<<"Replist is small, use direct LU instead."<<endl;
	  solve_LU_core();
	  return 0;
	  }
	  */
	select_omega();
	partition_circuit();
	// Only one block, use direct LU instead
	if( block_info.size() == 1 ){
		//clog<<"Block size = 1, use direct LU instead."<<endl;
		solve_LU_core();
		return 0;
	}
	
	// cm declared in circuit class
	// cholmod_common c, *cm;
	cm = &c;
	cholmod_start(cm);
	cm->print = 5;
	cm->final_ll = true;
	block_init();

	if(my_id==0) {
		//for(size_t i=0;i<block_info.size();i++)
			//clog<<i<<" "<<block_info[i].count<<endl;
		clog<<"num_blocks: "<<block_info.X_BLOCKS<<" / "<<block_info.Y_BLOCKS <<endl;
	}
	return 1;
	/*clog<<"e="<<EPSILON
	  <<"\to="<<OMEGA
	  <<"\tr="<<OVERLAP_RATIO
	  <<"\tb="<<MAX_BLOCK_NODES
	  <<"\tmode="<<MODE<<endl;
	  */
	
	/*int iter = 0;	
	double diff=0;
	bool successful = false;

	//clog<<"start and end_task for: "<<
	//my_id<<" "<<start_task<<" "<<end_task<<endl;
	float time=0;
	double t1, t2;
	t1= MPI_Wtime();
	while( iter < MAX_ITERATION ){
		diff = solve_iteration(my_id, num_procs, 
			mpi_class, L_n, num_blocks, 
			block_size, iter);
		iter++;
		//if(my_id ==0)
			//clog<<"iter, diff: "<<iter<<" "<<diff<<endl;
		if( diff < EPSILON ){
			successful = true;
			break;
		}
	}
	
	t2 = MPI_Wtime();
	time = t2-t1;
	if(my_id==0){
		clog<<"solve iteration time for "<<my_id<<" is "<<time<<endl;
		clog<<"# iter: "<<iter<<endl;
		get_voltages_from_block_LU_sol();
		get_vol_mergelist();
	}
	for(size_t i=0;i<block_info.size();i++){
		if(block_info[i].count > 0)
			block_info[i].free_block_cholmod(cm);
	}
	
	cholmod_finish(cm);
	return successful;*/
}

// TODO: add comment
void Circuit::node_voltage_init(){
	for(size_t i=0;i<block_info.size();i++){
		Block & block = block_info[i];
		for(size_t j=0;j<block.count;j++){
			block.xp[i] = VDD;
			//block.x[i] = VDD;
			block.nodes[i]->value = VDD;
		}
	}
}

double Circuit::modify_voltage(Block & block, double * x_old){
	double max_diff = 0.0;
	for(size_t i=0;i<block.count;i++){
		block.x_new[i] = (1-OMEGA) * x_old[i] + OMEGA * block.x_new[i];
		block.nodes[i]->value = block.x_new[i];
		double diff = fabs(x_old[i] - block.x_new[i]);
		if( diff > max_diff ) max_diff = diff;
	}
	return max_diff;
}

// stamp the matrix and solve
void Circuit::solve_LU_core(){
	size_t n = replist.size();	// replist dosn't contain ground node
	if( n == 0 ) return;		// No node
	//Vec b(n), x(n);
	cholmod_common c, *cm;
	cholmod_dense *b, *x;
	cm = &c;
	cholmod_start(cm);
	cm->print = 5;
	b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
	x = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
	double *bp, *xp;
	bp = static_cast<double *> (b->x);
	Matrix A;
	stamp_by_set(A, bp);
	make_A_symmetric(A, bp);
	//A.merge();

	A.set_row(replist.size());
	Algebra::solve_CK(A, x, b, cm, peak_mem, CK_mem);
	
	xp = static_cast<double *> (x->x);
	// Vec b contains result, copy it back to nodelist
	get_voltages_from_LU_sol(xp);
	get_vol_mergelist();
	cholmod_free_dense(&x, cm);
	cholmod_free_dense(&b, cm);
	cholmod_finish(&c);
}

// solve the node voltages using direct LU
void Circuit::solve_LU(){
	solve_init();
	solve_LU_core();
}

// given vector x that obtained from LU, set the value to the corresponding
// node in nodelist
void Circuit::get_voltages_from_LU_sol(double * x){
	for(size_t i=0;i<nodelist.size()-1;i++){
		Node * node = nodelist[i];
		size_t id = node->rep->rid;	// get rep's id in Vec
		double v = x[id];		// get its rep's value
		node->value = v;
	}
}

// compute value of mergelist nodes
void Circuit::get_vol_mergelist(){
	DIRECTION p, q;
	for(size_t i=0;i<mergelist.size();i++){
		Node * node = mergelist[i];
		// check direction
		if( node->nbr[WEST] != NULL ){
			p = WEST;
			q = EAST;
		}
		else{
			p = SOUTH;
			q = NORTH;
		}
		// assign vol value to node
		// left end node and its value
		double r1 = node->eqvr[p];
		double v1 = node->end[p]->value;
		//clog<<" left node: "<<r1<<" / "<<v1<<endl;
		//clog<<" left end "<<node->end[p]->name<<endl;
		// right end node and its value
		double r2 = node->eqvr[q];
		double v2 = node->end[q]->value;
		//clog<<"right node: "<<r2<<" / "<<v2<<endl;
		//clog<<"right end "<<node->end[q]->name<<endl;
		// value for node
		if(v1 > v2){
			node->value = v2 + (v1 - v2) * r2 / (r1 + r2);
		}
		else{
			node->value = v1 + (v2 - v1)  * r1 / (r1 + r2);
		}
		//clog<<" node "<<*node<<endl;
	}
}

// copy solution of block into circuit
void Circuit::get_voltages_from_block_LU_sol(){

	size_t block_id;
	for(size_t i=0;i<nodelist.size()-1;i++){
		Node * node = nodelist[i];
		if( node->is_mergeable() ) continue;
		block_id = node->rep->blocklist[0];
		Block &block = block_info[block_id];
		size_t id = node->rep->id_in_block[0];
		//if(node->name == "n7_9690300_60300") clog<<node->rep->name<<" "<<block.x_new[id]<<endl;
		//Vec &p = block.x;
		//double v = p[id];		// get its rep's value
		double v = block.x_new[id];
		node->value = v;
		//if(i<10) clog<<node->name<<" "<<node->rep->name<<" "<<node->value<<endl;
	}
}

// solve the circuit using PCG method
/*bool Circuit::solve_pcg(){
	solve_init();
	// remember to exclude ground node
	size_t n = nodelist.size()-1;

	// xk is for k-th iteration, and xk1 is for (k+1)-th iteration
	Vec b(n);
	Vec pk(n), pk1(n);
	Vec rk(n), rk1(n);
	Vec zk(n), zk1(n);
	Vec xk(n), xk1(n);
	double alpha, beta;
	bool successful = false;
	Matrix A, ML, MU, D;

	this->stamp_by_set(A, b);
	
	// initialize first iteration	
	init_precondition(A, ML, D, MU);
	copy_node_voltages(xk1, true);
	rk1 = b - A * xk1; 
	zk1 = compute_precondition(ML, D, MU, rk1);// solve M*z0=r0
	pk1 = zk1;
	int k=0; 
	while(k++ < MAX_ITERATION){
		// new values becomes old values
		zk = zk1;
		pk = pk1;
		xk = xk1;
		double diff=1.0; 

		alpha = (rk * zk) / (pk * A * pk);
		xk1 = xk + alpha * pk;
		rk1 = rk - alpha * (A * pk);
		diff = find_diff(rk1);
		if( diff < EPSILON ) {
			successful = true;
			break;
		}
		zk1 = compute_precondition(ML, D, MU, rk1);
		beta = (rk1 * zk1) / (rk * zk);
		pk1 = zk1 + beta * pk; 
	}
	if( successful ) copy_node_voltages(xk1, false);
	return successful;
}

// M1= (D/OMEGA+L); M2= (D/OMEGA+L'); D:diagonal
void Circuit::init_precondition(const Matrix &A, 
		Matrix &ML, Matrix & D, Matrix &MU){
	A.diagonal_split(ML, D, MU);
	ML *= 1/OMEGA;
	ML += D;

	MU *= 1/OMEGA;
	MU += D;
}

// solve M*z=r; 
// compute a preconditioner for matrix A = ML + D + MU
Vec Circuit::compute_precondition(const Matrix & ML, const Matrix & D, 
			const Matrix & MU, Vec &r){
	size_t n = nodelist.size()-1;
	Vec z1(n), z2(n), z(n);
		
	// solve (D/OMEGA+L)*z=r
	// D is the diagonal of M, L is the lower triangular
	//z1 = Algebra::solve(ML, r);	
	//Algebra::solve(ML, r, z1);	

	// solve (OMEGA/(2-OMEGA)*D^-1)*z=r: D is the diagonal of M,
	z2 = (2-OMEGA)/OMEGA * (D * z1);
	
	// solve (D/OMEGA+L')*z=r: D is the diagonal of M, L' is the upper 
	//z  = Algebra::solve(MU, r);
	//Algebra::solve(MU, r, z);
	return z;
}
*/

// 1. copy node voltages from the circuit to a Vec
//    from = true then copy circuit to x
//    else copy from x to circuit
// 2. map block voltage into global_index
void Circuit::copy_node_voltages_block(bool from){
	size_t id;
	if( from == true ){
		for(size_t i=0;i<replist.size();i++){
			Node *node = replist[i];
			const vector<size_t> &block_id = node->get_block_id();
			for(size_t j=0;j<block_id.size();j++){
				Block &block = block_info[block_id[j]];
				id = node->id_in_block[j];
				block.x_new[id] = replist[i]->value;
				//block.xp[id] = replist[i]->value;
				//block.x[id] = replist[i]->value;
				block.nodes[id] = replist[i];
			}
		}
	}
	else{
		for(size_t i=0;i<nodelist.size()-1;i++){
			Node * node = nodelist[i];
			const vector<size_t> &block_id = 
				node->rep->get_block_id();
			for(size_t j=0;j<block_id.size();j++){
				Block &block = block_info[block_id[j]];
				id = node->rep->id_in_block[j];
				node->value = block.xp[id];
				//Vec &p = block.x;
				//node->value = p[id];
			}
		}
	}
}

// stamp the net in each set, 
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_by_set(Matrix & A, double* b){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			for(it=ns.begin();it!=ns.end();++it){
				if( (*it) == NULL ) continue;
				assert( fzero((*it)->value) == false );
				stamp_resistor(A, (*it));
				//block_boundary_insert_net(ns[i]);
			}
			break;
		case CURRENT:
			for(it=ns.begin();it!=ns.end();++it)
				stamp_current(b, (*it));
			break;
		case VOLTAGE:
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_VDD(A, b, (*it));
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
}

void Circuit::make_A_symmetric(Matrix &A, double *b){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p, *q;

	for(it=ns.begin();it!=ns.end();it++){
		if( (*it) == NULL ) continue;
			assert( fzero((*it)->value) == false );
		// node a points to X node
		if((*it)->ab[0]->isX()){
			p = (*it)->ab[0]; q = (*it)->ab[1];
		}
		else if((*it)->ab[1]->isX()){
			p = (*it)->ab[1]; q = (*it)->ab[0];
		}
		else continue;
		size_t id = q->rep->rid;
		double G = 1.0 / (*it)->value;
		b[id] += p->value * G;
	}
}

void Circuit::make_A_symmetric_block(){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	vector<size_t> *p, *q;
	Node *nk, *nl;

	for(it=ns.begin();it!=ns.end();it++){
		if( (*it) == NULL ) continue;
			assert( fzero((*it)->value) == false );
	
		Net *net = *it;
		Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
		if(!nd[0]->isX() && !nd[1]->isX()) continue;
		double G;
		G = 1./net->value;
		vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
		vector<size_t >::const_iterator it_1;
		for(size_t j=0;j<2;j++){
			p = ls[j];
			q = ls[1-j];
			nk = nd[j];
		       	nl = nd[1-j];
			for(size_t i=0;i<p->size();i++){
				// find whether block_id in another list
				size_t block_id = (*p)[i];
				it_1 = find( (*q).begin(), (*q).end(), block_id);
				// 2 nodes in the same block
				if(it_1!=(*q).end() && !nk->isX()){
					Block &block = block_info[block_id];
					size_t k1 = nk->id_in_block[i];
					block.bp[k1] += G *(nl->value);
				}
			}
		}
	}
}
void Circuit::stamp_resistor(Matrix & A, Net * net){
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	if( !nk->isX() ) {
		A.push_back(k,k, G);
		if(!nl->isX() && l < k) // store lower triangular
			A.push_back(k,l,-G);
	}

	if( !nl->isX() ) {
		A.push_back(l,l, G);
		if(!nk->isX() && k < l) // store ower triangular
			A.push_back(l,k,-G);
	}
}

// stamp a current source
void Circuit::stamp_current(double * b, Net * net){
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground() && !nk->isX() ) { 
		size_t k = nk->rid;
		b[k] += -net->value;
	}
	if( !nl->is_ground() && !nl->isX() ) {
		size_t l = nl->rid;
		b[l] +=  net->value;
	}
}

// stamp a voltage source
void Circuit::stamp_VDD(Matrix & A, double * b, Net * net){
	// find the non-ground node
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	size_t id = X->rep->rid;
	A.push_back(id, id, 1.0);
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
	    south->type == CURRENT ){
		// this node connects to a VDD and a current
		assert( feqn(1.0, b[id]) ); // the current should be stamped
		b[id] = net->value;	    // modify it
	}
	else
		b[id] += net->value;
}

// =========== stamp block version of matrix =======

void Circuit::stamp_block_resistor(Net * net, Matrix * A){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};

	double G;	
	G = 1./net->value;

	vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
	vector<size_t>::const_iterator it;
	for(size_t j=0;j<2;j++){
		vector<size_t> *p = ls[j], *q = ls[1-j];
		Node *nk = nd[j], *nl = nd[1-j];
		for(size_t i=0;i<p->size();i++){
			// find whether block_id in another list
			size_t block_id = (*p)[i];
			it = find( (*q).begin(), (*q).end(), block_id);

			// block_id not in another list, net is boundary
			if( it == (*q).end() ){
				Block & blk = block_info[block_id];
				blk.boundary_netlist.push_back(net);
				if( !nk->isX() ) {
					// stamp value into block_ids
					size_t k1 = nk->id_in_block[i];
					Matrix &pk = A[block_id];	
					pk.push_back(k1,k1, G);
				}
			}
			// else 2 nodes belongs to the same block
			// stamp resistor
			else if( !nk->isX() ) {
				size_t k1 = nk->id_in_block[i];
				Matrix &pk = A[block_id];

				size_t j1 = it - (*q).begin();
				size_t l1 = nl->id_in_block[j1];

				pk.push_back(k1,k1, G);
				if(!nl->isX() && l1 < k1) // only store the lower triangular part
					pk.push_back(k1,l1,-G);
			}
		}// end of for k
	}// end of for j	
}

void Circuit::stamp_block_current(Net * net, Matrix * A){
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground() && !nk->isX() ) { 
		for(size_t i=0;i<nk->blocklist.size();i++){
			size_t block_idk = nk->blocklist[i];
			Block &block_k = block_info[block_idk];
			//Vec & pk = block_k.b;
			size_t k = nk->id_in_block[i];
			block_k.bp[k] += -net->value;
			//pk[k] += -net->value;
		}
	}
	if( !nl->is_ground() && !nl->isX() ) {
		for(size_t i=0;i<nl->blocklist.size();i++){
			size_t block_idl = nl->blocklist[i];
			Block & block_l = block_info[block_idl];
			//Vec & pl = block_l.b;
			size_t l = nl->id_in_block[i];
			block_l.bp[l] += net->value;
			//pl[l] +=  net->value;
		}
	}
}

void Circuit::stamp_block_VDD(Net * net, Matrix * A){
	// find the non-ground node
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	for(size_t i=0;i<X->rep->id_in_block.size();i++){
		size_t block_id = X->rep->blocklist[i];
		Block &block = block_info[block_id];	
		Matrix & p = A[block_id];
		//Vec & q = block.b;
		size_t id =X->rep->id_in_block[i];
		p.push_back(id, id, 1.0);
		Net * south = X->rep->nbr[SOUTH];
		if( south != NULL &&
	   	 south->type == CURRENT ){
		// this node connects to a VDD and a current
		// the current should be stamped
			//assert( feqn(1.0, q[id]) ); 
			assert( feqn(1.0, block.bp[id]) );
			block.bp[id] = net->value;
			//q[id] = net->value;	    // modify it
		}
		else{
			block.bp[id] += net->value;
			//q[id] += net->value;
		}
	}
}

// set the block_id in a node
// according to the following order
// where 0 is the original block the node should be 
// 8 1 2
// 7 0 3
// 6 5 4
//
void Circuit::set_blocklist(Node * nd){
	const double len_per_block_x = block_info.len_per_block_x;
	const double len_per_block_y = block_info.len_per_block_y;
	const double len_ovr_x = block_info.len_ovr_x;
	const double len_ovr_y = block_info.len_ovr_y;
	const size_t X_BLOCKS = block_info.X_BLOCKS;
	const size_t Y_BLOCKS = block_info.Y_BLOCKS;
	const long x = nd->pt.x - x_min;	// point relative location
	const long y = nd->pt.y - y_min;
	size_t bx0 = x / len_per_block_x;
	size_t by0 = y / len_per_block_y;
	long bx, by;		// block index
	double lx, ly, ux, uy;	// block bounding box

	const long dx[]={0, 0, 1, 1,  1,  0, -1, -1, -1};
	const long dy[]={0, 1, 1, 0, -1, -1, -1,  0, 1};

	// test whether the node is in one of the nine blocks
	for(int i=0;i<9;i++){
		bx = bx0 + dx[i];
		by = by0 + dy[i];

		// check the block index is valid
		if( bx < 0 || bx >= (long)X_BLOCKS ) continue;
		if( by < 0 || by >= (long)Y_BLOCKS ) continue;

		// compute block coordinate
		lx = bx * len_per_block_x - len_ovr_x;
		ly = by * len_per_block_y - len_ovr_y;
		ux = (bx+1) * len_per_block_x + len_ovr_x;
		uy = (by+1) * len_per_block_y + len_ovr_y;

		// check if the point is in the block
		if( !(x>=lx && x<=ux && y>=ly && y<=uy) ) continue;	

		size_t id = by * X_BLOCKS + bx;
		assert( id<X_BLOCKS* Y_BLOCKS );

		/*
		vector<size_t >::const_iterator it;
		it = find(nd->blocklist.begin(),nd->blocklist.end(), id);
		if(it!=nd->blocklist.end()){
			printf("id=%ld, i=%d\n", id,i);
			printf("xbase,ybase=%lf %lf\n", lx_base, ly_base);
			printf("lx,ly=%lf %lf\n", lx, ly);
			printf("ux,uy=%lf %lf\n", ux, uy);
			printf("x,y=%ld %ld\n", x, y);
			printf("bx,by=%ld %ld\n", bx, by);
			printf("xc,yc=%ld %ld\n", x_center, y_center);
			continue;
		}
		*/
		nd->blocklist.push_back(id);
		Block & block = block_info[id];
		nd->id_in_block.push_back(block.count++);
	}
}

void Circuit::get_parameters(
		double & epsilon,
		double & omega,
		double & overlap_ratio,
		size_t & max_block_nodes,
		int & mode){
	epsilon		= EPSILON;
	omega		= OMEGA;
	overlap_ratio	= OVERLAP_RATIO; 
	max_block_nodes	= MAX_BLOCK_NODES;
	mode		= MODE;
}

// default values of these parameters are at the begining of this file
void Circuit::set_parameters(
		double epsilon, 
		double omega, 
		double overlap_ratio,
		size_t max_block_nodes,
		int mode){
	EPSILON		= epsilon;
	OMEGA		= omega;
	OVERLAP_RATIO 	= overlap_ratio;
	MAX_BLOCK_NODES	= max_block_nodes;
	MODE		= mode;
}

// choose an appropriate omega for the circuit s.t.
// - node size (use replist)
// - type (c4 or wb)
// - number of layers
void Circuit::select_omega(){
	double omega=OMEGA;
	size_t num_nodes = replist.size();
	size_t num_layers = layers.size();
	if( num_nodes < 0.05e6 )
		omega=1.0;
	else if (num_nodes < 0.2e6 )
		omega = 1.1;
	else if (num_nodes < 0.3e6 )
		omega = 1.2;
	else if (num_nodes < 0.5e6 )
		omega = 1.3;
	else if (num_nodes < 1.2e6 )
		omega = 1.4;
	else
		omega = 1.5;

	if( circuit_type == WB && num_layers >= 8 ) omega += 0.2;

	if( circuit_type == C4 ) omega += 0.1;

	if( name == "GND" && num_nodes < 1.2e6) omega -= 0.1;

	if( omega >= 1.6 ) omega = 1.6;
	if( omega <= 1.0 ) omega = 1.0;

	OMEGA = omega;
}

// Randomly choose a number of sample nodes to monitor
void Circuit::get_samples(){
	size_t num_nodes = replist.size();
	srand(time(NULL));
	while(sample.size()<SAMPLE_NUM_NODE){
		int id = rand() % num_nodes;
		sample.push_back(replist[id]);
	}
}

bool Circuit::check_diverge() const{
	for(size_t i=0;i<SAMPLE_NUM_NODE;i++){
		double x = sample[i]->value;
		if(VDD > 0){
			if( x < 0.0 || x > VDD ) return true;
		}
		else
			if(x<0.0) return true;
	}
	return false;
}

Node * Circuit::merge_along_dir_one_pass(Node * start, DIRECTION dir, bool remove){
	double sum = 0.0;
	DIRECTION ops = get_opposite_dir(dir);
	Node * p = start;

	// traverse along the direction, sum the resistor value and set the node end
	while(1){
		p = p->get_nbr_node(dir);
		p->end[ops] = start;
		Net * net = p->nbr[ops];
		sum += net->value;
		p->eqvr[ops] = sum;
		if( remove ) {
			size_t id = net_id[net];
			net_set[RESISTOR][id] = NULL;
			delete net;
		}
		if( !p->is_mergeable() ) break;
	}

	return p;	// return end point
}

// merge a line along direction
void Circuit::merge_along_dir(Node * node, DIRECTION dir){
	// two pass traversal
	DIRECTION ops = get_opposite_dir(dir);
	node->end[dir] = merge_along_dir_one_pass(node, dir, false);
	Node * other = node->end[dir];
	other->end[ops] = node;
	merge_along_dir_one_pass(other, ops, true);
	//assert( ret == node );

	// add a new net between `node' and its end
	Net * net = new Net(RESISTOR, node->eqvr[dir], node, other);
	node->nbr[dir] = other->nbr[ops] = net;
	this->add_net(net);
}
