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
size_t Circuit::MAX_BLOCK_NODES =100000;//5500;
double Circuit::OMEGA = 1.2;
double Circuit::OVERLAP_RATIO = 0;
int    Circuit::MODE = 0;
const int MAX_ITERATION =  1000;
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

	// mpi relate	
	bd_x_g=NULL;
	internal_x_g = NULL;
	
	bd_x = NULL;
	internal_x = NULL;
	bd_base = NULL;
	internal_base = NULL;
	bd_base_gd = NULL;
	internal_base_gd = NULL;
	bd_base_g = NULL;
	internal_base_g = NULL;
	bd_size_g = NULL;
	internal_size_g = NULL;
	bd_dd_size = NULL;
	internal_dd_size = NULL;
	bd_dd_size_g = NULL;
	internal_dd_size_g = NULL;

	bd_size = 0;
	internal_size = 0;
	block_size=0;
}

// Trick: do not release memory to increase runtime
Circuit::~Circuit(){
	for(size_t i=0;i<nodelist.size();i++) delete nodelist[i];
	bd_nodelist_sw.clear();
	bd_nodelist_s.clear();
	bd_nodelist_se.clear();
	bd_nodelist_w.clear();
	bd_nodelist_e.clear();
	bd_nodelist_nw.clear();
	bd_nodelist_n.clear();
	bd_nodelist_ne.clear();

	internal_nodelist_sw.clear();
	internal_nodelist_s.clear();
	internal_nodelist_se.clear();
	internal_nodelist_w.clear();
	internal_nodelist_e.clear();
	internal_nodelist_nw.clear();
	internal_nodelist_n.clear();
	internal_nodelist_ne.clear();
	
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		for(it=ns.begin();it!=ns.end();++it)
			delete *it;
	}
	
	// mpi related variables
	delete [] bd_x_g;
	delete [] internal_x_g;
	delete [] bd_x;
	delete [] internal_x;
	delete [] bd_base;
	delete [] internal_base;
	delete [] bd_base_gd;
	delete [] internal_base_gd;
	delete [] bd_base_g;
	delete [] internal_base_g;
	delete [] bd_size_g;
	delete [] internal_size_g;
	delete [] bd_dd_size;
	delete [] internal_dd_size;
	delete [] bd_dd_size_g;
	delete [] internal_dd_size_g;
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
// sort nodelist
void Circuit::sort_nodes(){
	sort(nodelist.begin(), nodelist.end(), compare_node_ptr);
	// update node id mapping, 
	// NOTE: ground node will be the last
}

// sort 4 boudnary nodelist
void Circuit::sort_bd_nodes(int &my_id){
	sort(bd_nodelist_sw.begin(), bd_nodelist_sw.end(), 
			compare_node_ptr);

	sort(bd_nodelist_s.begin(), bd_nodelist_s.end(), 
			compare_node_ptr);
	
	sort(bd_nodelist_se.begin(), bd_nodelist_se.end(), 
			compare_node_ptr);

	sort(bd_nodelist_w.begin(), bd_nodelist_w.end(), 
			compare_node_ptr);

	sort(bd_nodelist_e.begin(), bd_nodelist_e.end(), 
			compare_node_ptr);

	sort(bd_nodelist_nw.begin(), bd_nodelist_nw.end(), 
			compare_node_ptr);
	
	sort(bd_nodelist_n.begin(), bd_nodelist_n.end(), 
			compare_node_ptr);

	sort(bd_nodelist_ne.begin(), bd_nodelist_ne.end(), 
			compare_node_ptr);
}

// sort 4 boudnary nodelist
void Circuit::sort_internal_nodes(int &my_id){
	sort(internal_nodelist_sw.begin(), 
		internal_nodelist_sw.end(), compare_node_ptr);

	sort(internal_nodelist_s.begin(), 
		internal_nodelist_s.end(), compare_node_ptr);
	
	sort(internal_nodelist_se.begin(), 
		internal_nodelist_se.end(), compare_node_ptr);

	sort(internal_nodelist_w.begin(), 
		internal_nodelist_w.end(), compare_node_ptr);
	
	sort(internal_nodelist_e.begin(), 
		internal_nodelist_e.end(), compare_node_ptr);

	sort(internal_nodelist_nw.begin(), 
		internal_nodelist_nw.end(), compare_node_ptr);
	
	sort(internal_nodelist_n.begin(), 
		internal_nodelist_n.end(), compare_node_ptr);

	sort(internal_nodelist_ne.begin(), 
		internal_nodelist_ne.end(), compare_node_ptr);
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
void Circuit::solve_init(int &my_id){
	sort_nodes();
	sort_bd_nodes(my_id);
	sort_internal_nodes(my_id);
	
	size_t size = nodelist.size() - 1;
	Node * p = NULL;
	Net * net = NULL;
	size_t nr = 0;
	size_t i=0;
	for(i=0, nr=0;i<size;i++){
		p=nodelist[i];
		net = p->nbr[TOP];
		
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
	
	block_info.count = nr;

	size_t n_merge = mergelist.size();
	size_t n_nodes = nodelist.size();
	size_t n_reps  = replist.size();
	double ratio = n_merge / (double) (n_merge + n_reps);
	
	net_id.clear();
}

// build up block info
// 1. Find block divide point
// 2. Set each node in replist into block
// 3. Compute block size
// 4. Insert boundary netlist into map
void Circuit::block_init(int &my_id, Matrix &A, MPI_CLASS &mpi_class){
	block_info.update_block_geometry(mpi_class);
	//block_info.allocate_resource(cm);
	block_info.allocate_resource_LU();
	copy_node_voltages_block();
	stamp_block_matrix(my_id, A, mpi_class);
}

// stamp the nets by sets, block version
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_block_matrix(int &my_id, Matrix &A, MPI_CLASS &mpi_class){	
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
				stamp_block_resistor(my_id, *it, A);
			}
			break;
		case CURRENT:
			for(it=ns.begin();it!=ns.end();++it){
				stamp_block_current(my_id, (*it), mpi_class);
			}
			break;
		case VOLTAGE:
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_block_VDD(my_id,(*it), A);
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	make_A_symmetric(block_info.b);
	
	A.set_row(block_info.count);
	if(block_info.count >0){
		//block_info.CK_decomp(A, cm);
		block_info.get_col_compressed(A);
		block_info.LU_decomposition();
	}
}

// 1. mark rep nodes into corresponding blocks
// 2. find block size of replist nodes
// 3. allocate rhs size of each block
// 4. find local block index for each node
void Circuit::find_block_size(MPI_CLASS &mpi_class){	
	// for each block, allocate resource
	block_info.allocate_resource_LU();
}

void Circuit::solve(int &orig_rank, int &my_id, int&num_procs, MPI_Comm &comm, MPI_CLASS &mpi_class){
	// each block is solved by IT
	solve_IT(orig_rank, my_id, num_procs, mpi_class, comm);
}

// solve Circuit
bool Circuit::solve_IT(int &orig_rank, int &my_id, int&num_procs, MPI_CLASS &mpi_class, MPI_Comm &comm){	
	double time=0;
	double t1, t2;
	
	total_blocks = mpi_class.num_blocks;

	// did not find any `X' node
	if( circuit_type == UNKNOWN )
		circuit_type = C4;
	if(my_id<num_procs){
		solve_init(my_id);
	}
	
	block_init(my_id, A, mpi_class);
	boundary_init(my_id, num_procs, comm);
	internal_init(my_id, num_procs, comm);
	// stores 4 boundary base into bd_base_gd
	MPI_Gather(bd_base, 8, MPI_INT, bd_base_gd, 8, MPI_INT,
		0, comm);
	
	MPI_Gather(internal_base, 8, MPI_INT, internal_base_gd,
		8, MPI_INT, 0, comm);
		
	int iter = 0;	
	double diff=0;
	bool successful = false;

	// before iteration, copy boundary nodes value to corresponding blocks
	assign_bd_internal_array(my_id);
	MPI_Gatherv(internal_x, internal_size, MPI_FLOAT, 
		internal_x_g, internal_size_g, 
		internal_base_g, MPI_FLOAT, 0, comm);
	
	// reorder boundary array according to nbrs
	if(my_id==0)	reorder_bd_x_g(mpi_class);
	time=0;
	t1= MPI_Wtime();
	while( iter < MAX_ITERATION ){
		diff = solve_iteration(my_id, iter, num_procs, mpi_class, comm);
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
	if(my_id==0 && mpi_class.color ==0) clog<<"iter time: "<<time<<endl;
	if(my_id==0 && mpi_class.color ==0) cout<<replist<<endl;
	if(my_id==0){
		//clog<<"solve iteration time is: "<<time<<endl;
		clog<<"# iter: "<<iter<<endl;
		//clog<<replist<<endl;
		get_voltages_from_block_LU_sol();
		//get_vol_mergelist();
	}

	/*if(block_info.count > 0)
		block_info.free_block_cholmod(cm);
	cholmod_finish(cm);*/
	return successful;
}

// solve blocks with mpi: multi-core
// One iteration during solving the circuit, for any block B:
// 1. update the righthand-side of the matrix of B
// 2. solve the matrix
// 3. update node voltages
// 4. track the maximum error of solution
double Circuit::solve_iteration(int &my_id, int &iter,
		int&num_procs, MPI_CLASS &mpi_class, MPI_Comm &comm){	
	float diff = .0;
	float diff_root=0;

	// 0 rank cpu will scatter all bd valuesfrom bd_x_g to bd_x
	MPI_Scatterv(bd_x_g, bd_size_g, 
			bd_base_g, MPI_FLOAT, bd_x, bd_size, 
			MPI_FLOAT, 0, comm);
	
	assign_bd_array();

	// new rhs store in bnewp or bnew
	block_info.update_rhs(my_id);
	
	// x_old stores old solution
	for(size_t j=0;j<block_info.count;j++)
		block_info.xold[j] = block_info.x[j];	

	if(block_info.count>0){
		block_info.solve_LU();
		//block_info.solve_CK(cm);
		//block_info.xp = static_cast<double *>(block_info.x_ck->x);
	}
	//if(my_id==0) clog<<block_info.x<<endl;
	diff = modify_voltage(my_id, block_info, 
			block_info.xold);

	assign_bd_internal_array(my_id);
	// 0 rank cpu will gather all the solution from bd_x
	// to bd_x_g
	MPI_Gatherv(internal_x, internal_size, MPI_FLOAT, 
		internal_x_g, internal_size_g, 
		internal_base_g, MPI_FLOAT, 0, comm);
	
	// reorder boundary array according to nbrs
	if(my_id==0){
		reorder_bd_x_g(mpi_class);
	}
	
	MPI_Reduce(&diff, &diff_root, 1, MPI_FLOAT, MPI_MAX, 0, comm);
	MPI_Bcast(&diff_root, 1, MPI_FLOAT, 0, comm);	
	
	//if(my_id==0) clog<<"iter, diff: "<<iter<<" "<<diff_root<<endl;
	return diff_root;
}

double Circuit::modify_voltage(int &my_id, Block &block, Vec & x_old){
	double max_diff = 0.0;
	//if(get_name()=="VDDA") OMEGA = 1.0;
	//else OMEGA = 1.15;
	OMEGA = 1.0;
	for(size_t i=0;i<block.count;i++){
		block.x[i] = (1-OMEGA)*x_old[i] + OMEGA*
			block.x[i];
		// update block nodes value
		block.nodes[i]->value = block.x[i];
		double diff = fabs(x_old[i] - block.x[i]);
		if( diff > max_diff ) max_diff = diff;
	}

	return max_diff;
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
/*void Circuit::get_vol_mergelist(){
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
}*/

// copy solution of block into circuit
void Circuit::get_voltages_from_block_LU_sol(){
	for(size_t i=0;i<nodelist.size()-1;i++){
		Node * node = nodelist[i];
		//if( node->is_mergeable() ) continue;
		size_t id = node->rep->rid;
		double v = block_info.x[id];
		node->value = v;
	}
}

// 1. copy node voltages from the circuit to a Vec
//    from = true then copy circuit to x
//    else copy from x to circuit
// 2. map block voltage into global_index
void Circuit::copy_node_voltages_block(){
	size_t id;
	// copy node voltages from nodelist
	for(size_t i=0;i<replist.size();i++){
		Node *node = replist[i];
		id = node->rid;
		block_info.x[id] = replist[i]->value;
		block_info.nodes[id] = replist[i];
	}
}

void Circuit::make_A_symmetric(Vec &b){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p, *q;

	for(it=ns.begin();it!=ns.end();it++){
		if( (*it) == NULL || (*it)->flag_bd ==true) 
			continue;
		assert( fzero((*it)->value) == false );
		Node *nd[] = {(*it)->ab[0]->rep, (*it)->ab[1]->rep};
		// node a points to X node
		if(nd[0]->isX()){
			p = nd[0]; q = nd[1];
		}
		else if(nd[1]->isX()){
			p = nd[1]; q = nd[0];
		}
		else continue;
		size_t id = q->rid;
		double G = 1.0 / (*it)->value;
		b[id] += p->value * G;
	}
}

// =========== stamp block version of matrix =======

void Circuit::stamp_block_resistor(int &my_id, Net * net, Matrix &A){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
	
	double G;	
	G = 1./net->value;
	
	if(net->flag_bd ==true)
		block_info.boundary_netlist.push_back(net);

	for(size_t j=0;j<2;j++){
		Node *nk = nd[j], *nl = nd[1-j];
		// if boundary net
		if(net->flag_bd ==true){
			if(nk->flag_bd ==false && !nk->isX()){
				// stamp value into block_ids
				size_t k1 = nk->rid;
				A.push_back(k1,k1, G);
			}
		}
		// else internal net
		else if( !nk->isX() ) {
			size_t k1 = nk->rid;
			size_t l1 = nl->rid;
			A.push_back(k1,k1, G);
			if(!nl->isX() )//&& l1 < k1) // only store the lower triangular part
				A.push_back(k1,l1,-G);
		}
	}// end of for j	
}

void Circuit::stamp_block_current(int &my_id, Net * net, MPI_CLASS &mpi_class){
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;

	// only stamp for internal node
	if( !nk->is_ground() && !nk->isX() && nk->flag_bd ==false) {
		size_t k = nk->rid;
		block_info.b[k] += -net->value;
		//pk[k] += -net->value;
	}
	if( !nl->is_ground() && !nl->isX() && nl->flag_bd ==false) {
		size_t l = nl->rid;
		block_info.b[l] += net->value;
		//pl[l] +=  net->value;
	}
}

void Circuit::stamp_block_VDD(int &my_id, Net * net, Matrix &A){
	// find the non-ground node
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];

	if(X->rep->flag_bd ==true) return;
	// do stamping for internal node
	long id =X->rep->rid;
	A.push_back(id, id, 1.0);
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
			south->type == CURRENT ){
		// this node connects to a VDD and a current
		// the current should be stamped
		//assert( feqn(1.0, q[id]) ); 
		assert( feqn(1.0, block_info.b[id]) );
		block_info.b[id] = net->value;
		//q[id] = net->value;	    // modify it
	}
	else{
		block_info.b[id] += net->value;
		//q[id] += net->value;
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

/*Node * Circuit::merge_along_dir_one_pass(Node * start, DIRECTION dir, bool remove){
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
	//clog<<"newly added net: "<<*net<<endl;
	this->add_net(net);
}*/

void Circuit::boundary_init(int &my_id, int &num_procs, MPI_Comm &comm){
	assign_bd_base(my_id);
	assign_bd_dd_size(my_id);

	// allocate boundary_nodelist size
	bd_dd_size_g = new int [8*num_procs];	
	
	MPI_Gather(bd_dd_size, 8, MPI_INT, bd_dd_size_g, 
			8, MPI_INT, 0, comm);
	
	bd_x = new float[bd_size];

	// assign bd_size_g value
	bd_size_g = new int[num_procs];
	MPI_Gather(&bd_size, 1, MPI_INT, bd_size_g, 1, MPI_INT,
			0, comm);

	total_size = 0;
	if(my_id ==0){
		for(int i=0;i<total_blocks;i++)
			total_size += bd_size_g[i];
	}

	bd_x_g = new float[total_size];
	bd_base_g = new int [num_procs];
	//if(my_id != 0) return;
	int base = 0;
	for(int i=0;i<num_procs;i++){
		bd_base_g[i] = base;
		base += bd_size_g[i];
	}
	bd_base_gd = new int [8*num_procs];
	for(int i=0;i<8*num_procs;i++){
		bd_base_gd[i] = 0;
	}
}

void Circuit::assign_bd_base(int &my_id){
	// boundary base along e, w, n, s directions
	bd_base = new int [8];

	// assign 4 base values for array bd_x
	int  base = 0;
	bd_base[0] = base;
	base += bd_nodelist_sw.size();
	bd_base[1] = base;
	base += bd_nodelist_s.size();
	bd_base[2] = base;
	base += bd_nodelist_se.size();
	bd_base[3] = base;
	base += bd_nodelist_w.size();
	bd_base[4] = base;
	base += bd_nodelist_e.size();
	bd_base[5] = base;
	base += bd_nodelist_nw.size();
	bd_base[6] = base;
	base += bd_nodelist_n.size();
	bd_base[7] = base;
}

void Circuit::assign_bd_dd_size(int &my_id){
	bd_dd_size = new int [8];

	bd_dd_size[0] = bd_nodelist_sw.size();
	bd_dd_size[1] = bd_nodelist_s.size();
	bd_dd_size[2] = bd_nodelist_se.size();
	bd_dd_size[3] = bd_nodelist_w.size();	
	bd_dd_size[4] = bd_nodelist_e.size();
	bd_dd_size[5] = bd_nodelist_nw.size();
	bd_dd_size[6] = bd_nodelist_n.size();
	bd_dd_size[7] = bd_nodelist_ne.size();

	bd_size = 0;
	for(int i=0;i<8;i++)
		bd_size += bd_dd_size[i];
}

void Circuit::assign_internal_base(int &my_id){
	// boundary base along e, w, n, s directions
	internal_base = new int [8];

	// assign 4 base values for array bd_x
	int  base = 0;
	internal_base[0] = base;
	base += internal_nodelist_sw.size();
	internal_base[1] = base;
	base += internal_nodelist_s.size();
	internal_base[2] = base;
	base += internal_nodelist_se.size();
	internal_base[3] = base;
	base += internal_nodelist_w.size();
	internal_base[4] = base;
	base += internal_nodelist_e.size();
	internal_base[5] = base;
	base += internal_nodelist_nw.size();
	internal_base[6] = base;
	base += internal_nodelist_n.size();
	internal_base[7] = base;
}

// assign bd internal arrays
void Circuit::assign_bd_internal_array(int &my_id){
	//clog<<"sw: "<<endl;
	assign_bd_internal_array_dir(internal_base[0], 
		internal_nodelist_sw, internal_x);
	//clog<<"s: "<<endl;
	assign_bd_internal_array_dir(internal_base[1], 
		internal_nodelist_s, internal_x);
	//clog<<"se: "<<endl;
	assign_bd_internal_array_dir(internal_base[2], 
		internal_nodelist_se, internal_x);
	//clog<<"w: "<<endl;
	assign_bd_internal_array_dir(internal_base[3], 
		internal_nodelist_w, internal_x);
	//clog<<"e: "<<endl;
	assign_bd_internal_array_dir(internal_base[4], 
		internal_nodelist_e, internal_x);
	//clog<<"nw: "<<endl;
	assign_bd_internal_array_dir(internal_base[5], 
		internal_nodelist_nw, internal_x);
	//clog<<"n: "<<endl;
	assign_bd_internal_array_dir(internal_base[6], 
		internal_nodelist_n, internal_x);
	//clog<<"ne: "<<endl;
	assign_bd_internal_array_dir(internal_base[7], 
		internal_nodelist_ne, internal_x);
}

// assign 4 boundary internal nodes value, store
// them in array bd_x
void Circuit::assign_bd_internal_array_dir(int &base, NodePtrVector & list, float *internal_x){
	Node *nd;
	for(size_t i=0;i<list.size();i++){
		nd = list[i]->rep;
		internal_x[base+i] = nd->value;
		//clog<<base+i<<" "<<*nd<<endl;
	}
}

// assign 4 boundary nodes value from receiving
// array bd_x
void Circuit::assign_bd_array(){
	int base =0;
	base = bd_base[0];
	assign_bd_array_dir(base, bd_nodelist_sw);
	
	base = bd_base[1];
	assign_bd_array_dir(base, bd_nodelist_s);	
	
	base = bd_base[2];
	assign_bd_array_dir(base, bd_nodelist_se);	

	base = bd_base[3];
	assign_bd_array_dir(base, bd_nodelist_w);

	base = bd_base[4];
	assign_bd_array_dir(base, bd_nodelist_e);

	base = bd_base[5];
	assign_bd_array_dir(base, bd_nodelist_nw);

	base = bd_base[6];
	assign_bd_array_dir(base, bd_nodelist_n);

	base = bd_base[7];
	assign_bd_array_dir(base, bd_nodelist_ne);
}

// input: internal_x
// output: bd_x
// func: assign internal_x to corresponding bd_x
void Circuit::reorder_bd_x_g(MPI_CLASS &mpi_class){
	// total_blocks to reorder
	// each block boundary is with e, w, n, s bd
	int X_BLOCKS = mpi_class.X_BLOCKS;
	int Y_BLOCKS = mpi_class.Y_BLOCKS;

	int base_p = 0; // stores current block base
	int base_q = 0; // stores nbr block base
	int base_glo_p = 0;
	int base_glo_q = 0;
	int by = 0;
	int bx = 0;
	int bid_nbr = 0;
	int size = 0;

	for(int i=0;i<total_blocks;i++){
		by = i / X_BLOCKS;
		bx = i % X_BLOCKS;
		//clog<<"reorder for "<<i<<" th block."<<endl<<endl;
		// compute block base for i
		base_glo_p = bd_base_g[i];

		// south west nbr dir
		base_p = base_glo_p+bd_base_gd[8*i];
		if(by>=1 && bx>=1){
			bid_nbr = i - X_BLOCKS - 1;
			//clog<<"east bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+7];
			size = bd_dd_size_g[8*i];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}
		
		// south nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+1];
		if(by>=1){
			bid_nbr = i - X_BLOCKS;
			//clog<<"south bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+6];
			size = bd_dd_size_g[8*i+1];
			//clog<<size<<endl;
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}

		// south east nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+2];
		if(by>=1 && bx<X_BLOCKS-1){
			bid_nbr = i - X_BLOCKS + 1;
			//clog<<"east bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+5];
			size = bd_dd_size_g[8*i+2];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}

		// west nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+3];
		if(bx>=1){
			bid_nbr = i - 1;
			//clog<<"west bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+4];
			size = bd_dd_size_g[8*i+3];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}

		// east nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+4];
		//clog<<"base_glo_p: "<<base_glo_p<<" base_p: "<<base_p<<endl;
		if(bx<X_BLOCKS-1){
			bid_nbr = i+1;
			//clog<<"east bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			//clog<<"base_glo_q: "<<base_glo_q<<" ";
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+3];
			//clog<<"base_q: "<<base_q<<endl;
			size = bd_dd_size_g[8*i+4];
			//clog<<"size: "<<size<<endl;
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}
		
		// north west nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+5];
		if(by<Y_BLOCKS-1 && bx>=1){
			bid_nbr = i + X_BLOCKS - 1;
			//clog<<"east bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+2];
			size = bd_dd_size_g[8*i+5];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}

		// north nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+6];
		if(by<Y_BLOCKS-1){
			bid_nbr = i + X_BLOCKS;
			//clog<<"north bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr+1];
			size = bd_dd_size_g[8*i+6];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}

		// north east nbr dir
		base_p = base_glo_p+bd_base_gd[8*i+7];
		if(by<Y_BLOCKS-1 && bx<X_BLOCKS-1){
			bid_nbr = i + X_BLOCKS + 1;
			//clog<<"east bid: "<<bid_nbr<<endl;
			// compute block base
			base_glo_q = internal_base_g[bid_nbr];
			// find east nbr block base
			base_q = base_glo_q+internal_base_gd[8*bid_nbr];
			size = bd_dd_size_g[8*i+7];
			for(int j=0;j<size;j++){
				bd_x_g[base_p+j] = internal_x_g[base_q+j];
				//clog<<base_p+j<<" "<<base_q+j<<" "<<internal_x_g[base_q+j]<<endl;
			}
		}
	}
}

void Circuit::internal_init(int &my_id, int &num_procs, MPI_Comm &comm){
	assign_internal_base(my_id);

	// allocate boundary_nodelist size
	internal_dd_size = new int [8];
	internal_dd_size[0] = internal_nodelist_sw.size();
	internal_dd_size[1] = internal_nodelist_s.size();
	internal_dd_size[2] = internal_nodelist_se.size();
	internal_dd_size[3] = internal_nodelist_w.size();	
	internal_dd_size[4] = internal_nodelist_e.size();
	internal_dd_size[5] = internal_nodelist_nw.size();
	internal_dd_size[6] = internal_nodelist_n.size();
	internal_dd_size[7] = internal_nodelist_ne.size();

	internal_dd_size_g = new int [8*num_procs];
	
	MPI_Gather(internal_dd_size, 8, MPI_INT, 
		internal_dd_size_g, 8, MPI_INT, 0, comm);

	internal_size = internal_nodelist_sw.size()+
		  internal_nodelist_s.size()+
		  internal_nodelist_se.size()+
		  internal_nodelist_w.size()+
		  internal_nodelist_e.size()+
		  internal_nodelist_nw.size()+
		  internal_nodelist_n.size()+
		  internal_nodelist_ne.size();

	internal_x = new float[internal_size];

	// assign bd_size_g value
	internal_size_g = new int[num_procs];
	MPI_Gather(&internal_size, 1, MPI_INT, 
			internal_size_g, 1, MPI_INT,
			0, comm);
	
	total_internal_size = 0;
	if(my_id ==0){
		for(int i=0;i<total_blocks;i++)
			total_internal_size += internal_size_g[i];
	}

	internal_x_g = new float[total_internal_size];
	internal_base_g = new int [num_procs];
	//if(my_id != 0) return;
	int base = 0;
	for(int i=0;i<num_procs;i++){
		internal_base_g[i] = base;
		base += internal_size_g[i];
	}
	internal_base_gd = new int [8*num_procs];
	for(int i=0;i<8*num_procs;i++){
		internal_base_gd[i] = 0;
	}
}

// assign dir boundary node value
void Circuit::assign_bd_array_dir(int &base, NodePtrVector &list){
	Node *p;
	for(size_t i=0;i<list.size();i++){
		p = list[i]->rep;
		p->value = bd_x[base+i];
	}
}
