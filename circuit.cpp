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
//#include "umfpack.h"
#include "circuit.h"
#include "util.h"
#include "algebra.h"
#include "node.h"

#include "mpi.h"
using namespace std;

double Circuit::EPSILON = 1e-3;
size_t Circuit::MAX_BLOCK_NODES =100000;//5500;
double Circuit::OMEGA = 1.2;
double Circuit::OVERLAP_RATIO = 0;
int    Circuit::MODE = 0;
const int MAX_ITERATION = 1000;
const int SAMPLE_INTERVAL = 5;
const size_t SAMPLE_NUM_NODE = 10;
const double MERGE_RATIO = 0.3;
int Circuit::NUM_BLOCKS_X = 2;
int Circuit::NUM_BLOCKS_Y = 1;
int Circuit::DEBUG=1;

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
	// delete node
	for(size_t i=0;i<nodelist.size();i++) delete nodelist[i];
	// delete node
	for(size_t i=0;i<replist.size();i++) delete replist[i];
	// delete node
	for(size_t i=0;i<mergelist.size();i++) delete mergelist[i];
	nodelist.clear();
	replist.clear();
	mergelist.clear();
	block_vec.clear();
	for(size_t i=0;i<bd_netlist.size();i++)
    		delete bd_netlist[i];
    	bd_netlist.clear();

	// delete bd nodes
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
	
	// delete nets
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		for(it=ns.begin();it!=ns.end();++it)
			delete *it;
	}

	// free Li, Lx and so on
	/*delete [] Lx;
	delete [] Lp;
	delete [] Lnz;
	delete [] Li;*/
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
	//clog<<"* UF_long size = "<< sizeof(UF_long)<<endl;
	clog<<"* Max nodelist = "<<(size_t)nodelist.max_size()<<endl;
	clog<<"****            END              ****"<<endl<<endl;
}

// functor to be used in STL sort function
// order: y > x > z > flag 
// input: two node a, b
// return true if a < b, false o/w
// note that ground note are put to last
#if DEBUG
bool compare_node_ptr(const Node * a, const Node * b){
	if( a->is_ground() ) return false;
	if (b->is_ground() ) return true;

	if( a->pt.y == b->pt.y ){
		if( a->pt.x == b->pt.x ){
			if( a->pt.z == b->pt.z ){
				return (a->isS() > b->isS());
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
#endif
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

void Circuit::print_matlab(Matrix A){
	// uncomment this if want to output to a file
	//freopen("output.txt","w",stdout);

	// don't output ground node
	for(size_t i=0;i<A.size();i++){
		printf("%ld %ld %.5e\n", A.Ti[i]+1, A.Tj[i]+1, A.Tx[i]);
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
	if(my_id==0)
		clog<<"total num_nodes: "<<nodelist.size()<<endl;
	sort_bd_nodes(my_id);
	sort_internal_nodes(my_id);
	
	/*if(my_id==3){
		clog<<"sw: "<<internal_nodelist_sw<<endl;
		clog<<"s: "<<internal_nodelist_s<<endl;
		clog<<"se: "<<internal_nodelist_se<<endl;
		clog<<"w: "<<internal_nodelist_w<<endl;
		clog<<"e: "<<internal_nodelist_e<<endl;
		clog<<"nw: "<<internal_nodelist_nw<<endl;
		clog<<"n: "<<internal_nodelist_n<<endl;
		clog<<"ne: "<<internal_nodelist_ne<<endl;
		clog<<endl;
	}

	if(my_id==3){
		clog<<"sw: "<<bd_nodelist_sw<<endl;
		clog<<"s: "<<bd_nodelist_s<<endl;
		clog<<"se: "<<bd_nodelist_se<<endl;
		clog<<"w: "<<bd_nodelist_w<<endl;
		clog<<"e: "<<bd_nodelist_e<<endl;
		clog<<"nw: "<<bd_nodelist_nw<<endl;
		clog<<"n: "<<bd_nodelist_n<<endl;
		clog<<"ne: "<<bd_nodelist_ne<<endl;
	}*/

	size_t size = nodelist.size() - 1;
	Node * p = NULL;
	size_t nr = 0;
	size_t i=0;
	for(i=0, nr=0;i<size;i++){
		p=nodelist[i];

		// test if it can be merged
		/*if( p->is_mergeable() ){
			mergelist.push_back(p);
			continue;
		}*/
		
		Net * net = p->nbr[TOP];
		//merge_node(p);
		
		// find the VDD value
		if( p->isS()==Y ) VDD = p->get_value();
		
		// test short circuit
		if( p->isS() !=Y && // Y must be representative 
		    net != NULL &&
		    fzero(net->value) ){
			// TODO: ensure ab[1] is not p itself
			assert( net->ab[1] != p );
			p->rep = net->ab[1]->rep;
		} // else the representative is itself

		// push the representatives into list
		if( p->rep == p ) {
			replist.push_back(p);
			//rep_id[p] = nr; // save the id
			p->rid = nr;
			++nr;
		}	
	}// end of for i

	/*if(my_id==0){
		cout<<endl;
		for(int i=0;i<nodelist.size()-1;i++)
			cout <<*nodelist[i]<<" "<<nodelist[i]->rid<<endl;
		//clog<<endl;
	}*/
	

	// if(nr >=0)
		// block_info.count = nr;

	size_t n_merge = mergelist.size();
	size_t n_nodes = nodelist.size();
	size_t n_reps  = replist.size();
	double ratio = n_merge / (double) (n_merge + n_reps);
	
	/*clog<<"mergeable  "<<n_merge<<endl;
	clog<<"replist    "<<n_reps <<endl;
	clog<<"nodelist   "<<n_nodes<<endl;
	clog<<"ratio =    "<<ratio  <<endl;*/

	net_id.clear();
	//clog<<my_id<<" "<<block_info.count<<endl;
}

// build up block info
// 1. Find block divide point
// 2. Set each node in replist into block
// 3. Compute block size
// 4. Insert boundary netlist into map
void Circuit::block_init(int &my_id, MPI_CLASS &mpi_class){
	// assign nodes into blocks and sort
	assign_block_nodes(my_id);
	assign_block_nets(my_id);
	for(size_t i=0;i<block_vec.size();i++){
		block_vec[i]->allocate_resource();
		// if(my_id==0)
			// clog<<endl<<" block: "<<i<<endl;
		// copy_node_voltages_block();
		block_vec[i]->stamp_matrix(my_id, mpi_class);
	}
}
#if DEBUG
// stamp the nets by sets, block version
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_block_matrix(int &my_id, Matrix &A, MPI_CLASS &mpi_class){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			if(my_id==0)
				clog<<"resis net. "<<ns.size()<<endl;
			for(it=ns.begin();it!=ns.end();++it){
				Net * net = *it;
				if( net == NULL ) continue;
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
			if(my_id==0)
				clog<<"VDD net: "<<ns.size()<<endl;
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_block_VDD(my_id,(*it), A);
			}
			break;
		case CAPACITANCE:
			break;
		case INDUCTANCE:
			if(my_id==0)
				clog<<"induc net: "<<ns.size()<<endl;
			for(size_t i=0;i<ns.size();i++)
				stamp_inductance_dc(A, ns[i], my_id);
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	//if(my_id==0)
		//clog<<A<<endl;
	/*if(my_id==0){
		for(int i=0;i<10;i++)
			clog<<"b origin: "<<i<<" "<<block_info.bp[i]<<endl;
	}*/
	 make_A_symmetric(block_info.bp, my_id);
	
	A.set_row(block_info.count);
	//if(my_id==0)
		//check_matrix(A);
		//cout<<"before CK_decomp. "<<endl;
	if(block_info.count >0){
		block_info.CK_decomp(A, cm);
		//if(cm->status ==1)
			//clog<<" non SPD: "<<my_id<<" "<<cm->status<<endl;
	}
	//if(my_id==0)
		//clog<<"after CK_decomp. "<<endl;
}
#endif
// stamp the nets by sets, block version
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_block_matrix_tr(int &my_id, Matrix &A, MPI_CLASS &mpi_class, Tran &tran){	
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			for(it=ns.begin();it!=ns.end();++it){
				Net * net = *it;
				if( net == NULL ) continue;
				assert( fzero(net->value) == false );
				stamp_block_resistor_tr(my_id, *it, A);
			}
			break;
		case CURRENT:
			break;
		case VOLTAGE:
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_block_VDD_tr(my_id,(*it), A);
			}
			break;
		case CAPACITANCE:
			for(size_t i=0;i<ns.size();i++)
				stamp_capacitance_tr(A, ns[i], tran, my_id);
			break;
		case INDUCTANCE:
			for(size_t i=0;i<ns.size();i++)
				stamp_inductance_tr(A, ns[i], tran, my_id);
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	//if(my_id==0)
		//clog<<A<<endl;
	/*if(my_id==0){
		for(int i=0;i<10;i++)
			clog<<"b origin: "<<i<<" "<<block_info.bp[i]<<endl;
	}*/	
}

void Circuit::solve(int &my_id, int&num_procs, MPI_CLASS &mpi_class, Tran &tran){
	// each block is solved by IT
	solve_IT(my_id, num_procs, mpi_class, tran);
	//clog<<my_id<<" finish solve: "<<endl;
}

// solve Circuit
bool Circuit::solve_IT(int &my_id, int&num_procs, MPI_CLASS &mpi_class, Tran &tran){
	double time=0;
	double t1, t2;
	
	/*cm = &c;
	cholmod_start(cm);
	cm->print = 5;*/

	total_blocks = mpi_class.X_BLOCKS *mpi_class.Y_BLOCKS;

	// did not find any `X' node
	if( circuit_type == UNKNOWN )
		circuit_type = C4;
	
	if(mpi_class.block_size>0){
		solve_init(my_id);
	}
	// update bd info of circuit and block vec
	update_geometry(my_id, mpi_class);
	// build boundary netlist of circuit
	build_bd_netlist();
	block_init(my_id, mpi_class);
	// clog<<"block_init. "<<endl;
	//return true;

	boundary_init(my_id, num_procs);

	internal_init(my_id, num_procs);
	
	bool successful = false;
	/*if(my_id==0){
		for(size_t i=0;i<block_vec.size();i++){
		clog<<"block: "<<i<<endl;
		clog<<"DC matrix: "<<block_vec[i]->A<<endl;
		}
	}*/
	//get_voltages_from_block_LU_sol();	
	solve_DC(num_procs, my_id, mpi_class);
	/*if(my_id==0)
		cout<<nodelist<<endl;
		for(size_t i=0;i<block_vec.size();i++){
			for(size_t j=0;j<block_vec[i]->count;j++)
				cout<<"b: "<<*block_vec[i]->replist[j]<<endl;
		}*/
	// return true;
	// then sync
	MPI_Barrier(MPI_COMM_WORLD);

	// return 0;
//#if 0
	for(size_t i=0;i<block_vec.size();i++){
		block_vec[i]->reset_array(block_vec[i]->bp);
		block_vec[i]->reset_array(block_vec[i]->bnewp);
	}
	
	/***** solve tran *********/
	// link transient nodes
	link_ckt_nodes(tran, my_id);
	for(size_t i=0;i<block_vec.size();i++){
		block_vec[i]->stamp_matrix_tr(my_id, mpi_class, tran);
		
	// stamp_block_matrix_tr(my_id, A, mpi_class, tran);	
		block_vec[i]->make_A_symmetric_tr(my_id, tran);	   

		block_vec[i]->stamp_current_tr(my_id, time);
		

		// if(my_id==0)
			//clog<<block_vec[i]->A<<endl;

   		block_vec[i]->CK_decomp();
		// Algebra::CK_decomp(A, block_info.L, cm);
   	/*Lp = static_cast<int *>(block_info.L->p);
   	Lx = static_cast<double*> (block_info.L->x);
   	Li = static_cast<int*>(block_info.L->i) ;
   	Lnz = static_cast<int *>(block_info.L->nz); */
   		block_vec[i]->clear_A(); //A.clear();
		// bnewp = bp
		block_vec[i]->copy_vec(block_vec[i]->bnewp,
				block_vec[i]->bp);
	}

   /*********** the following 2 parts can be implemented with pthreads ***/
   // build id_map immediately after transient factorization
   size_t n = replist.size();
#if 0
   id_map = new int [n];
   cholmod_build_id_map(CHOLMOD_A, block_info.L, cm, id_map);

   temp = new double [n];
   // then substitute all the nodes rid
   for(size_t i=0;i<n;i++){
	int id = id_map[i];
	replist[id]->rid = i;
	temp[i] = block_info.bp[i];
   }

   for(size_t i=0;i<n;i++)
	block_info.bp[i] = temp[id_map[i]];
   for(size_t i=0;i<n;i++)
        temp[i] = block_info.xp[i];
   for(size_t i=0;i<n;i++)
        block_info.xp[i] = temp[id_map[i]];
   delete [] temp;
   delete [] id_map;
#endif   
   /*for(size_t i=0;i<replist.size();i++){
	block_info.bnewp[i] = block_info.bp[i];
   }*/
  
   // set Geq for induc and capac
   set_eq_induc(tran);
   set_eq_capac(tran);

   //if(my_id==0)
	   //clog<<"before modify_rhs_tr_0. "<<endl;
   // already push back cap and induc into set_x and b
   for(size_t i=0;i<block_vec.size();i++){
   	block_vec[i]->modify_rhs_tr_0(block_vec[i]->bnewp, block_vec[i]->xp, my_id);
  	
   }
   
   //if(my_id==0)
	   //clog<<"after modify_rhs_tr_0. "<<endl;
#if 0
   // push rhs node into node_set b
   for(size_t i=0;i<n;i++){
	   if(block_info.bnewp[i] !=0)
		   pg.node_set_b.push_back(i);
   }

   // push back all nodes in output list
   vector<size_t>::iterator it;
   size_t id;
   for(size_t i=0;i<tran.nodes.size();i++){
	   if(tran.nodes[i].node == NULL) continue;
	   if(!tran.nodes[i].node->rep->is_ground()){
		   id = tran.nodes[i].node->rep->rid;
		   it = find(pg.node_set_x.begin(), pg.node_set_x.end(), 
				   id);
		   if(it == pg.node_set_x.end()){
			   pg.node_set_x.push_back(id);
		   }
	   }
   }
   
   // then push back boundary nodes into node_se_x
   push_bd_nodes(pg, my_id); 
   // get path_b, path_x, len_path_b, len_path_x
   build_path_graph();
 
   s_col_FFS = new int [len_path_b];
   s_col_FBS = new int [len_path_x];
   find_super();
#endif

   //if(my_id==0)
	   //clog<<"before first time step. "<<endl;
   //for(size_t i=0;i<replist.size();i++)
   // solve_eq_sp(block_info.xp, block_info.bnewp);
   /*if(my_id==0)
	   cout<<nodelist<<endl;
   for(size_t i=0;i<block_vec.size();i++){
	   for(size_t j=0;j<block_vec[i]->count;j++)
		   cout<<"b: "<<*block_vec[i]->replist[j]<<endl;
   }*/
   
   solve_tr_step(num_procs, my_id, mpi_class);
   // if(my_id==0)
	  //  cout<<endl<<" "<<nodelist<<endl;

   //save_tr_nodes(tran, xp);
   for(size_t i=0;i<block_vec.size();i++)
	save_ckt_nodes(tran, block_vec[i]->xp);

   time += tran.step_t;
   MPI_Barrier(MPI_COMM_WORLD);

   // return 0;
   int iter = 0;
   //if(my_id==0)
	  // clog<<"after first time step. "<<endl;
   //for(; time <= tran.tot_t; time += tran.step_t){
   while(time <= tran.tot_t && iter < 1){
	// bnewp[i] = bp[i];
	for(size_t i=0;i<block_vec.size();i++){
		block_vec[i]->copy_vec(block_vec[i]->bnewp, block_vec[i]->bp);

      		block_vec[i]->stamp_current_tr_1(time);
      		// get the new bnewp
      		block_vec[i]->modify_rhs_tr(block_vec[i]->bnewp, block_vec[i]->xp);
	}

      //if(my_id==0)
	      //clog<<" ===== step: ===== "<<my_id<<" "<<time<<endl;
      // need to add bcast function for processors
      // need to be modified into block version
      //solve_eq_sp(block_info.xp, block_info.bnewp);
      // solution stored in block_info.xp
      //if(replist.size()>0)
      solve_tr_step(num_procs, my_id, mpi_class);

      if(my_id==0)
	  cout<<endl<<" "<<nodelist<<endl;
      //save_tr_nodes(tran, xp);

      for(size_t i=0;i<block_vec.size();i++)
      	save_ckt_nodes(tran, block_vec[i]->xp);
      time += tran.step_t;
      // sync in the end of each time step
      MPI_Barrier(MPI_COMM_WORLD);
      //clog<<"after time-t step barrier. "<<my_id<<" "<<time<<endl;
      iter ++;
   }

   save_ckt_nodes_to_tr(tran);
   release_ckt_nodes(tran);
   /*delete [] s_col_FFS;
   delete [] s_col_FBS;*/
// #endif
#if 0
	/////////// release resources
	if(block_info.count > 0)
		block_info.free_block_cholmod(cm);
	//if(my_id==0) clog<<"free block info. "<<endl;
	cholmod_finish(cm);
	//clog<<"cholmod finish. "<<my_id<<endl;
#endif
        // MPI_Barrier(MPI_COMM_WORLD);
	return successful;
}
// solve blocks with mpi: multi-core
// One iteration during solving the circuit, for any block B:
// 1. update the righthand-side of the matrix of B
// 2. solve the matrix
// 3. update node voltages
// 4. track the maximum error of solution
double Circuit::solve_iteration_tr(int &my_id, int &iter,
		int&num_procs, MPI_CLASS &mpi_class){	
	float diff = .0;
	float diff_root=0;

	// 0 rank cpu will scatter all bd valuesfrom bd_x_g to bd_x
	if(iter >0){
		MPI_Scatterv(bd_x_g, bd_size_g, 
			bd_base_g, MPI_FLOAT, bd_x, bd_size, 
			MPI_FLOAT, 0, MPI_COMM_WORLD);

		assign_bd_array(my_id);
	}
	if(iter == 0){
		reset_bd_array(my_id);
		reset_replist(my_id);
		for(size_t i=0;i<block_vec.size();i++){
			block_vec[i]->copy_array(block_vec[i]->bnewp_temp, block_vec[i]->bnewp);
		}
	}
	
	for(size_t i=0;i<block_vec.size();i++){
		block_vec[i]->update_rhs(
			block_vec[i]->bnewp_temp, 
			block_vec[i]->bnewp, my_id);
	}

	if(iter==0){
		for(size_t i=0;i<block_vec.size();i++){
			block_vec[i]->reset_array(block_vec[i]->x_old);
		}
	}
	else{
		// x_old stores old solution
		for(size_t i=0;i<block_vec.size();i++){
			block_vec[i]->copy_array(block_vec[i]->x_old, block_vec[i]->xp);
		}
	}

	// new rhs store in bnewp
	for(size_t i=0;i<block_vec.size();i++){
		// if(my_id==0)
			// cout<<endl<<"before update rhs. "<<endl;	
		//block_info.solve_CK(cm);
		if(block_vec[i]->count >0){
			block_vec[i]->solve_CK_tr();
			/*if(my_id==0)
				for(size_t j=0;j<block_vec[i]->count;j++)
					cout<<"j, xp: "<<j<<" "<<block_vec[i]->xp[j]<<endl;*/
		}
		//solve_eq_sp(block_info.xp, block_info.bnewp);
			//block_info.xp = static_cast<double *>(block_info.x_ck->x);
		double local_diff = 
			block_vec[i]->modify_voltage(my_id);
		if(local_diff > diff)
			diff = local_diff;

	}
	/*if(my_id==1)
		for(int i=0;i<replist.size();i++)
			clog<<"xp: "<<i<<" "<<*replist[i]<<" "<<endl;*/

	// clog<<"diff: "<<my_id<<" "<<diff<<endl;
	assign_bd_internal_array(my_id);

	// 0 rank cpu will gather all the solution 
	// from bd_x to bd_x_g
	MPI_Gatherv(internal_x, internal_size, MPI_FLOAT, 
		internal_x_g, internal_size_g, 
		internal_base_g, MPI_FLOAT, 0, 
		MPI_COMM_WORLD);
	
	// reorder boundary array according to nbrs
	if(my_id==0){
		reorder_bd_x_g(mpi_class);
	}

	MPI_Reduce(&diff, &diff_root, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&diff_root, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);	
	return diff_root;
}
// solve blocks with mpi: multi-core
// One iteration during solving the circuit, for any block B:
// 1. update the righthand-side of the matrix of B
// 2. solve the matrix
// 3. update node voltages
// 4. track the maximum error of solution
double Circuit::solve_iteration(int &my_id, int &iter,
		int&num_procs, MPI_CLASS &mpi_class){	
	float diff = .0;
	float diff_root=0;

	// 0 rank cpu will scatter all bd valuesfrom bd_x_g to bd_x
	MPI_Scatterv(bd_x_g, bd_size_g, 
			bd_base_g, MPI_FLOAT, bd_x, bd_size, 
			MPI_FLOAT, 0, MPI_COMM_WORLD);
			
	assign_bd_array(my_id);

	/*if(my_id==1) {
		for(int i=0;i<replist.size();i++)
			clog<<i<<" "<<*replist[i]<<" "<<block_info.bp[i]<<endl;
	clog<<endl;
	}*/

	// new rhs store in bnewp and solve
	for(size_t i=0; i < block_vec.size();i++){
		

		block_vec[i]->solve_CK_DC(my_id);
		double local_diff = 
			block_vec[i]->modify_voltage(my_id);
		//if(my_id==0)
			//clog<<"local_diff: "<<local_diff<<endl;
		if(local_diff > diff)
			diff = local_diff;
	}
	// if(my_id==0)
		//clog<<"diff: "<<diff<<endl;
	assign_bd_internal_array(my_id);
	// 0 rank cpu will gather all the solution from bd_x
	// to bd_x_g
	MPI_Gatherv(internal_x, internal_size, MPI_FLOAT, 
		internal_x_g, internal_size_g, 
		internal_base_g, MPI_FLOAT, 0, 
		MPI_COMM_WORLD);
	
	// reorder boundary array according to nbrs
	if(my_id==0){
		reorder_bd_x_g(mpi_class);
	}
	
	MPI_Reduce(&diff, &diff_root, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(&diff_root, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);	
	
	//if(my_id==0) clog<<"iter, diff: "<<iter<<" "<<diff_root<<endl;
	return diff_root;
}
#if DEBUG
double Circuit::modify_voltage(int &my_id, Block &block, double * x_old){
	double max_diff = 0.0;
	//if(get_name()=="VDDA") OMEGA = 1.0;
	//else OMEGA = 1.15;
	OMEGA = 1.0;
	for(size_t i=0;i<block.count;i++){
		block.xp[i] = (1-OMEGA)*x_old[i] + OMEGA*
			block.xp[i];
		// update block nodes value
		block.nodes[i]->rep->value = block.xp[i];
		double diff = fabs(x_old[i] - block.xp[i]);
		if( diff > max_diff ) max_diff = diff;
	}

	return max_diff;
}
#endif
// given vector x that obtained from LU, set the value to the corresponding
// node in nodelist
void Circuit::get_voltages_from_LU_sol(double * x){
#if DEBUG
	for(size_t i=0;i<nodelist.size()-1;i++){
		Node * node = nodelist[i];
		size_t id = node->rep->rid;	// get rep's id in Vec
		double v = x[id];		// get its rep's value
		node->value = v;
	}
#endif
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
#if DEBUG
	for(size_t i=0;i<block_vec.size();i++)
		copy_voltages_rep();
	/*for(size_t i=0;i<nodelist.size()-1;i++){
		Node * node = nodelist[i];
		//if( node->is_mergeable() ) continue;
		size_t id = node->rep->rid;
		double v = block_info.xp[id];
		node->value = v;
		// node->rep->value = v;
	}*/
#endif
}

#if DEBUG
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
		block_info.xp[id] = replist[i]->value;
		block_info.nodes[id] = replist[i];
	}
}
#endif
void Circuit:: release_ckt_nodes(Tran &tran){
   for(size_t j=0;j<ckt_nodes.size();j++){
         ckt_nodes[j].node = NULL;
   }
}


void Circuit::make_A_symmetric(double *b, int &my_id){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p=NULL, *q=NULL, *r =NULL;

	for(it=ns.begin();it!=ns.end();it++){
           if( (*it) == NULL ) continue;
           assert( fzero((*it)->value) == false );
           if(!((*it)->ab[0]->rep->isS()==X || (*it)->ab[1]->rep->isS()==X)) continue;

           // node p points to X node
           if((*it)->ab[0]->rep->isS()==X && ((*it)->ab[0]->rep->nbr[TOP]!=NULL && 
                (*it)->ab[0]->rep->nbr[TOP]->type==INDUCTANCE)){
              p = (*it)->ab[0]->rep; q = (*it)->ab[1]->rep;
           } 
           else if((*it)->ab[1]->rep->isS()==X && ((*it)->ab[1]->rep->nbr[TOP]!=NULL && 
                (*it)->ab[1]->rep->nbr[TOP]->type==INDUCTANCE)){
              p = (*it)->ab[1]->rep; q = (*it)->ab[0]->rep;
           }

           r = p->nbr[TOP]->ab[0]->rep;
           if(r->isS()!=Y) 
              r = p->nbr[TOP]->ab[1]->rep;
	   /*if(my_id==0){
		   clog<<"modify rhs net: "<<*(*it)<<endl;
		   clog<<"p, q and r: "<<*p<<" "<<*q<<" "<<*r<<endl;
	   }*/
           size_t id = q->rid;
           double G = 1.0 / (*it)->value;
           
           b[id] += r->value * G;
        }
}

// make A symmetric for tran
void Circuit::make_A_symmetric_tr(int &my_id, Tran &tran){
#if DEBUG
	int type = INDUCTANCE;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p=NULL, *q=NULL, *r =NULL;

	for(it=ns.begin();it!=ns.end();it++){
           if( (*it) == NULL ) continue;
           assert( fzero((*it)->value) == false );
           if(!((*it)->ab[0]->rep->isS()==Y || (*it)->ab[1]->rep->isS()==Y)) continue;
           // node p points to X node
           if((*it)->ab[0]->rep->isS()==Y){
              p = (*it)->ab[0]->rep; q = (*it)->ab[1]->rep;
           } 
           else if((*it)->ab[1]->rep->isS()==Y){
              p = (*it)->ab[1]->rep; q = (*it)->ab[0]->rep;
           }           

           size_t id = q->rid;
	   double G = tran.step_t / ((*it)->value*2);
           
           //b[id] += p->value * G;
           block_info.bp[id] += block_info.xp[p->rid] *G;
	}
#endif
}

// =========== stamp block version of matrix =======
void Circuit::stamp_block_resistor(int &my_id, Net * net, Matrix &A){
#if DEBUG
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
	
	double G;	
	G = 1./net->value;
	
	if(net->flag_bd ==1)
		block_info.boundary_netlist.push_back(net);

	//if(my_id==0 && net->flag_bd ==1)
		//clog<<endl<<*net<<" "<<net->flag_bd<<endl;
	int count = 0;
	for(size_t j=0;j<2;j++){
		Node *nk = nd[j], *nl = nd[1-j];

		// if boundary net
		if(net->flag_bd ==1){
			//if(my_id==0)
				//clog<<"bd net: "<<*net<<endl;
			if(!nl->is_ground() && 
				nk->flag_bd ==0 && 
				nk->isS()!=Y &&
				(nk->nbr[TOP]==NULL || 
				 nk->nbr[TOP]->type!=INDUCTANCE)){
				// stamp value into block_ids
				size_t k1 = nk->rid;
				A.push_back(k1,k1, G);
				//if(my_id==0)
					//clog<<"nk, nl: "<<*nk<<" "<<nk->rid<<" "<<*nl<<" "<<nl->rid<<" "<<nk->is_ground()<<" "<<nl->is_ground()<<endl;

			}
		}
		// else internal net
		else if( nk->isS()!=Y && nl->isS()!=Y) {
			/*if(my_id==0){
				cout<<nk->flag_bd<<" "<<nl->flag_bd<<endl;
				//cout<<"internal net: "<<*net<<endl<<endl;
			}*/
				//clog<<"nk, nl: "<<*nk<<" "<<nk->rid<<" "<<*nl<<" "<<nl->rid<<" "<<nk->is_ground()<<" "<<nl->is_ground()<<endl;
			size_t k1 = nk->rid;
			size_t l1 = nl->rid;
			if( !nk->is_ground()&&  
          			(nk->nbr[TOP]== NULL ||
				 nk->nbr[TOP]->type != INDUCTANCE)){
				/*if(my_id==0)// && nk->isS()!=Z && nk->isS()!=X && nl->isS()!=Z && nl->isS()!=Z)
				{	//cout<<*net<<endl;
					//cout<<*nk<<" "<<k1<<endl;
					cout<<"push + ("<<k1<<","<<k1<<","<<G<<")"<<endl;
				}*/
				A.push_back(k1,k1, G);
				// count ++;
			}
			/*else{
				if(my_id==0)
					cout<<"except net: "<<*net<<" "<<*nk<<endl;
			}*/
			if(!nk->is_ground() && 
				!nl->is_ground() && 
				l1 < k1 &&
				(nl->nbr[TOP] ==NULL ||
				 nl->nbr[TOP]->type != INDUCTANCE)){ // only store the lower triangular part{
				//if(my_id==0) cout<<"push - ("<<k1<<","<<l1<<","<<-G<<")"<<endl;
				//if(nk->isS()!= X && nl->isS()!=X && nk->isS()!=Z && nl->isS()!=Z){
				
				//if(my_id==0) cout<<"push -("<<*nk<<" "<<*nl<<" ("<<k1<<","<<l1<<","<<-G<<")"<<endl;
				A.push_back(k1,l1,-G);
				//}
			}
		}
	}// end of for j	
	/*if(net->flag_bd ==0 && count <2 && my_id==0){
		clog<<nd[0]->rid<<" "<<nd[1]->rid<<endl;
		clog<<endl<<count<<endl;
		clog<<*net<<" stamp error. "<<endl;
	}*/
#endif
}

// only stamp resistor node connected to inductance
void Circuit::stamp_block_resistor_tr(int &my_id, Net * net, Matrix &A){
	// skip boundary nets
	if(net->flag_bd ==1)
		return;

	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
	double G;	
	G = 1./net->value;

	for(size_t j=0;j<2;j++){
		Node *nk = nd[j], *nl = nd[1-j];

		// skip resistor layer net
		if(nk->isS()!=X && nl->isS()!=X) continue;
		if(nk->isS()!=X) continue;

		size_t k1 = nk->rid;
		size_t l1 = nl->rid;
		if( !nk->is_ground()&& 
          		(nk->nbr[TOP]!= NULL &&
			 nk->nbr[TOP]->type == INDUCTANCE))
			A.push_back(k1,k1, G);
		// also need to push back connection
		if(!nl->is_ground() 
		    &&(nl->nbr[TOP] ==NULL ||
		    nl->nbr[TOP]->type != INDUCTANCE))
			if(l1 < k1){
				//if(my_id==0) clog<<"push ("<<k1<<","<<l1<<","<<-G<<")"<<endl;
				A.push_back(k1,l1,-G);
			}
			else if(l1>k1){
				//if(my_id==0) clog<<"push ("<<l1<<","<<k1<<","<<-G<<")"<<endl;
				A.push_back(l1, k1, -G);
			}
	}// end of for j	
}

// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void Circuit::modify_rhs_l_tr_0(Net *net, double *rhs, double *x, int &my_id){
	//clog<<"l net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	if(nk->isS() !=X){ 
		swap<Node*>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
	size_t k = nk->rid;
	size_t l = nl->rid;
	double Ieq = 0;

	double i_t = 0;
	double temp = 0;
	//temp = tran.step_t / (2*net->value) * 
		//(nl->value - nk->value);
	//temp = tran.step_t / (2*net->value)*(x[l] - x[k]);
	temp = net->value *(x[l] - x[k]);
	//if(nk->value != x[k] || nl->value != x[l])
	   //clog<<"k, l, x_k, x_l: "<<nk->value<<" "<<nl->value<<" "<<
	     //x[k]<<" "<<x[l]<<endl;

	//clog<<"delta_t/2L, nl-nk, temp: "<<tran.step_t / (2*net->value)<<" "<<(nl->value-nk->value)<<" "<<temp<<endl;
	
	Net *r = nk->nbr[BOTTOM];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to X node
	if(a->isS()!=X) {
		swap<Node*>(a, b);
		swap<Node*>(r->ab[0], r->ab[1]);
	}
	size_t id_a = a->rid;
	size_t id_b = b->rid;
	i_t = (x[id_a] - x[id_b]) / r->value;
	//i_t = (a->value - b->value) / r->value;
        //if(b->value != x[id_b] || a->value != x[id_a])
	   //clog<<"a, b, x_a, x_b: "<<a->value<<" "<<b->value<<" "<<
	     //x[id_a]<<" "<<x[id_b]<<endl;

	//clog<<"resiste r: "<<*r<<endl;
	//clog<<*a<<" "<<*b<<endl;
	//clog<<"a, b, r, i_t: "<<a->value<<" "<<b->value<<" "<<
		//r->value<<" "<<i_t<<endl;
       
        // push inductance nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
        //clog<<*b<<" "<<id_b<<endl;
//#if 0
        pg.node_set_x.push_back(k);
        pg.node_set_x.push_back(id_b);
//#endif
	Ieq  = i_t + temp;
	if(nk->isS() !=Y && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD circuit
	}
	if(nl->isS()!=Y && !nl->is_ground()){
		 rhs[l] += -Ieq; // VDD circuit
	}
}



// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void Circuit::modify_rhs_l_tr(Net *net, double *rhs, double *x){
	//clog<<"l net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	//if(nk->isS() !=X){ 
		//swap<Node*>(nk, nl);
	//}
	size_t k = nk->rid;
	size_t l = nl->rid;
	double Ieq = 0;

	double i_t = 0;
	double temp = 0;
	//temp = tran.step_t / (2*net->value) * 
		//(nl->value - nk->value);
	//temp = tran.step_t / (2*net->value)*(x[l] - x[k]);
	temp = net->value *(x[l] - x[k]);	
	//if(nk->value != x[k] || nl->value != x[l])
	   //clog<<"k, l, x_k, x_l: "<<nk->value<<" "<<nl->value<<" "<<
	     //x[k]<<" "<<x[l]<<endl;

	//clog<<"delta_t/2L, nl-nk, temp: "<<tran.step_t / (2*net->value)<<" "<<(nl->value-nk->value)<<" "<<temp<<endl;
	
	Net *r = nk->nbr[BOTTOM];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to X node
	//if(a->isS()!=X) {
		//swap<Node*>(a, b);
	//}
	size_t id_a = a->rid;
	size_t id_b = b->rid;
	i_t = (x[id_a] - x[id_b]) / r->value;
	//i_t = (a->value - b->value) / r->value;
        //if(b->value != x[id_b] || a->value != x[id_a])
	   //clog<<"a, b, x_a, x_b: "<<a->value<<" "<<b->value<<" "<<
	     //x[id_a]<<" "<<x[id_b]<<endl;

	//clog<<"resiste r: "<<*r<<endl;
	//clog<<*a<<" "<<*b<<endl;
	//clog<<"a, b, r, i_t: "<<a->value<<" "<<b->value<<" "<<
		//r->value<<" "<<i_t<<endl;
       
        // push inductance nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
        //clog<<*b<<" "<<id_b<<endl;
 #if 0
        if(iter==0){
           pg.node_set_x.push_back(k);
           pg.node_set_x.push_back(id_b);
        }
#endif
	Ieq  = i_t + temp;
	//clog<<"Ieq: "<<Ieq<<endl;
	if(nk->isS() !=Y && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD circuit
		//clog<<*nk<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=Y && !nl->is_ground()){
		 rhs[l] += -Ieq; // VDD circuit
		//clog<<*nl<<" "<<rhs[l]<<endl;
	}
}

void Circuit::stamp_block_current(int &my_id, Net * net, MPI_CLASS &mpi_class){
#if DEBUG
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;

	// only stamp for internal node
	if( !nk->is_ground() && nk->isS()!=Y && nk->flag_bd ==0) {
		size_t k = nk->rid;

		block_info.bp[k] += -net->value;
		//pk[k] += -net->value;
	}
	if( !nl->is_ground() && nl->isS()!=Y && nl->flag_bd ==0) {
		size_t l = nl->rid;

		block_info.bp[l] += net->value;
		//if(my_id==1) clog<<"bk: "<<l<<" "<<block_info.bp[l]<<endl;
		//pl[l] +=  net->value;
	}
#endif
}

void Circuit::stamp_block_VDD(int &my_id, Net * net, Matrix &A){
#if DEBUG
	// find the non-ground node
	Node * X = net->ab[0];
	//if(my_id==0) clog<<"net: "<<*net<<endl;
	if( X->is_ground() ) X = net->ab[1];

	if(X->rep->flag_bd ==1) return;
	// do stamping for internal node
	long id =X->rep->rid;

	//if(my_id==0) clog<<" stamp net: "<<*net<<endl;
	A.push_back(id, id, 1.0);
	//if(my_id==0) clog<<"push ("<<id<<", "<<id<<", 1.0)"<<endl;
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
			south->type == CURRENT ){
		// this node connects to a VDD and a current
		// the current should be stamped
		//assert( feqn(1.0, q[id]) ); 
		assert( feqn(1.0, block_info.bp[id]) );
		block_info.bp[id] = net->value;
		//q[id] = net->value;	    // modify it
	}
	else{
		block_info.bp[id] += net->value;
		//q[id] += net->value;
	}
#endif
}

// stamp a voltage source
void Circuit::stamp_block_VDD_tr(int &my_id, Net * net, Matrix &A){
#if DEBUG
	// find the non-ground node
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];

	if(X->rep->flag_bd ==1) return;
	// do stamping for internal node
	long id =X->rep->rid;
	// A.push_back(id, id, 1.0);
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
	    south->type == CURRENT ){
		// this node connects to a VDD and a current
		// the current should be stamped
		//assert( feqn(1.0, q[id]) ); 
		assert( feqn(1.0, block_info.bp[id]) );
		block_info.bp[id] = net->value;
		//q[id] = net->value;	    // modify it
	}
	else{
		block_info.bp[id] += net->value;
		//q[id] += net->value;
	}
#endif
}

// all cores stamp dc inductance
void Circuit::stamp_inductance_dc(Matrix & A, Net * net, int &my_id){
#if DEBUG
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	//if(my_id==0)
		//clog<<"induc net: "<<*net<<endl;
	if( nk->isS()!=Y && !nk->is_ground()&&nk->flag_bd==0) {
		A.push_back(k,k, 1);
		//if(my_id==0)
		//clog<<"push ("<<*nk<<" "<<k<<" "<<k<<" 1)"<<endl;
		// general stamping
		if(!nl->is_ground())
		// A.push_back(k,l,-1);
		// make it symmetric
			block_info.bp[k] = block_info.bp[l];
		//clog<<"("<<k<<" "<<k<<" "<<1<<")"<<endl;
		//clog<<"("<<k<<" "<<l<<" "<<-1<<")"<<endl;
	}

	if( nl->isS() !=Y && !nl->is_ground() &&nl->flag_bd==0) {
		A.push_back(l,l, 1);
		//if(my_id==0)
		//clog<<"push ("<<*nl<<" "<<l<<" "<<l<<" 1)"<<endl;
		if(!nk->is_ground())
		// general stamping
		// A.push_back(l,k,-1);
		block_info.bp[l] = block_info.bp[k];

		//clog<<"("<<l<<" "<<l<<" "<<1<<")"<<endl;
		//clog<<"("<<l<<" "<<k<<" "<<-1<<")"<<endl;
	}
#endif
}

// stamp capacitance Geq = 2C/delta_t
void Circuit::stamp_capacitance_tr(Matrix &A, Net *net, Tran &tran, int &my_id){
	//clog<<"net: "<<*net<<endl;
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	// Geq = 2*C / delta_t
	Geq = (2*net->value) / tran.step_t;
	//net->value = Geq;
	//clog<<"C delta_t Geq: "<<net->value<<" "<<tran.step_t<<" "<<Geq<<endl;
	// Ieq = i(t) + 2*C / delta_t * v(t)

	if( nk->isS()!=Y  && !nk->is_ground() && nk->flag_bd==0) {
		A.push_back(k,k, Geq);
		//clog<<"("<<k<<" "<<k<<" "<<Geq<<")"<<endl;
		if(!nl->is_ground()&& k > l){
			A.push_back(k,l,-Geq);
			//clog<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=Y && !nl->is_ground() && nl->flag_bd==0) {
		A.push_back(l,l, Geq);
		//clog<<"("<<l<<" "<<l<<" "<<Geq<<")"<<endl;
		if(!nk->is_ground()&& l > k){
			A.push_back(l,k,-Geq);
			//clog<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
	}
}

// stamp inductance Geq = delta_t/(2L)
void Circuit::stamp_inductance_tr(Matrix & A, Net * net, Tran &tran, int &my_id){
	//clog<<"net: "<<*net<<endl;
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	// Geq = delta_t / (2*L)
	Geq = tran.step_t / (2*net->value);
	//net->value = Geq;

	if( nk->isS()!=Y  && !nk->is_ground() && nk->flag_bd==0) {
		// -1 is to clear formal inserted 1 at (k,k)
		A.push_back(k,k, Geq-1);
		//clog<<"("<<k<<" "<<k<<" "<<Geq-1<<")"<<endl;
		//clog<<nl->isS()<<endl;
		if(!nl->is_ground()&& nl->isS()!=Y && k>l){
			A.push_back(k,l,-Geq);
		        //clog<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=Y && !nl->is_ground() && nl->flag_bd==0) {
		// -1 is to clear formal inserted 1 at (l,l)
		A.push_back(l,l, Geq-1);
		//clog<<"("<<l<<" "<<l<<" "<<Geq-1<<")"<<endl;
		if(!nk->is_ground() && nk->isS()!=Y && l>k){
			A.push_back(l,k,-Geq);
			//clog<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
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
	//clog<<"newly added net: "<<*net<<endl;
	this->add_net(net);
}

void Circuit::boundary_init(int &my_id, int &num_procs){
	assign_bd_base(my_id);

	assign_bd_dd_size(my_id);

	// allocate boundary_nodelist size
	bd_dd_size_g = new int [8*num_procs];	
	
	//if(my_id==0)
		//clog<<"before gathering bd_dd_size. "<<my_id<<endl;
	MPI_Gather(bd_dd_size, 8, MPI_INT, bd_dd_size_g, 
			8, MPI_INT, 0, MPI_COMM_WORLD);

	
	bd_x = new float[bd_size];

	// assign bd_size_g value
	bd_size_g = new int[num_procs];
	MPI_Gather(&bd_size, 1, MPI_INT, bd_size_g, 1, MPI_INT,
			0, MPI_COMM_WORLD);

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
		internal_nodelist_sw, internal_x, my_id);
	//clog<<"s: "<<endl;
	assign_bd_internal_array_dir(internal_base[1], 
		internal_nodelist_s, internal_x, my_id);
	//clog<<"se: "<<endl;
	assign_bd_internal_array_dir(internal_base[2], 
		internal_nodelist_se, internal_x, my_id);
	//clog<<"w: "<<endl;
	assign_bd_internal_array_dir(internal_base[3], 
		internal_nodelist_w, internal_x, my_id);
	//clog<<"e: "<<endl;
	assign_bd_internal_array_dir(internal_base[4], 
		internal_nodelist_e, internal_x, my_id);
	//clog<<"nw: "<<endl;
	assign_bd_internal_array_dir(internal_base[5], 
		internal_nodelist_nw, internal_x, my_id);
	//clog<<"n: "<<endl;
	assign_bd_internal_array_dir(internal_base[6], 
		internal_nodelist_n, internal_x, my_id);
	//clog<<"ne: "<<endl;
	assign_bd_internal_array_dir(internal_base[7], 
		internal_nodelist_ne, internal_x, my_id);
}

// assign 4 boundary internal nodes value, store
// them in array bd_x
void Circuit::assign_bd_internal_array_dir(int &base, NodePtrVector & list, float *internal_x, int &my_id){
	Node *nd;
	for(size_t i=0;i<list.size();i++){
		nd = list[i]->rep;
		internal_x[base+i] = nd->value;
		//if(my_id==0)
		//clog<<base+i<<" "<<*nd<<endl;
	}
}

// assign 4 boundary nodes value from receiving
// array bd_x
void Circuit::assign_bd_array(int &my_id){
	int base =0;
	base = bd_base[0];
	assign_bd_array_dir(base, bd_nodelist_sw, my_id);
	
	base = bd_base[1];
	assign_bd_array_dir(base, bd_nodelist_s, my_id);	
	
	base = bd_base[2];
	assign_bd_array_dir(base, bd_nodelist_se, my_id);	

	base = bd_base[3];
	assign_bd_array_dir(base, bd_nodelist_w, my_id);

	base = bd_base[4];
	assign_bd_array_dir(base, bd_nodelist_e, my_id);

	base = bd_base[5];
	assign_bd_array_dir(base, bd_nodelist_nw, my_id);

	base = bd_base[6];
	assign_bd_array_dir(base, bd_nodelist_n, my_id);

	base = bd_base[7];
	assign_bd_array_dir(base, bd_nodelist_ne, my_id);
}

// reset 4 boundary nodes value to 0
void Circuit::reset_bd_array(int &my_id){
	size_t n = bd_nodelist_sw.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_sw[i]->rep->value = 0;
	
	n = bd_nodelist_s.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_s[i]->rep->value = 0;

	n = bd_nodelist_se.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_se[i]->rep->value = 0;

	n = bd_nodelist_w.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_w[i]->rep->value = 0;

	n = bd_nodelist_e.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_e[i]->rep->value = 0;

	n = bd_nodelist_nw.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_nw[i]->rep->value = 0;

	n = bd_nodelist_n.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_n[i]->rep->value = 0;

	n = bd_nodelist_ne.size();
	for(size_t i=0;i<n;i++)
		bd_nodelist_ne[i]->rep->value = 0;
}

void Circuit::reset_replist(int &my_id){
	for(size_t i=0;i<replist.size();i++){
		if(replist[i]->isS()!=Y)
			replist[i]->value = 0;
	}
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
			//clog<<"notrhwest bid: "<<bid_nbr<<endl;
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
			//clog<<"north east bid: "<<bid_nbr<<endl;
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

void Circuit::internal_init(int &my_id, int &num_procs){
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
		internal_dd_size_g, 8, MPI_INT, 0, 
		MPI_COMM_WORLD);

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
			0, MPI_COMM_WORLD);
	
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
void Circuit::assign_bd_array_dir(int &base, NodePtrVector &list, int &my_id){
	Node *p;
	for(size_t i=0;i<list.size();i++){
		p = list[i]->rep;
		p->value = bd_x[base+i];
		// if(my_id==0)
			// clog<<"bd_array: "<<*p<<endl;
	}
}

// assign value back to transient nodes
void Circuit:: save_ckt_nodes(Tran &tran, double *x){
   size_t id=0;
   for(size_t j=0;j<ckt_nodes.size();j++){
	 //cout<<"nodes: "<<ckt_nodes[j].node->name<<endl;
         id = ckt_nodes[j].node->rep->rid;
	 //cout<<"value: "<<x[id]<<endl;
         ckt_nodes[j].value.push_back(x[id]);
      }
}

void Circuit::save_ckt_nodes_to_tr(Tran &tran){
	int index = 0;
	double time = 0;
	int j=0;
	double value = 0;
	for(size_t i=0;i<ckt_nodes.size();i++){
		index = ckt_nodes[i].flag;
		//cout<<"ckt_nodes, index: "<<ckt_nodes[i].node->name
			//<<" "<<ckt_nodes[i].flag<<endl;
		//tran.nodes[index].node = ckt_nodes[i];
		j=0;
		for(time = 0; time < tran.tot_t; time+=tran.step_t){
			value = ckt_nodes[i].value[j];
			//cout<<"time, value: "<<time<<" "<<value<<endl;
			tran.nodes[index].value.push_back(value);
			j++;
		}
	}	
}

// link transient nodes with nodelist
void Circuit:: link_ckt_nodes(Tran &tran, int &my_id){ 
   Node_TR_PRINT nodes_temp;
   for(size_t i=0;i<nodelist.size();i++){
      for(size_t j=0;j<tran.nodes.size();j++){
         if(nodelist[i]->name == tran.nodes[j].name){
	    // record the index in tran.nodes
	    nodes_temp.flag = j;
	    //if(my_id==3)
	    	//clog<<"tran.nodes, index: "<<
		//nodelist[i]->name<<" "<<nodes_temp.flag<<endl;
	    nodes_temp.node = nodelist[i];
	    ckt_nodes.push_back(nodes_temp);
            // record the id for ckt_nodes
            break;
         }
      }
   }
}

// stamp transient current values into rhs
void Circuit::stamp_current_tr(int &my_id, double &time){
	NetPtrVector & ns = net_set[CURRENT];
	for(size_t i=0;i<ns.size();i++)
		stamp_current_tr_net(ns[i], time, my_id);
}

// stamp transient current values into rhs
void Circuit::stamp_current_tr_1(double *bp, double *b, double &time){
	NetPtrVector & ns = net_set[CURRENT];
	for(size_t i=0;i<ns.size();i++)
		stamp_current_tr_net_1(bp, b, ns[i], time);
}

void Circuit::stamp_current_tr_net_1(double *bp, double * b, Net * net, double &time){
	double diff = 0;
	double current = net->value;
	current_tr(net, time);
	//if(time / 1e-11 == 22)
	//cout<<"time, co, cn: "<<time<<" "<<current<<" "<<net->value<<endl;
	// only stamps when net got a different current
	if(current != net->value){
		diff = net->value - current;
		
		//cout<<"time, old, new, diff: "<<time <<" "<<current<<" "<<net->value<<" "<<diff<<endl;
		//cout<<"net: "<<*net;
		//clog<<"current: "<<current<<endl;
		Node * nk = net->ab[0]->rep;
		Node * nl = net->ab[1]->rep;
		if( !nk->is_ground()&& nk->isS()!=Y) { 
			size_t k = nk->rid;
			//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
			//clog<<"time, k, b bef: "<<time<<" "<<k<<" "<<b[k]<<endl;
			b[k] += -diff;//current;
			bp[k] = b[k];
			//clog<<"time, k, b: "<<time <<" "<<k<<" "<<b[k]<<endl;
		}
		if( !nl->is_ground() && nl->isS()!=Y) {
			size_t l = nl->rid;
			//clog<<"time, l, b bef: "<<time<<" "<<l<<" "<<b[l]<<endl;
			//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
			b[l] +=  diff;// current;
			bp[l] = b[l];
			//clog<<"time, l, b: "<<time<<" "<<l<<" "<<b[l]<<endl;
		}
	}
}

void Circuit::stamp_current_tr_net(Net * net, double &time, int &my_id){
#if DEBUG
	current_tr(net, time);
	//clog<<"net: "<<*net<<endl;
	//clog<<"current: "<<current<<endl;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground()&& nk->isS()!=Y) { 
		size_t k = nk->rid;
		//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
		block_info.bp[k] += -net->value;//current;
		//clog<<"time, k, b: "<<time<<" "<<k<<" "<<b[k]<<endl;
	}
	if( !nl->is_ground() && nl->isS()!=Y) {
		size_t l = nl->rid;
		//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
		block_info.bp[l] +=  net->value;// current;
		//clog<<"time, l, b: "<<time<<" "<<l<<" "<<b[l]<<endl;
	}
#endif
}

// decide transient step current values
void Circuit::current_tr(Net *net, double &time){
	double slope = 0;
	double Tr = net->tr[3]; // Tr
	double PW = Tr + net->tr[5]; // PW
	double Tf = PW + net->tr[4]; // Tf
	double t_temp = time - net->tr[2]; // TD
	double t = fmod(t_temp, net->tr[6]); // Period
	if(time <= net->tr[2])// TD
		net->value = net->tr[0];// V1
	else if(t > 0 && t<= Tr){
		slope = (net->tr[1] - net->tr[0]) / 
			(net->tr[3]);
		net->value = net->tr[0] + t*slope;
	}
	else if(t > Tr && t<= PW)
		net->value = net->tr[1];
	else if(t>PW && t<=Tf){
		slope = (net->tr[0]-net->tr[1])/(net->tr[4]);
		net->value = net->tr[1] + slope*(t-PW);
	}
	else
		net->value = net->tr[0];
	//return current;
}

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void Circuit::modify_rhs_c_tr_0(Net *net, double * rhs, double *x, int &my_id){
	double i_t = 0;
	double temp = 0;
	double Ieq = 0;
	//if(my_id==0)
	//clog<<"c net: "<<*net<<" net->flag_bd: "<<net->flag_bd<< endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
        // nk point to Z node
        if(nk->isS() != Z){
		swap<Node *>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
	//if(my_id==0)
	//clog<<"nk, nl: "<<*nk<<" "<<*nl<<endl;
	size_t k = nk->rid;
	size_t l = nl->rid;
	//if(my_id==0)
	//clog<<"k, l: "<<k<<" "<<l<<" "<<nk->flag_bd<<" "<<nl->flag_bd<<endl;

	Net *r = nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node
	if(a->isS()!=Z) {
		swap<Node *>(a, b);
		swap<Node*>(r->ab[0], r->ab[1]);
	}
	//clog<<"a, b: "<<*a<<" "<<*b<<endl;

	size_t id_a = a->rid;
	size_t id_b = b->rid;
	//i_t = (b->value - a->value) / r->value;
	i_t = (x[id_b] - x[id_a]) / r->value;
	//if(b->value != x[id_b] || a->value != x[id_a])
	   //cout<<"a, b, x_a, x_b: "<<a->value<<" "<<b->value<<" "<<
	     //x[id_a]<<" "<<x[id_b]<<endl;
	//clog<<"i_t: "<<i_t<<endl;
	//temp = 2*net->value / tran.step_t * 
		//(nk->value - nl->value);
       
        // push 2 nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
//#if 0
        pg.node_set_x.push_back(k);
        if(!nl->is_ground()) {
              //clog<<*nl<<" "<<l<<endl;
           pg.node_set_x.push_back(l);
        }
        else if(!b->is_ground()){
              //clog<<*b<<" "<<id_b<<endl;
           pg.node_set_x.push_back(id_b);
        }
//#endif
	//if(my_id==0)
	//	clog<<"before calc temp. "<<endl;
	if(nk->is_ground())
	 //temp = 2*net->value/tran.step_t*(0-x[l]);
	 temp = net->value *(-x[l]);
        else if(nl->is_ground()){
         //temp = 2*net->value/tran.step_t *(x[k]);
	 temp = net->value *x[k];
        }
        else
         //temp = 2*net->value/tran.step_t *(x[k] - x[l]);
	 temp = net->value *(x[k]-x[l]);
	//if(nk->value != x[k] || nl->value != x[l])
	   //cout<<"k, l, x_k, x_l: "<<nk->value<<" "<<nl->value<<" "<<
	     //x[k]<<" "<<x[l]<<endl;
	//clog<<"nk-nl "<<(nk->value - nl->value)<<" "<<2*net->value/tran.step_t<<" "<<temp<<endl;
	
	Ieq  = (i_t + temp);
	//if(my_id==0)
	//clog<< "Ieq is: "<<Ieq<<endl;
	//clog<<"Geq is: "<<2*net->value / tran.step_t<<endl;
	if(!nk->is_ground()&& nk->isS()!=Y){
		 rhs[k] += Ieq;	// for VDD circuit
		 //if(my_id==1)
		    // clog<<k<<" "<<*nk<<" rhs +: "<<rhs[k]<<endl;
	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		 rhs[l] += -Ieq; 
		 //if(my_id==1)
		    // clog<<l<<" "<<*nl<<" rhs +: "<<rhs[l]<<endl;
	}
	//if(my_id==0)
	//	clog<<"finish 1 net. "<<endl;
}

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void Circuit::modify_rhs_c_tr(Net *net, double * rhs, double *x){
	double temp = 0;
	//clog<<"c net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
        
	// nk point to Z node
	size_t k = nk->rid;
	size_t l = nl->rid;

	Net *r = nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node

	size_t id_a = a->rid;
	size_t id_b = b->rid;
	double i_t = (x[id_b] - x[id_a]) / r->value;
	
	if(nk->is_ground())
	 temp = net->value *(-x[l]);
        else if(nl->is_ground())
	 temp = net->value *x[k];
        else
	 temp = net->value *(x[k]-x[l]);
	
	double Ieq  = i_t + temp;
	if(!nk->is_ground()&& nk->isS()!=Y){
		 rhs[k] += Ieq;	// for VDD circuit
	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		 rhs[l] -= Ieq; 
	}
}

void Circuit::set_eq_induc(Tran &tran){
	NetPtrVector &ns = net_set[INDUCTANCE];
	for(size_t i=0;i<ns.size();i++)
		ns[i]->value = tran.step_t /(2*ns[i]->value);
}

void Circuit::set_eq_capac(Tran &tran){
	NetPtrVector &ns = net_set[CAPACITANCE];
	for(size_t i=0;i<ns.size();i++)
		ns[i]->value = 2*ns[i]->value/tran.step_t;
}

// update rhs by transient nets
void Circuit::modify_rhs_tr_0(double * b, double *x, int &my_id){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){	
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_c_tr_0(ns[i], b, x, my_id);
			}
		}
		else if(type == INDUCTANCE){
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_l_tr_0(ns[i], b, x, my_id);	
			}
		}
	}
}

// update rhs by transient nets
void Circuit::modify_rhs_tr(double * b, double *x){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){
			for(size_t i=0;i<ns.size();i++)
				modify_rhs_c_tr(ns[i], b, x);
		}
		else if(type == INDUCTANCE){
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_l_tr(ns[i], b, x);
			}
		}
	}
}

void Circuit::parse_path_table(){
   // build up nodelist info
      Node_G *node;
      for(size_t i=0;i<replist.size();i++){
         node = new Node_G();
         node->value = i;
         pg.nodelist.push_back(node);
      }
}

void Circuit::build_path_graph(){
   build_FFS_path();

   build_FBS_path();

   // only keep the 2 paths, switch from List into array
   len_path_b = pg.path_FFS.get_size();
   len_path_x = pg.path_FBS.get_size();

   path_b = new int[len_path_b];
   path_x = new int [len_path_x];
   
   Node_G *nd;
   nd = pg.path_FFS.first;
   for(int i=0;i<len_path_b;i++){
      path_b[i] = nd->value;
      if(nd->next != NULL)
         nd = nd->next;
   }
   pg.path_FFS.destroy_list();

   nd = pg.path_FBS.first;
   for(int i=0;i<len_path_x;i++){
      path_x[i] = nd->value;
      if(nd->next != NULL)
         nd = nd->next;
   }
   pg.path_FBS.destroy_list();

   pg.nodelist.clear();
   pg.node_set_b.clear();
   pg.node_set_x.clear();
}

void Circuit::build_FFS_path(){
   parse_path_table(); 

   set_up_path_table();

   find_path(pg.node_set_b, pg.path_FFS);

   pg.path_FFS.assign_size();

   for(size_t i=0;i<replist.size();i++)
      pg.nodelist[i]->flag = 0;
}

void Circuit::build_FBS_path(){
  pg.nodelist.clear();
  parse_path_table();
  set_up_path_table();
  find_path(pg.node_set_x, pg.path_FBS);
   pg.path_FBS.assign_size();
}

void Circuit::set_up_path_table(){
#if DEBUG
	size_t n = block_info.L->n;
   //int *Lp, *Li, *Lnz;
   int p, lnz, s, e;
   //Lp = static_cast<int *> (L->p);
   //Lnz = (int *)L->nz;
   //Li = static_cast<int *> (L->i);
   for(size_t i=0;i<n;i++){
      p = Lp[i];
      lnz = Lnz[i];

      s = Li[p];
      e = s;
      if(lnz >1) 
         e = Li[p+1];

      if(s<e)
         pg.nodelist[s]->next = pg.nodelist[e];
   }
#endif
}

bool compare_Node_G(const Node_G *nd_1, const Node_G *nd_2){
   return (nd_1->value < nd_2->value);
 }

void Circuit::find_path(vector<size_t> &node_set, List_G &path){
   Node_G* ne = pg.nodelist[pg.nodelist.size()-1];
   vector <Node_G *> insert_list;
   sort(node_set.begin(), node_set.end());
   if(node_set.size() == 0) return;
   
   // build up the first path start with id = min 
   int id = node_set[0];
   do{
      path.add_node(pg.nodelist[id]);
      pg.nodelist[id]->flag =1;
      if(pg.nodelist[id]->next == NULL) break;
      pg.nodelist[id] = pg.nodelist[id]->next;
   }while(pg.nodelist[id]->value != ne->value);
   path.add_node(ne);

   for(size_t i=0; i<node_set.size();i++){
      int id = node_set[i];
      if(pg.nodelist[id]->flag == 1) continue;
      // stops at first place where flag = 0
      // clog<<"stops at i: "<<i<<endl;
      do{
         if(pg.nodelist[id]->flag ==0){
            insert_list.push_back(pg.nodelist[id]);
            pg.nodelist[id]->flag =1;
         }
         if(pg.nodelist[id]->next == NULL || 
           pg.nodelist[id]->next->flag ==1)
            break;
         // else clog<<"next node: "<<*pg.nodelist[id]->next;
         pg.nodelist[id] = pg.nodelist[id]->next;
      }while(pg.nodelist[id]->value != ne->value); 
   }

   //clog<<"insert_list.size: "<<insert_list.size()<<endl;
   sort(insert_list.begin(), insert_list.end(), compare_Node_G);
   //for(int i=0;i<insert_list.size();i++)
      //clog<<"i, insert: "<<i<<" "<<*insert_list[i]<<endl;

   //clog<<"path: "<<&path<<endl;
   // p is the old pointer to the list
   // will be updated into new one
   Node_G *q=NULL;
   Node_G *p=NULL;
//#if 0
   for(size_t k=0;k<insert_list.size();k++){
      if(k ==0) p = path.first;
      else p = q;
      q = path.insert_node(insert_list[k], p);
   }
//#endif

   //clog<<"path: "<<&path<<endl;
   //clog<<endl;
   insert_list.clear();
}

 // find super node columns for path_b and path_x
void Circuit::find_super(){
#if DEBUG
	int p, lnz;
    int j, k;
    // FFS loop
    for(k=0;k<len_path_b;){
       j = path_b[k];
       p = Lp[j];
       lnz = Lnz[j];
       if (lnz < 4 || path_b[k+1]!=j+1 || lnz != Lnz [j+1] + 1
         || Li [p+1] != j+1){
          s_col_FFS[k]= 1;
          k++;
       }
       else if (path_b[k+2]!= j+2 || lnz != Lnz [j+2] + 2
         || Li [p+2] != j+2){
         s_col_FFS[k]=2;
         k+=2;
       }
       else{
          s_col_FFS[k]=3;
          k+=3;
       }
    }
    //FBS loop
    for(k=len_path_x-1;k>=0;){
       j = path_x[k];
       p = Lp[j];
       lnz = Lnz[j];
       if (j < 4 || path_x[k-1]!=j-1||lnz != Lnz [j-1] - 1
         || Li [Lp [j-1]+1] != j){
         s_col_FBS[k]=1;
         k--;
       }
       else if (path_x[k-2] != j-2 ||lnz != Lnz [j-2]-2 ||
         Li [Lp [j-2]+2] != j){
         s_col_FBS[k]=2;
         k-=2;
       }
       else{
          s_col_FBS[k]=3;
          k-=3;
       }
    }
#endif
 }
 
 void Circuit::solve_eq_sp(double *X, double *bnewp){
#if DEBUG
	 int p, q, r, lnz, pend;
    int j, k, n = block_info.L->n ;
    for(int i=0;i<n;i++){
       X[i] = bnewp[i];
    }
    // FFS solve
    for(k=0; k < len_path_b;){
       j = path_b[k];
       //for (j = 0 ; j < n ; ){
       /* get the start, end, and length of column j */
       p = Lp [j] ;
       lnz = Lnz [j] ;
       pend = p + lnz ;
 
       if (s_col_FFS[k]==1)//lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
       {
 
          /* -------------------------------------------------------------- */
          /* solve with a single column of L */
          /* -------------------------------------------------------------- */
          double y = X [j] ;
          if(block_info.L->is_ll == true){
             X[j] /= Lx [p] ;
          }
          for (p++ ; p < pend ; p++)
          {
             X [Li [p]] -= Lx [p] * y ;
          }
          k++ ;  /* advance to next column of L */
 
       }
       else if (s_col_FFS[k]==2)//lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
       {
          {
             double y [2] ;
             q = Lp [j+1] ;
             if(block_info.L->is_ll == true){
                y [0] = X [j] / Lx [p] ;
                y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
                X [j  ] = y [0] ;
                X [j+1] = y [1] ;
             }
 
             else{
                y [0] = X [j] ;
                y [1] = X [j+1] - Lx [p+1] * y [0] ;
                X [j+1] = y [1] ;
             }
             for (p += 2, q++ ; p < pend ; p++, q++)
             {
                X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] ;
             }
 
          }
          k += 2 ;           /* advance to next column of L */
 
       }
       else
       {
          {
             double y [3] ;
             q = Lp [j+1] ;
             r = Lp [j+2] ;
             //#ifdef LL
             if(block_info.L->is_ll == true){
                y [0] = X [j] / Lx [p] ;
                y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
                y [2] = (X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1]) / Lx [r] ;
                X [j  ] = y [0] ;
                X [j+1] = y [1] ;
                X [j+2] = y [2] ;
             }
 
             else{
                y [0] = X [j] ;
                y [1] = X [j+1] - Lx [p+1] * y [0] ;
                y [2] = X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1] ;
                X [j+1] = y [1] ;
                X [j+2] = y [2] ;
             }
             for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
             {
                X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] + Lx [r] * y [2] ;
             }
          }
             // marched to next 2 columns of L
          k += 3;
       }
    }
    // FBS solve
    for(k = len_path_x - 1; k >=0;){
       j = path_x[k];
       //for(j = n-1; j >= 0; ){
 
       /* get the start, end, and length of column j */
       p = Lp [j] ;
       lnz = Lnz [j] ;
       pend = p + lnz ;
 
       /* find a chain of supernodes (up to j, j-1, and j-2) */
 
        if (s_col_FBS[k]==1)//j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
       {
 
          /* -------------------------------------------------------------- */
          /* solve with a single column of L */
          /* -------------------------------------------------------------- */
 
          double d = Lx [p] ;
          if(block_info.L->is_ll == false)
             X[j] /= d ;
          for (p++ ; p < pend ; p++)
          {
             X[j] -= Lx [p] * X [Li [p]] ;
          }
          if(block_info.L->is_ll == true)
             X [j] /=  d ;
          k--;
       }
       else if (s_col_FBS[k]==2)//lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j)
       {
          {
             double y [2], t ;
             q = Lp [j-1] ;
             double d [2] ;
             d [0] = Lx [p] ;
             d [1] = Lx [q] ;
             t = Lx [q+1] ;
             if(block_info.L->is_ll == false){
                y [0] = X [j  ] / d [0] ;
                y [1] = X [j-1] / d [1] ;
             }
             else{
                y [0] = X [j  ] ;
                y [1] = X [j-1] ;
             }
             for (p++, q += 2 ; p < pend ; p++, q++)
             {
                int i = Li [p] ;
                y [0] -= Lx [p] * X [i] ;
                y [1] -= Lx [q] * X [i] ;
             }
             if(block_info.L->is_ll == true){
                y [0] /= d [0] ;
                y [1] = (y [1] - t * y [0]) / d [1] ;
             }
             else
                y [1] -= t * y [0] ;
             X [j  ] = y [0] ;
             X [j-1] = y [1] ;
          }
          k -= 2;
       }
       else
       {
          {
             double y [3], t [3] ;
             q = Lp [j-1] ;
             r = Lp [j-2] ;
             double d [3] ;
             d [0] = Lx [p] ;
             d [1] = Lx [q] ;
             d [2] = Lx [r] ;
             t [0] = Lx [q+1] ;
             t [1] = Lx [r+1] ;
             t [2] = Lx [r+2] ;
             if(block_info.L->is_ll == false){
                y [0] = X [j]   / d [0] ;
                y [1] = X [j-1] / d [1] ;
                y [2] = X [j-2] / d [2] ;
             }
             else{
                y [0] = X [j] ;
                y [1] = X [j-1] ;
                y [2] = X [j-2] ;
             }
             for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
             {
                int i = Li [p] ;
                y [0] -= Lx [p] * X [i] ;
                y [1] -= Lx [q] * X [i] ;
                y [2] -= Lx [r] * X [i] ;
             }
             if(block_info.L->is_ll == true){
                y [0] /= d [0] ;
                y [1] = (y [1] - t [0] * y [0]) / d [1] ;
                y [2] = (y [2] - t [2] * y [0] - t [1] * y [1]) / d [2] ;
             }
             else{
                y [1] -= t [0] * y [0] ;
                y [2] -= t [2] * y [0] + t [1] * y [1] ;
             }
             X [j-2] = y [2] ;
             X [j-1] = y [1] ;
             X [j  ] = y [0] ;
          }
          k -= 3;
       }
    }
#endif
 }

// push the 8 set of internal bd nodes into node_set_x
void Circuit::push_bd_nodes(Path_Graph &pg, int&my_id){
	vector<size_t>::iterator it;
	size_t id;
	// sw direction
	size_t n_sw = internal_nodelist_sw.size();
	for(size_t i=0;i<n_sw;i++){
		id = internal_nodelist_sw[i]->rep->rid;
		it = find(pg.node_set_x.begin(), pg.node_set_x.end(), id);
		if(it == pg.node_set_x.end())
			pg.node_set_x.push_back(id);
	}

	// s direction
	size_t n_s = internal_nodelist_s.size();
	for(size_t i=0;i<n_s;i++){
		id = internal_nodelist_s[i]->rep->rid;
		it = find(pg.node_set_x.begin(), pg.node_set_x.end(), id);
		if(it == pg.node_set_x.end())
			pg.node_set_x.push_back(id);
	}

	// se direction
	size_t n_se = internal_nodelist_se.size();
	for(size_t i=0;i<n_se;i++){
		id = internal_nodelist_se[i]->rep->rid;
		it = find(pg.node_set_x.begin(), pg.node_set_x.end(), id);
		if(it == pg.node_set_x.end())
			pg.node_set_x.push_back(id);
	}

	// w direction
	size_t n_w = internal_nodelist_w.size();
	for(size_t i=0;i<n_w;i++){
		id = internal_nodelist_w[i]->rep->rid;
		it = find(pg.node_set_x.begin(), pg.node_set_x.end(), id);
		if(it == pg.node_set_x.end())
			pg.node_set_x.push_back(id);
	}

	// e direction
	size_t n_e = internal_nodelist_e.size();
	for(size_t i=0;i<n_e;i++){
		id = internal_nodelist_e[i]->rep->rid;
		it = find(pg.node_set_x.begin(), pg.node_set_x.end(), id);
		if(it == pg.node_set_x.end())
			pg.node_set_x.push_back(id);
	}

	// nw direction
	size_t n_nw = internal_nodelist_nw.size();
	for(size_t i=0;i<n_nw;i++){
		id = internal_nodelist_nw[i]->rep->rid;
		it = find(pg.node_set_x.begin(), pg.node_set_x.end(), id);
		if(it == pg.node_set_x.end())
			pg.node_set_x.push_back(id);
	}

	// n direction
	size_t n_n = internal_nodelist_n.size();
	for(size_t i=0;i<n_n;i++){
		id = internal_nodelist_n[i]->rep->rid;
		it = find(pg.node_set_x.begin(), pg.node_set_x.end(), id);
		if(it == pg.node_set_x.end())
			pg.node_set_x.push_back(id);
	}

	// ne direction
	size_t n_ne = internal_nodelist_ne.size();
	for(size_t i=0;i<n_ne;i++){
		id = internal_nodelist_ne[i]->rep->rid;
		it = find(pg.node_set_x.begin(), pg.node_set_x.end(), id);
		if(it == pg.node_set_x.end())
			pg.node_set_x.push_back(id);
	}
}

// solve transient version
bool Circuit::solve_tr_step(int &num_procs, int &my_id, MPI_CLASS &mpi_class){
	int iter = 0;	
	double diff=0;
	bool successful = false;
	double time=0;
	double t1, t2;		

	t1= MPI_Wtime();
	while( iter < MAX_ITERATION ){
		diff = solve_iteration_tr(my_id, iter, num_procs, mpi_class);
		iter++;
		// if(my_id ==0)
			//clog<<"iter, diff: "<<iter<<" "<<diff<<endl;
		if( diff < EPSILON ){
			successful = true;
			break;
		}
	}
	//if(my_id ==0)
		//clog<<"iter: "<<iter<<endl;

	t2 = MPI_Wtime();
	time = t2-t1;
	return successful;
}

void Circuit::solve_DC(int &num_procs, int &my_id, MPI_CLASS &mpi_class){
	// stores 4 boundary base into bd_base_gd
	MPI_Gather(bd_base, 8, MPI_INT, bd_base_gd, 8, MPI_INT,
		0, MPI_COMM_WORLD);
	
	MPI_Gather(internal_base, 8, MPI_INT, internal_base_gd,
		8, MPI_INT, 0, MPI_COMM_WORLD);
		
	int iter = 0;	
	double diff=0;
	bool successful = false;

	// before iteration, copy boundary nodes value to corresponding blocks
	assign_bd_internal_array(my_id);
	MPI_Gatherv(internal_x, internal_size, 
		MPI_FLOAT, 
		internal_x_g, internal_size_g, 
		internal_base_g, MPI_FLOAT, 0, 
		MPI_COMM_WORLD);
	
	// reorder boundary array according to nbrs
	if(my_id==0)	reorder_bd_x_g(mpi_class);

	double time=0;
	double t1= MPI_Wtime();
	while( iter < MAX_ITERATION ){
		diff = solve_iteration(my_id, iter, num_procs, mpi_class);
		iter++;
		// if(my_id ==0)
			// clog<<"iter, diff: "<<iter<<" "<<diff<<endl;
		if( diff < EPSILON ){
			successful = true;
			break;
		}
	}
	double t2 = MPI_Wtime();
	time = t2-t1;
	/*if(my_id==0){
		clog<<"# iter: "<<iter<<endl;
	}*/
	// get_voltages_from_block_LU_sol();
	// if(my_id==0)
		// clog<<nodelist<<endl;
}

// check if matrix is SPD or not
void Circuit::check_matrix(Matrix &A){
	size_t count = 0;
	A.merge();
	clog<<A.get_row()<<endl;
	for(size_t i=0;i<A.size();i++){
		cout<<A.Ti[i]+1<<" "<<A.Tj[i]+1<<" "<<A.Tx[i]<<endl;
		//if(A.Ti[i]!=A.Tj[i])
			//cout<<A.Tj[i]+1<<" "<<A.Ti[i]+1<<" "<<A.Tx[i]<<endl;
		//if(A.Ti[i]==A.Tj[i] && A.Tx[i]<10)
			//cout<<"diagonal:  "<<A.Ti[i]<<" "<<A.Tj[i]<<" "<<A.Tx[i]<<endl;
	}
	return;
	for(size_t i=0;i<A.size();i++){
		if(A.Ti[i]==847 || A.Tj[i]==847)
			cout<<"all: "<<A.Ti[i]<<" "<<A.Tj[i]<<" "<<A.Tx[i]<<endl;
	}
	return;
	for(size_t i=0;i<A.size();i++){
		size_t row = A.Ti[i];
		size_t col = A.Tj[i];
		if(A.Ti[i] != A.Tj[i]){// && count<1){
			count++;
			if(A.Tx[i]>=0)
			cout<<"pos: "<<A.Tx[i]<<endl;
			double row_sum=0;
			double col_sum=0;
			for(size_t j=0;j<A.size();j++){
				if(A.Ti[j]==row && A.Tj[j]==row)
					row_sum += A.Tx[j];
				if(A.Ti[j]==col && A.Tj[j]==col)
					col_sum += A.Tx[j];

			}
			clog<<"row / col sum: "<<row_sum<<" "<<col_sum<<endl;
			if(A.Tx[i]+row_sum<0){
				cout<<"error data: "<<A.Tx[i]<<" "<<row_sum<<" "<<A.Ti[i]<<" "<<A.Tj[i]<<endl;
			}
			if(A.Tx[i]+col_sum<0){
				cout<<"error data: "<<A.Tx[i]<<" "<<col_sum<<" "<<A.Ti[i]<<" "<<A.Tj[i]<<endl;
			}
		}
	}
}

// update block 4 corners
void Circuit::update_geometry(int my_id, MPI_CLASS &mpi_class){
	// compute the geometrical information for the blocks
	x_min = mpi_class.block_geo[0];
	y_min = mpi_class.block_geo[1];
	x_max = mpi_class.block_geo[2];
	y_max = mpi_class.block_geo[3];

	num_blocks = NUM_BLOCKS_X * NUM_BLOCKS_Y;
	block_vec.clear();
	for(int i=0; i<num_blocks;i++){
		Block *temp_block = new Block();
		block_vec.push_back(temp_block);
	}

	double delta_x = (x_max - x_min) / NUM_BLOCKS_X;
	double delta_y = (y_max - y_min) / NUM_BLOCKS_Y;
	double base_x = x_min;
	double base_y = y_min;
	double new_base_x;
	double new_base_y;
	// range for a block is lx <= x < ux and
	// ly <= y < uy
	for(size_t j=0;j<NUM_BLOCKS_Y;j++){
		new_base_y = base_y + j*delta_y;
		for(size_t i=0;i<NUM_BLOCKS_X;i++){
			size_t id_block = j*NUM_BLOCKS_X + i;
			new_base_x = base_x + i*delta_x;
			block_vec[id_block]->lx = new_base_x;
			block_vec[id_block]->ly = new_base_y;  
			block_vec[id_block]->ux = new_base_x + delta_x;
			block_vec[id_block]->uy = new_base_y + delta_y;

			// if(my_id==0)
				// clog<<block_vec[id_block]->lx<<" "<<block_vec[id_block]->ux<<" "<<block_vec[id_block]->ly<<" "<<block_vec[id_block]->uy<<endl;
		}
	}
}

// assign circuit nodes into blocks
void Circuit::assign_block_nodes(int my_id){
	Node *nd = NULL;
	for(size_t j=0;j<block_vec.size();j++){
		block_vec[j]->nd_GND = nd;
	}

	for(size_t i=0;i<replist.size();i++){
		nd = replist[i];
		for(size_t j=0;j<block_vec.size();j++){
			// if a node belongs to 
			// some block
			if(block_vec[j]->node_in_block(nd)){
				block_vec[j]->replist.push_back(nd);
			}
		}
	}
	// sort internal nodes of blocks
	for(size_t i=0;i<block_vec.size();i++){
		block_vec[i]->count = block_vec[i]->replist.size();
		block_vec[i]->sort_nodes();
		for(size_t j=0;j<block_vec[i]->replist.size();j++)
			// clog<<my_id<<" "<<i<<" "<<j<<" "<<*block_vec[i]->replist[j]<<endl;
		block_vec[i]->build_nd_IdMap();
	}
}

// assign circuit nets into blocks
void Circuit::assign_block_nets(int my_id){
	Net *net;
	Node *na, *nb;
	int net_flag = 0;
	// first handle internal nets of ckt
	for(int type= 0;type <= NUM_NET_TYPE; type++){	
		NetList & ns = net_set[type];
		for(size_t i=0;i<ns.size();i++){
			net = ns[i];
			// skip boundary net
			if(net->flag_bd == 1)
				continue;
			// clog<<"net: "<<*net<<endl;
			net_flag = 0;
			for(int j=0;j<block_vec.size();j++){
				net_flag = block_vec[j]->net_in_block(net);
				// if(my_id==0)
					// cout<<"block, net, flag: "<<i<<" "<<j<<" "<<*net<<" "<<net_flag<<endl;
				if(net_flag == 2)
					block_vec[j]->net_set[type].push_back(net);
				else if(net_flag ==1){
					block_vec[j]->bd_netlist.push_back(net);
				}
			}
		}
	}
	// clog<<"first stage of nets. "<<endl;
	int type  =RESISTOR;
	// then handle boundary nets of ckt
	for(size_t i=0;i<bd_netlist.size();i++){
		net = bd_netlist[i];
		net_flag = 0;
		for(int j=0;j<block_vec.size();j++){
			net_flag = block_vec[j]->net_in_block(net);

			// if(my_id==0)
				// cout<<"block, net, flag: "<<i<<" "<<j<<" "<<*net<<" "<<net_flag<<endl;
			if(net_flag == 2)
				block_vec[j]->net_set[type].push_back(net);
			else if(net_flag ==1)
				block_vec[j]->bd_netlist.push_back(net);
		}
	}
}

// build boundary netlist for circuit
void Circuit::build_bd_netlist(){
	int type = RESISTOR;
	NetPtrVector & ns = net_set[type];
	Net *net = NULL;
	for(size_t i=0;i<ns.size();i++){
		net = ns[i];
		if(net->flag_bd ==1)
			bd_netlist.push_back(net);
	}
}
