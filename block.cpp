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
//#include "umfpack.h"

Block::Block(size_t _count):
        b_ck(NULL),
	b_new_ck(NULL),
	x_old(NULL),
	//bp(NULL),
	//bnewp(NULL),
	//xp(NULL),
	x_ck(NULL),	
	count(_count),
	nodes(NULL),
	lx(-1.0), ly(-1.0),
	ux(-1.0), uy(-1.0){}

Block::~Block(){
    /* for(size_t i=0;i<boundary_netlist.size();i++)
    	delete boundary_netlist[i];
    boundary_netlist.clear();
    */
    
    A.clear();
    // free Li, Lx and so on
    delete [] Lx;
    delete [] Lp;
    delete [] Lnz;
    delete [] Li;

    nodelist.clear();
    netlist.clear();
    delete [] net_set;
    delete [] x_old;
    delete [] xp;
    delete [] bnewp;
    delete [] bp;
    delete [] bnewp_temp;
}

void Block::free_block_cholmod(cholmod_common *cm){
    cholmod_free_factor(&L, cm);
    cholmod_free_dense(&b_ck, cm);
    cholmod_free_dense(&b_new_ck, cm);
    cholmod_free_dense(&x_ck, cm);
    cholmod_free_dense(&bnew_temp, cm);
}

void Block::CK_decomp(Matrix & A, cholmod_common *cm){
	Algebra::CK_decomp(A, L, cm);
}

void Block::solve_CK(cholmod_common *cm){
	x_ck = cholmod_solve(CHOLMOD_A, L, b_new_ck, cm);
	//cholmod_solve_new(CHOLMOD_A, L, b_new_ck, x_ck, cm);
}

void Block::solve_CK_tr(cholmod_common *cm){
	x_ck = cholmod_solve(CHOLMOD_A, L, bnew_temp, cm);
	//cholmod_solve_new(CHOLMOD_A, L, b_new_ck, x_ck, cm);
}

// start cm and allocate resources
void Block::allocate_resource(){
	cm = &c;
	cholmod_start(cm);
	cm->print = 5;

	if( count == 0 ) return;
	nodes = new Node *[count];

	b_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	x_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	bp = static_cast<double*>(b_ck->x);
	xp = static_cast<double*>(x_ck->x);
	b_new_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	bnewp = static_cast<double*>(b_new_ck->x);
	bnew_temp = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	bnewp_temp = static_cast<double*>(bnew_temp->x);

	x_old = new double [count];
}

// update rhs of each block with its boundary netlist
void Block::update_rhs(double *bnewp, double *bp, int &my_id){
	size_t size = boundary_netlist.size();
	size_t k=0, l=0;

	int temp = 2;
	//b_new = b;
	for(size_t i=0;i<count;i++){
		bnewp[i] = bp[i];
		//if(my_id==1)
			//clog<<"i, bp: "<<i<<" "<<bnewp[i]<<endl;
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
			if(a->isS()!=Y){
				//if(my_id==temp){
					//clog<<k<<" "<<G<<" "<<*b<<endl;
				//}
				bnewp[k] += G * b->value;
			}
		}
		else if(b->flag_bd ==0){
			l = b->rid;
			if(b->isS()!=Y){
				//if(my_id==temp){
					//clog<<l<<" "<<G<<" "<<*a<<endl;
				//}
				bnewp[l] += G * a->value;
			}
		}
	} // end of for i
}

/////////////////////////////////////////////////////////////////
// methods for BlockInfo
void Block::sort_nodes(){
	sort(replist.begin(), replist.end(), compare_node_ptr);
}

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

// judge whether a node is within a block
bool Block::node_in_block(Node *nd){
	long x = nd->pt.x;
	long y = nd->pt.y;
	// if a node belongs to some block
	if(x>=lx && x <ux && 
		y>=ly && y<uy){
		return true;
	}
	return false;
}

// judge whether a net is within a block
// 2: internal net of a block
// 1: boundary net of a block
// 0: outside net of a block
int Block::net_in_block(Net *net){
	Node *na, *nb;
	na = net->ab[0];
	nb = net->ab[1];
	bool flag_a = false;
	bool flag_b = false;
	if(!na->is_ground()){
		flag_a = node_in_block(na);
	}
	if(!nb->is_ground()){
		flag_b = node_in_block(nb);
	}
	if(flag_a == true && flag_b == true)
		return 2;
	if(flag_a == true || flag_b == true)
		return 1;
	return 0;
}
