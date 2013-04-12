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

    replist.clear();
    bd_netlist.clear();
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
	size_t size = bd_netlist.size();
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
		Net * net = bd_netlist[i];
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

void Block::sort_nodes(){
	sort(replist.begin(), replist.end(), compare_node_ptr);
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

// 1. copy node voltages from the circuit to a Vec
//    from = true then copy circuit to x
//    else copy from x to circuit
// 2. map block voltage into global_index
void Block::copy_node_voltages_block(){
	size_t id;
	// copy node voltages from nodelist
	for(size_t i=0;i<replist.size();i++){
		Node *node = replist[i];
		id = node->rid;
		xp[id] = replist[i]->value;
	}
}

void Block::stamp_matrix(int &my_id, MPI_CLASS &mpi_class){
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
	 make_A_symmetric(bp, my_id);
	
	A.set_row(count);
	//if(my_id==0)
		//check_matrix(A);
		//cout<<"before CK_decomp. "<<endl;
	if(count >0){
		CK_decomp(A, cm);
		//if(cm->status ==1)
			//clog<<" non SPD: "<<my_id<<" "<<cm->status<<endl;
	}
	//if(my_id==0)
		//clog<<"after CK_decomp. "<<endl;
}

// =========== stamp block version of matrix =======
void Block::stamp_block_resistor(int &my_id, Net * net, Matrix &A){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
	
	double G;	
	G = 1./net->value;
	
	// if(net->flag_bd ==1)
		// bd_netlist.push_back(net);

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
}

void Block::stamp_block_current(int &my_id, Net * net, MPI_CLASS &mpi_class){
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;

	// only stamp for internal node
	if( !nk->is_ground() && nk->isS()!=Y && nk->flag_bd ==0) {
		size_t k = nk->rid;

		bp[k] += -net->value;
		//pk[k] += -net->value;
	}
	if( !nl->is_ground() && nl->isS()!=Y && nl->flag_bd ==0) {
		size_t l = nl->rid;

		bp[l] += net->value;
		//if(my_id==1) clog<<"bk: "<<l<<" "<<block_info.bp[l]<<endl;
		//pl[l] +=  net->value;
	}
}

void Block::stamp_block_VDD(int &my_id, Net * net, Matrix &A){
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
		assert( feqn(1.0, bp[id]) );
		bp[id] = net->value;
		//q[id] = net->value;	    // modify it
	}
	else{
		bp[id] += net->value;
		//q[id] += net->value;
	}
}

// all cores stamp dc inductance
void Block::stamp_inductance_dc(Matrix & A, Net * net, int &my_id){
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
			bp[k] = bp[l];
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
		bp[l] = bp[k];

		//clog<<"("<<l<<" "<<l<<" "<<1<<")"<<endl;
		//clog<<"("<<l<<" "<<k<<" "<<-1<<")"<<endl;
	}
}

void Block::make_A_symmetric(double *b, int &my_id){
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

void Block::copy_array(double *x_old, double *xp){
	for(size_t j=0;j<count;j++)
		x_old[j] = xp[j];	
}

// update rhs and solve CK
void Block::solve_CK_DC(int my_id){
	// new rhs store in bnewp and solve
	update_rhs(bnewp, bp, my_id);
	// x_old stores old solution
	copy_array(x_old, xp);

	if(count<=0)
		return;
	x_ck = cholmod_solve(CHOLMOD_A, L, b_new_ck, cm);
	// solve_CK(cm);
	xp = static_cast<double *>(x_ck->x);
	/*if(my_id==0 && iter<3)
	for(int i=0;i<replist.size();i++)
		clog<<"xp: "<<i<<" "<<*replist[i]<<" "<<block_info.xp[i]<<endl;*/
}

double Block::modify_voltage(int &my_id){
	double max_diff = 0.0;
	//if(get_name()=="VDDA") OMEGA = 1.0;
	//else OMEGA = 1.15;
	double OMEGA = 1.0;
	for(size_t i=0;i<count;i++){
		xp[i] = (1-OMEGA)*x_old[i] + OMEGA*
			xp[i];
		// update block nodes value
		replist[i]->value = xp[i];
		double diff = fabs(x_old[i] - xp[i]);
		if( diff > max_diff ) max_diff = diff;
	}
	return max_diff;
}

// build block node id map
void Block::build_nd_IdMap(){
	pair<Node *, size_t> id_pair;
	for(size_t i=0;i<replist.size();i++){
		id_pair.first = replist[i];
		id_pair.second = i;
		nd_IdMap.insert(id_pair);
	}
}
