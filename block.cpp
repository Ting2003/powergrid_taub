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
    nd_IdMap.clear();
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

void Block::solve_CK_tr(){
	x_ck = cholmod_solve(CHOLMOD_A, L, bnew_temp, cm);
	xp = static_cast<double *>(x_ck->x);
	// for(size_t i=0;i<count;i++)
		// clog<<"i, xp: "<<i<<" "<<xp[i]<<endl;
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

	int temp = 0;
	// bnewp  = bp
	copy_array(bnewp, bp);
	//b_new = b;
	/*for(size_t i=0;i<count;i++){
		// bnewp[i] = bp[i];
		if(my_id==temp)
			clog<<"i, bp: "<<i<<" "<<bnewp[i]<<endl;
	}*/

	// for each net in this block
	for(size_t i=0;i<size;i++){
		Net * net = bd_netlist[i];
		// if(my_id==temp)
			// cout<<"net: "<<*net<<endl;
		double G = 1.0/net->value;

		Node * a = net->ab[0]->rep;
		Node * b = net->ab[1]->rep;
		// if(my_id==0)
			// cout<<"a, b: "<<*a<<" "<<*b<<endl;
	
		// if a is inside block
		//if(a->flag_bd == 0){
		if(node_in_block(a)){
			k = nd_IdMap[a];//a->rid;
			if(a->isS()!=Y){
				/*if(my_id==temp){
					cout<<"a inside block: "<<k<<" "<<*a<<" "<<G<<" "<<*b<<endl;
				}*/
				bnewp[k] += G * b->value;
				// if(my_id==temp)
					// cout<<"k, bnewp: "<<k<<" "<<bnewp[k]<<endl;
			}
		}
		else if(node_in_block(b)){
		//if(b->flag_bd ==0){
			l = nd_IdMap[b];//b->rid;
			if(b->isS()!=Y){
				/*if(my_id==temp){
					cout<<"b inside block: "<<l<<" "<<*b<<" "<<G<<" "<<*a<<endl;
				}*/
				bnewp[l] += G * a->value;
				// if(my_id==temp)
					// cout<<"l, bnewp: "<<l<<" "<<bnewp[l]<<endl;
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
	if(na->is_ground() && flag_b == true)
		return 2;
	if(nb->is_ground() && flag_a == true)
		return 2;
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
		id = nd_IdMap[node];//node->rid;
		xp[id] = replist[i]->value;
	}
}

void Block::stamp_matrix(int &my_id, MPI_CLASS &mpi_class){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			// if(my_id==0)
				//clog<<"resis net. "<<ns.size()<<endl;
			for(it=ns.begin();it!=ns.end();++it){
				Net * net = *it;
				if( net == NULL ) continue;
				assert( fzero(net->value) == false );
				stamp_resistor(my_id, *it);
			}
			break;
		case CURRENT:
			for(it=ns.begin();it!=ns.end();++it){
				stamp_current(my_id, (*it), mpi_class);
			}
			break;
		case VOLTAGE:
			// if(my_id==0)
				// clog<<"VDD net: "<<ns.size()<<endl;
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_VDD(my_id,(*it));
			}
			break;
		case CAPACITANCE:
			break;
		case INDUCTANCE:
			// if(my_id==0)
				// clog<<"induc net: "<<ns.size()<<endl;
			for(size_t i=0;i<ns.size();i++)
				stamp_inductance_dc(ns[i], my_id);
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	// then stamp boundary nets
	for(size_t i=0;i<bd_netlist.size();i++){
		stamp_bd_net(my_id, bd_netlist[i]);
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
		// A.clear();
		//if(cm->status ==1)
			//clog<<" non SPD: "<<my_id<<" "<<cm->status<<endl;
	}
	//if(my_id==0)
		//clog<<"after CK_decomp. "<<endl;
}

// =========== stamp block version of matrix =======
void Block::stamp_resistor(int &my_id, Net * net){
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
		//if(my_id==0&& j==0)
			// clog<<"net: "<<*net<<endl;
		// if boundary net of circuit
		if(net->flag_bd ==1){
			// if(my_id==0)
				// clog<<"bd net: "<<*net<<endl;
			if(!nl->is_ground() && 
				nk->flag_bd ==0 && 
				nk->isS()!=Y &&
				(nk->nbr[TOP]==NULL || 
				 nk->nbr[TOP]->type!=INDUCTANCE)){
				// stamp value into block_ids
				size_t k1 = nd_IdMap[nk]; //nk->rid;
				A.push_back(k1,k1, G);
				// if(my_id==0)
					// clog<<"nk, nl: "<<*nk<<" "<<nk->rid<<" "<<*nl<<" "<<nl->rid<<" "<<nk->is_ground()<<" "<<nl->is_ground()<<endl;

			}
		}
		// else internal net
		else if( nk->isS()!=Y && nl->isS()!=Y) {
			/*if(my_id==0){
				cout<<nk->flag_bd<<" "<<nl->flag_bd<<endl;
				//cout<<"internal net: "<<*net<<endl<<endl;
			}*/
				//clog<<"nk, nl: "<<*nk<<" "<<nk->rid<<" "<<*nl<<" "<<nl->rid<<" "<<nk->is_ground()<<" "<<nl->is_ground()<<endl;
			size_t k1 = nd_IdMap[nk];//nk->rid;
			size_t l1 = nd_IdMap[nl];//nl->rid;
			if( !nk->is_ground()&&  
          			(nk->nbr[TOP]== NULL ||
				 nk->nbr[TOP]->type != INDUCTANCE)){
				// if(my_id==0)// && nk->isS()!=Z && nk->isS()!=X && nl->isS()!=Z && nl->isS()!=Z)
				// {	//cout<<*net<<endl;
					//cout<<*nk<<" "<<k1<<endl;
					// clog<<"push + ("<<k1<<","<<k1<<","<<G<<")"<<endl;
				// }
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
				// if(my_id==0) cout<<"push -("<<*nk<<" "<<*nl<<" ("<<k1<<","<<l1<<","<<-G<<")"<<endl;
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

void Block::stamp_current(int &my_id, Net * net, MPI_CLASS &mpi_class){
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;

	// if(my_id==0)
		// clog<<"cur net: "<<*net<<endl;
	// only stamp for internal node
	if( !nk->is_ground() && nk->isS()!=Y && nk->flag_bd ==0) {
		size_t k = nd_IdMap[nk];//nk->rid;

		bp[k] += -net->value;

		// if(my_id==0) clog<<"bk: "<<k<<" "<<bp[k]<<endl;
		//pk[k] += -net->value;
	}
	if( !nl->is_ground() && nl->isS()!=Y && nl->flag_bd ==0) {
		size_t l = nd_IdMap[nl];//nl->rid;

		bp[l] += net->value;
		// if(my_id==0) clog<<"bk: "<<l<<" "<<bp[l]<<endl;
		//pl[l] +=  net->value;
	}
}

void Block::stamp_VDD(int &my_id, Net * net){
	// find the non-ground node
	Node * X = net->ab[0];
	//if(my_id==0) clog<<"net: "<<*net<<endl;
	if( X->is_ground() ) X = net->ab[1];

	if(X->rep->flag_bd ==1) return;
	// do stamping for internal node
	long id =nd_IdMap[X->rep];//X->rep->rid;

	// if(my_id==0) clog<<" stamp net: "<<*net<<endl;
	A.push_back(id, id, 1.0);
	// if(my_id==0) clog<<"push ("<<id<<", "<<id<<", 1.0)"<<endl;
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

	// if(my_id==0) clog<<"id, bp: "<<id<<", "<<bp[id]<<endl;
}

// all cores stamp dc inductance
void Block::stamp_inductance_dc(Net * net, int &my_id){
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];// nl->rid;
	G = 1./net->value;
	// if(my_id==0)
		// clog<<"induc net: "<<*net<<endl;
	if( nk->isS()!=Y && !nk->is_ground()&&nk->flag_bd==0) {
		A.push_back(k,k, 1);
		// if(my_id==0)
			// clog<<"push ("<<*nk<<" "<<k<<" "<<k<<" 1)"<<endl;
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
		// if(my_id==0)
			// clog<<"push ("<<*nl<<" "<<l<<" "<<l<<" 1)"<<endl;
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
           size_t id = nd_IdMap[q];// q->rid;
           double G = 1.0 / (*it)->value;
           
           b[id] += r->value * G;
	   // if(my_id==0)
		   // clog<<"id, new b: "<<id<<" "<<b[id]<<endl;
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
	// x_old = xp
	copy_array(x_old, xp);

	if(count<=0)
		return;
	x_ck = cholmod_solve(CHOLMOD_A, L, b_new_ck, cm);
	// solve_CK(cm);
	xp = static_cast<double *>(x_ck->x);
	/*if(my_id==0)
	for(int i=0;i<replist.size();i++)
		clog<<"xp: "<<i<<" "<<*replist[i]<<" "<<xp[i]<<endl;
	*/	
}

double Block::modify_voltage(int &my_id){
	double max_diff = 0.0;
	//if(get_name()=="VDDA") OMEGA = 1.0;
	//else OMEGA = 1.15;
	double OMEGA = 1.0;
	for(size_t i=0;i<count;i++){
		xp[i] = (1-OMEGA)*x_old[i] + OMEGA*
			xp[i];
		double vol_old = replist[i]->value;
		// update block nodes value
		replist[i]->value = xp[i];
		double diff = fabs(replist[i]->value - vol_old);
		if( diff > max_diff ) max_diff = diff;
		// if(my_id==0)
			// clog<<"rep, diff: "<<*replist[i]<<" "<<max_diff<<endl;
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

void Block::stamp_matrix_tr(int &my_id, MPI_CLASS &mpi_class, Tran &tran){	
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			for(it=ns.begin();it!=ns.end();++it){
				Net * net = *it;
				if( net == NULL ) continue;
				assert( fzero(net->value) == false );
				stamp_resistor_tr(my_id, *it);
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
				stamp_VDD_tr(my_id,(*it));
			}
			break;
		case CAPACITANCE:
			for(size_t i=0;i<ns.size();i++)
				stamp_capacitance_tr(ns[i], tran, my_id);
			break;
		case INDUCTANCE:
			for(size_t i=0;i<ns.size();i++)
				stamp_inductance_tr(ns[i], tran, my_id);
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

// only stamp resistor node connected to inductance
void Block::stamp_resistor_tr(int &my_id, Net * net){
	// no boundary nets here
	// if(net->flag_bd ==1)
		// return;

	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
	double G;	
	G = 1./net->value;

	for(size_t j=0;j<2;j++){
		Node *nk = nd[j], *nl = nd[1-j];

		// skip resistor layer net
		if(nk->isS()!=X && nl->isS()!=X) continue;
		if(nk->isS()!=X) continue;

		size_t k1 = nd_IdMap[nk];//nk->rid;
		size_t l1 = nd_IdMap[nl];//nl->rid;

		// if(my_id==0) clog<<"net: "<<*net<<endl;
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

// stamp a voltage source
void Block::stamp_VDD_tr(int &my_id, Net * net){
	// find the non-ground node
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];

	if(X->rep->flag_bd ==1) return;
	// do stamping for internal node
	long id = nd_IdMap[X->rep];//X->rep->rid;
	// A.push_back(id, id, 1.0);
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

// stamp capacitance Geq = 2C/delta_t
void Block::stamp_capacitance_tr(Net *net, Tran &tran, int &my_id){
	//clog<<"net: "<<*net<<endl;
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;
	// Geq = 2*C / delta_t
	Geq = (2*net->value) / tran.step_t;
	//net->value = Geq;
	//clog<<"C delta_t Geq: "<<net->value<<" "<<tran.step_t<<" "<<Geq<<endl;
	// Ieq = i(t) + 2*C / delta_t * v(t)

	if( nk->isS()!=Y  && !nk->is_ground() && node_in_block(nk)){//nk->flag_bd==0) {
		A.push_back(k,k, Geq);
		//clog<<"("<<k<<" "<<k<<" "<<Geq<<")"<<endl;
		if(!nl->is_ground()&& k > l){
			A.push_back(k,l,-Geq);
			//clog<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=Y && !nl->is_ground() && node_in_block(nl)){//nl->flag_bd==0) {
		A.push_back(l,l, Geq);
		//clog<<"("<<l<<" "<<l<<" "<<Geq<<")"<<endl;
		if(!nk->is_ground()&& l > k){
			A.push_back(l,k,-Geq);
			//clog<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
	}
}

// stamp inductance Geq = delta_t/(2L)
void Block::stamp_inductance_tr(Net * net, Tran &tran, int &my_id){
	//clog<<"net: "<<*net<<endl;
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;
	// Geq = delta_t / (2*L)
	Geq = tran.step_t / (2*net->value);
	//net->value = Geq;

	if( nk->isS()!=Y  && !nk->is_ground() && node_in_block(nk)){//nk->flag_bd==0) {
		// -1 is to clear formal inserted 1 at (k,k)
		A.push_back(k,k, Geq-1);
		//clog<<"("<<k<<" "<<k<<" "<<Geq-1<<")"<<endl;
		//clog<<nl->isS()<<endl;
		if(!nl->is_ground()&& nl->isS()!=Y && k>l){
			A.push_back(k,l,-Geq);
		        //clog<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=Y && !nl->is_ground() && node_in_block(nl)){//nl->flag_bd==0) {
		// -1 is to clear formal inserted 1 at (l,l)
		A.push_back(l,l, Geq-1);
		//clog<<"("<<l<<" "<<l<<" "<<Geq-1<<")"<<endl;
		if(!nk->is_ground() && nk->isS()!=Y && l>k){
			A.push_back(l,k,-Geq);
			//clog<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
	}
}

// reset b and bnew
void Block::reset_array(double *bp){
	for(size_t i=0;i<count;i++){
		bp[i] = 0;
	}
}

// make A symmetric for tran
void Block::make_A_symmetric_tr(int &my_id, Tran &tran){
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

           size_t id = nd_IdMap[q];//q->rid;
	   double G = tran.step_t / ((*it)->value*2);
           
           //b[id] += p->value * G;
           bp[id] += xp[nd_IdMap[p]]*G;//p->rid] *G;
	}
}

// stamp transient current values into rhs
void Block::stamp_current_tr(int &my_id, double &time){
	NetPtrVector & ns = net_set[CURRENT];
	for(size_t i=0;i<ns.size();i++)
		stamp_current_tr_net(ns[i], time, my_id);
}

void Block::stamp_current_tr_net(Net * net, double &time, int &my_id){
	current_tr(net, time);
	/*if(my_id==0){
		clog<<"net: "<<*net<<endl;
		// clog<<"current: "<<current<<endl;
	}*/
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground()&& nk->isS()!=Y) { 
		size_t k = nd_IdMap[nk];//nk->rid;
		//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
		bp[k] += -net->value;//current;
		// if(my_id==0)
			// clog<<"time, k, b: "<<time<<" "<<k<<" "<<bp[k]<<endl;
	}
	if( !nl->is_ground() && nl->isS()!=Y) {
		size_t l = nd_IdMap[nl];//nl->rid;
		//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
		bp[l] +=  net->value;// current;
		// if(my_id==0)
			// clog<<"time, l, b: "<<time<<" "<<l<<" "<<bp[l]<<endl;
	}
}

void Block::clear_A(){
	A.clear();
}

void Block::CK_decomp(){
	Algebra::CK_decomp(A, L, cm);
}

void Block::copy_vec(double *bnewp, double *bp){
	for(size_t i=0;i<count;i++){
		bnewp[i] = bp[i];
	}
}

// update rhs by transient nets
void Block::modify_rhs_tr_0(double * b, double *x, int &my_id){
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

// stamp transient current values into rhs
void Block::stamp_current_tr_1(double &time){
	NetPtrVector & ns = net_set[CURRENT];
	for(size_t i=0;i<ns.size();i++)
		stamp_current_tr_net_1(ns[i], time);
}

void Block::stamp_current_tr_net_1(Net * net, double &time){
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
			size_t k = nd_IdMap[nk];//nk->rid;
			//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
			//clog<<"time, k, b bef: "<<time<<" "<<k<<" "<<b[k]<<endl;
			bnewp[k] += -diff;//current;
			bp[k] = bnewp[k];
			//clog<<"time, k, b: "<<time <<" "<<k<<" "<<b[k]<<endl;
		}
		if( !nl->is_ground() && nl->isS()!=Y) {
			size_t l = nd_IdMap[nl];//nl->rid;
			//clog<<"time, l, b bef: "<<time<<" "<<l<<" "<<b[l]<<endl;
			//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
			bnewp[l] +=  diff;// current;
			bp[l] = bnewp[l];
			//clog<<"time, l, b: "<<time<<" "<<l<<" "<<b[l]<<endl;
		}
	}
}

// update rhs by transient nets
void Block::modify_rhs_tr(double * b, double *x){
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

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void Block::modify_rhs_c_tr(Net *net, double * rhs, double *x){
	double temp = 0;
	//clog<<"c net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
        
	// nk point to Z node
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;

	Net *r = nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node

	size_t id_a = nd_IdMap[a];//a->rid;
	size_t id_b = nd_IdMap[b];//b->rid;
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

// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void Block::modify_rhs_l_tr(Net *net, double *rhs, double *x){
	//clog<<"l net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	//if(nk->isS() !=X){ 
		//swap<Node*>(nk, nl);
	//}
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;
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
	size_t id_a = nd_IdMap[a];//a->rid;
	size_t id_b = nd_IdMap[b];//b->rid;
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

// decide transient step current values
void Block::current_tr(Net *net, double &time){
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
void Block::modify_rhs_c_tr_0(Net *net, double * rhs, double *x, int &my_id){
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
	size_t k = nd_IdMap[nk];//nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;
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

	size_t id_a = nd_IdMap[a];//a->rid;
	size_t id_b = nd_IdMap[b];//b->rid;
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
#if 0
        pg.node_set_x.push_back(k);
        if(!nl->is_ground()) {
              //clog<<*nl<<" "<<l<<endl;
           pg.node_set_x.push_back(l);
        }
        else if(!b->is_ground()){
              //clog<<*b<<" "<<id_b<<endl;
           pg.node_set_x.push_back(id_b);
        }
#endif
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
// Ieq = i(t) + delta_t / (2*L) *v(t)
void Block::modify_rhs_l_tr_0(Net *net, double *rhs, double *x, int &my_id){
	//clog<<"l net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	if(nk->isS() !=X){ 
		swap<Node*>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
	size_t k = nd_IdMap[nk];//]nk->rid;
	size_t l = nd_IdMap[nl];//nl->rid;
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
	size_t id_a = nd_IdMap[a];//a->rid;
	size_t id_b = nd_IdMap[b];//b->rid;
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
        pg.node_set_x.push_back(k);
        pg.node_set_x.push_back(id_b);
#endif
	Ieq  = i_t + temp;
	if(nk->isS() !=Y && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD circuit
		// if(my_id==0)
			// cout<<"k, Ieq, rhs: "<<k<<" "<<Ieq<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=Y && !nl->is_ground()){
		 rhs[l] += -Ieq; // VDD circuit

		// if(my_id==0)
			// cout<<"l, Ieq, rhs: "<<l<<" "<<Ieq<<" "<<rhs[l]<<endl;
	}
}

void Block::stamp_bd_net(int my_id, Net *net){
	if(net->type != RESISTOR)
		return;

	Node *na= NULL;
	Node *nb = NULL;

	na = net->ab[0];
	nb = net->ab[1];
	size_t id;
	if(node_in_block(na)){
		id = nd_IdMap[na];
		A.push_back(id, id, 1.0/net->value);
		// if(my_id==0)
			// clog<<"bd net: "<<id<<" "<<id<<" "<<1.0/net->value<<endl;
	}
	else if(node_in_block(nb)){
		id = nd_IdMap[nb];
		A.push_back(id, id, 1.0/net->value);
		// if(my_id==0)
			// clog<<"bd net: "<<id<<" "<<id<<" "<<1.0/net->value<<endl;
	}
}
