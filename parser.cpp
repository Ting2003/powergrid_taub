// ----------------------------------------------------------------//
// Filename : parser.cpp
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// implementation file of parser.h
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include "util.h"
#include "parser.h"
#include "mpi.h"
using namespace std;

// store the pointer to circuits
Parser::Parser(vector<Circuit*> * ckts):n_layer(0), p_ckts(ckts),
	layer_in_ckt(vector<int>(MAX_LAYER)){
}

Parser::~Parser(){ }

// _X_n2_19505_20721 
// X:
// n2: layer 2
// 19505 20721: coordinate
void Parser::extract_node(char * str, Node & nd){
	//static Node gnd(string("0"), Point(-1,-1,-1));
	if( str[0] == '0' ) {
		nd.name="0";
		nd.pt.set(-1,-1,-1);
		return;
	}

	long z, y, x;
	int flag = -1;
	char * chs;
	char * saveptr;
	char l[MAX_BUF];
	strcpy(l, str);
	const char * sep = "_n";
	chs = strtok_r(l, sep, &saveptr); // initialize
	// for transient, 'Y' is the VDD source node
	if( chs[0] == 'X' || chs[0]== 'Y' || chs[0]== 'Z' ){
		flag = chs[0]-'X';
		chs = strtok_r(NULL, sep, &saveptr);
	}
	
	z = atol(chs);
	chs = strtok_r(NULL, sep, &saveptr);
	x = atol(chs);
	chs = strtok_r(NULL, sep, &saveptr);
	y = atol(chs);
	
	nd.name.assign(str);
	nd.pt.set(x,y,z);
	nd.flag = flag;
	nd.rid = -1;
}

// given a line, extract net and node information
void Parser::insert_net_node(char * line, int &my_id, MPI_CLASS &mpi_class){
	char *chs, *saveptr;
	const char* sep = " (),\n";

	static char sname[MAX_BUF];
	static char sa[MAX_BUF];
	static char sb[MAX_BUF];
	static Node nd[2];
	Node * nd_ptr[2];	// this will be set to the two nodes found
	double value;
	int count=0;

	sscanf(line, "%s %s %s %lf", sname, sa, sb, &value);

	extract_node(sa, nd[0]);
	extract_node(sb, nd[1]);

	int layer;
	if( nd[0].is_ground() ) 
		layer = nd[1].get_layer();
	else
		layer = nd[0].get_layer();

	int ckt_id = layer_in_ckt[layer];
	Circuit * ckt = (*p_ckts)[ckt_id];
	
	for(int i=0;i<2;i++){
		if ( nd[i].is_ground() ){
			nd_ptr[i] = ckt->nodelist[0]; // ground node
		}
		else if ( (nd_ptr[i] = ckt->get_node(nd[i].name) ) == NULL ){
			// create new node and insert
			nd_ptr[i] = new Node(nd[i]); // copy constructor
			nd_ptr[i]->rep = nd_ptr[i];  // set rep to be itself
			count = cpr_nd_block(nd_ptr[i], mpi_class.block_geo, my_id);
			if (count==1) // internal node
				ckt->add_node(nd_ptr[i]);
			else
				nd_ptr[i]->flag_bd = 1;

			if( nd_ptr[i]->isS()==Y)	     // determine circuit type
				ckt->set_type(WB);
		}
	}

	// if there is a cap with it, need to modify
	if((line[0]=='r' || line[0] =='R') && 
		!(nd[0].pt.x == nd[1].pt.x && 
		nd[0].pt.y == nd[1].pt.y)){
		// add node into bd and internal vector
		add_node_inter(nd_ptr[0], nd_ptr[1], 
			mpi_class, ckt, my_id);
	}
	
	NET_TYPE net_type = RESISTOR;
	// find net type
	switch(sname[0]){
	case 'r': // resistor
	case 'R':
		net_type = RESISTOR;
		break;
	case 'v': // VDD
	case 'V':
		net_type = VOLTAGE;
		break;
	case 'i': // current
	case 'I':
		net_type = CURRENT;
		break;
	case 'c': // capacitance
	case 'C':
		net_type = CAPACITANCE;
		break;
	case 'l':
	case 'L':
		net_type = INDUCTANCE;
		break;
	default:
		report_exit("Invalid net type!\n");
		break;
	}

	// create a Net
	Net * net = new Net(net_type, value, nd_ptr[0], nd_ptr[1]);

	if(net_type == CURRENT){
		net->tr = new double [7];
		// assign pulse paramter for pulse input
		chs = strtok_r(line, sep, &saveptr);
		for(int i=0;i<3;i++)
			chs = strtok_r(NULL, sep, &saveptr);
		if(chs != NULL)
			chs = strtok_r(NULL, sep, &saveptr);
		if(chs != NULL){
			chs = strtok_r(NULL, sep, &saveptr);
			// V1
			net->tr[0] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// V2
			net->tr[1] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// TD
			net->tr[2] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// Tr
			net->tr[3] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// Tf
			net->tr[4] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// PW
			net->tr[5] = atof(chs);
			chs = strtok_r(NULL, sep, &saveptr);
			// Period
			net->tr[6] = atof(chs);
		}
	}
	// trick: when the value of a resistor via is below a threshold,
	// treat it as a 0-voltage via
	//if( Circuit::MODE == (int)IT ) {
		try_change_via(net);
	//}

	// insert this net into circuit
	ckt->add_net(net);

	// IMPORTANT: set the relationship between node and net
	// update node voltage if it is an X node
	// set node to be X node if it connects to a voltage source
	update_node(net);
}

// Given a net with its two nodes, update the connection information
// for thet two nodes
void Parser::update_node(Net * net){
	// first identify their connection type:
	// 1. horizontal/vertical   2. via/VDD 3. current
	//
	// swap *a and *b so that a is:
	// WEST   for horizontal
	// SOUTH  for vertical
	// BOTTOM for via / XVDD
	// ground node for CURRENT
	Node *a=net->ab[0], *b=net->ab[1];
	//cout<<"setting "<<net->name<<" nd1="<<nd1->name<<" nd2="<<nd2->name<<endl;
	if(net->type == CAPACITANCE){
		// make sure a is the Z node, b is ground
		if(a->isS() != Z) swap<Node*>(a,b);
		// only needs single dir nbr net for index
		// TOP for resistance
		// BOTTOM for capacitance
		a->set_nbr(BOTTOM, net);
	}
	else if(net->type == INDUCTANCE){
		// a is Y, b is X
		if(a->isS()== X) swap<Node*>(a,b);
		a->set_nbr(BOTTOM, net);
		b->set_nbr(TOP, net);
	}
	// resistance type net with special nodes
	else if(net->type == RESISTOR && a->isS()!= -1 
		&& b->isS() == -1){
		if(a->isS()== X){
			a->set_nbr(BOTTOM, net);
			b->set_nbr(TOP, net);
		}
		else if(a->isS() == Z){
			// bottom for z node has been taken by 
			// current
			a->set_nbr(TOP, net);
		}
	}
	else if(net->type == RESISTOR && b->isS()!= -1 && a->isS() ==-1){
		if(b->isS()== X){
			b->set_nbr(BOTTOM, net);
			a->set_nbr(TOP, net);
		}
		else if(b->isS() == Z){
			b->set_nbr(TOP, net);
		}
	}
	else if( a->get_layer() == b->get_layer() && a->pt != b->pt ){
		// horizontal or vertical resistor in the same layer
		int layer = a->get_layer();
		if( a->pt.y == b->pt.y ){// horizontal
			if(a->pt.x > b->pt.x) swap<Node*>(a,b);
			a->set_nbr(EAST, net);
			b->set_nbr(WEST, net);
			Circuit::layer_dir[layer] = HR;
		}
		else if( a->pt.x == b->pt.x ){// vertical
			if(a->pt.y > b->pt.y) swap<Node*>(a,b);
			a->set_nbr(NORTH, net);
			b->set_nbr(SOUTH, net);
			Circuit::layer_dir[layer] = VT;
		}
		else
			report_exit("Diagonal net\n");
	}
	else if( //fzero(net->value) && 
		 !a->is_ground() &&
		 !b->is_ground() ){// this is Via (Voltage or Resistor )
		if( a->get_layer() > b->get_layer() ) swap<Node*>(a,b);
		a->set_nbr(TOP, net);
		b->set_nbr(BOTTOM, net);
	}
	else if (net->type == VOLTAGE){// Vdd Voltage
		// one is X node, one is ground node
		// Let a be X node, b be another
		if( a->is_ground() ) swap<Node*>(a,b);
		a->flag = Y;		// set a to be X node
		a->set_nbr(TOP, net);	// X -- VDD -- Ground
		a->set_value(net->value);
	}
	else{// if( net->type == CURRENT ){// current source
		// let a be ground node
		if( !a->is_ground() ) swap<Node*>(a,b);
		b->set_nbr(BOTTOM, net);
	}
}

// parse the file and create circuits
int Parser::create_circuits(vector<CKT_LAYER > &ckt_name_info){
	int layer, n_circuit=0;

	string prev_ckt_name("");
	string name_string;
	Circuit * p_last_circuit=NULL;
	// now read filename.info to create circuits (they are SORTED)
	for(size_t i=0;i<ckt_name_info.size();i++){
		name_string = ckt_name_info[i].name;
		layer = ckt_name_info[i].layer;
		//cout<<name_string<<":"<<layer<<endl;
		// compare with previous circuit name 
		//name_string.assign(name);
		if( prev_ckt_name == "" ||
		    name_string != prev_ckt_name ){
			Circuit * circuit = new Circuit(name_string);
			(*p_ckts).push_back(circuit);
			++n_circuit;
			prev_ckt_name = name_string;
			p_last_circuit = circuit;
		}

		p_last_circuit->layers.push_back(layer);

		// note that initial size may not be accurate
		if( layer > (int)layer_in_ckt.size()-1 ) 
			layer_in_ckt.resize(layer+10); // 10 can be a arbitrary num.

		layer_in_ckt[layer] = n_circuit-1; // map layer id to circuit id
		this->n_layer++;
	}
	
	// now we know the correct number of layers
	layer_in_ckt.resize(this->n_layer);
	Circuit::layer_dir.resize(this->n_layer);

	return n_circuit;
}

// parse the file
// Note: the file will be parsed twice
// the first time is to find the layer information
// and the second time is to create nodes
void Parser::parse(int &my_id, char * filename, MPI_CLASS &mpi_class, Tran &tran, int num_procs){	
	int MPI_Vector;
	int count =2;
	int lengths[2] = {10, 1};
	MPI_Aint offsets[2] = {0, sizeof(char)*10};
	MPI_Datatype types[3]={MPI_CHAR, MPI_INT};
	MPI_Type_struct(count, lengths, offsets, types, 
			&MPI_Vector);
	MPI_Type_commit(&MPI_Vector);

	this->filename = filename;

	// processor 0 will extract layer info
	// and bcast it into other processor
	vector<CKT_LAYER >ckt_name_info;
	if(my_id==0)
		extract_layer(my_id, ckt_name_info, mpi_class, tran);
	// broadcast info for transient 
	
	int ckt_name_info_size = ckt_name_info.size();
	
	MPI_Bcast(&ckt_name_info_size, 1, MPI_INT, 0,
			MPI_COMM_WORLD);
	if(my_id!=0) ckt_name_info.resize(ckt_name_info_size);

	MPI_Bcast(&ckt_name_info[0], ckt_name_info_size, 
			MPI_Vector, 0, MPI_COMM_WORLD);

	// first time parse:
	create_circuits(ckt_name_info);

	build_block_geo(my_id, mpi_class, tran, num_procs);

	MPI_Barrier(MPI_COMM_WORLD);

	// temporary comment second parse	
	second_parse(my_id, mpi_class, tran, num_procs);
}

void Parser::build_block_geo(int &my_id, MPI_CLASS &mpi_class, Tran &tran, int num_procs){
	// block info in cpu 0
	mpi_class.geo = new float[mpi_class.X_BLOCKS *mpi_class.Y_BLOCKS *4];
	// block info in other cpus
	mpi_class.block_geo = new float [4];

	// update block geometry
	if(my_id==0) set_block_geometry(mpi_class.geo, mpi_class);
	int total_blocks = 4 * mpi_class.X_BLOCKS * mpi_class.Y_BLOCKS;
	MPI_Bcast(mpi_class.geo, total_blocks, MPI_FLOAT, 0, MPI_COMM_WORLD);

	MPI_Bcast(&mpi_class.len_ovr_x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mpi_class.len_ovr_y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(mpi_class.geo, 4, MPI_FLOAT, 
		mpi_class.block_geo, 4, MPI_FLOAT, 
		0, MPI_COMM_WORLD);
	
	if(my_id==0){
		net_to_block(mpi_class.geo, mpi_class, tran, num_procs, my_id);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//if(my_id==3)
		//clog<<"lx, ly, ux, uy: "<<mpi_class.block_geo[0]<<" "<<mpi_class.block_geo[1]<<" "<<mpi_class.block_geo[2]<<" "<<mpi_class.block_geo[3]<<endl;
	// stores the geo boundary line for internal node:
	// e, w, n, s
	mpi_class.block_geo_origin = new float [4];
	
	mpi_class.set_geo_origin(mpi_class);
}

void Parser::second_parse(int &my_id, MPI_CLASS &mpi_class, Tran &tran, int num_procs){
	char  buff[100];
	FILE *f = NULL;

	int color = 0;
	int block_size = mpi_class.block_size;
	//if(block_size==0) return;
	//for(int i=0;i<block_size;i++){
		sprintf(buff, "./INPUT_FILE/netlist_%d%d.txt", color, my_id);
		f = fopen(buff, "r");
		if(f==NULL) report_exit("Input file not exist!\n");
	//}
	
	char line[MAX_BUF];
	char type;

	tran.nodes.clear();
	while( fgets(line, MAX_BUF,f)!=NULL){
		type = line[0];
		//clog<<line<<endl;
		switch(type){
			case 'r': // resistor
			case 'R':
			case 'v': // VDD
			case 'V':
			case 'i': // current
			case 'I':
			case 'c':
			case 'C':
			case 'l':
			case 'L':
				insert_net_node(line, my_id, mpi_class);
				break;
			case '.': // parse tran nodes: need to write
				 block_parse_dots(line, tran, my_id); 
				 break;
			case '*': // comment
			case ' ':
			case '\n':
				break;
			default:
				printf("Unknown input line: ");
				report_exit(line);
				break;
		}
	}
	fclose(f);
	// release map_node resource
	for(size_t i=0;i<(*p_ckts).size();i++){
		Circuit * ckt = (*p_ckts)[i];
		ckt->map_node.clear();
	}
}// end of parse

int Parser::get_num_layers() const{ return n_layer; }

// done by core 0
void Parser::parse_dot(char *line, Tran &tran){
	char *chs;
	char *saveptr;
	char sname[MAX_BUF];
	Node_TR_PRINT item;
	const char *sep = "= v() \n";
	switch(line[1]){
		case 't': // transient steps
			/*sscanf(line, "%s %lf %lf", sname, 
				&tran.step_t, &tran.tot_t);
			tran.isTran = 1; // do transient ana;*/
			//clog<<"step: "<<tran.step_t<<" tot: "<<tran.tot_t<<endl;
			break;
		case 'w': // output length
			/*chs = strtok_r(line, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			tran.length = atoi(chs);*/
			//clog<<"out len: "<<tran.length<<endl;
			break;
		case 'p': // print
			chs = strtok_r(line, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			while(chs != NULL){
				chs = strtok_r(NULL, sep, &saveptr);
				if(chs == NULL) break;
				item.name = chs;
				tran.nodes.push_back(item);
				// disribute nodes into cores
			};
			break;
		default: 
			break;
	}
}

int Parser::extract_layer(int &my_id, 
		vector<CKT_LAYER > &ckt_layer_info,
		MPI_CLASS &mpi_class,
		Tran &tran){
	char line[MAX_BUF];
	char word[MAX_BUF];
	string word_s;
	char name[10];
	static char sname[MAX_BUF];
	static char sa[MAX_BUF];
	static char sb[MAX_BUF];
	static Node nd[2];
	double value;
	int i=0;
	long x_max=0;
	long x_min=0;
	long y_max=0;
	long y_min=0;

	// only processor 0 will extract layer info
	if(my_id!=0) return 0;

	FILE *f;
	f = fopen(filename, "r");
	if(f==NULL) report_exit("Input file not exist!\n");

	CKT_LAYER ckt_name_layer;

	while(fgets(line, MAX_BUF, f)!=NULL){
		if(line[0]=='*'){
			// copy the entire line into stringstream
			stringstream ss;
			ss<< line;
			//if(my_id==0) clog<<"line: "<<line<<endl;
			int word_count = 0;
			while(ss.getline(word, 10, ' ')){
				if(word_count ==1){
					word_s = word;
					if(word_s !="layer:"){
						break;
					}
				}
				if(word_count==2){
					stringstream ss_1;
					ss_1<<word;
					int vdd_count=0;
					while(ss_1.getline(name, 10,',')){
						if(vdd_count==1){
							// extract ckt->name
							strcpy(ckt_name_layer.name, name);
						}
						vdd_count++;
					}
				}
				// extract layer number
				if(word_count==4){
					ckt_name_layer.layer = atoi(word);
					ckt_layer_info.push_back(ckt_name_layer);
					//fprintf(fp, "%s \n", word);
				}
				word_count++;
			}
		}
		else if(line[0]!='.'){
			// find grid boundary x and y
			sscanf(line, "%s %s %s %lf", sname,sa,sb, &value);
			if( sa[0] != '0' ) 
				extract_node(sa, nd[0]);

			if( sb[0] != '0' ) 
				extract_node(sb, nd[1]);

			for(i=0;i<2;i++){
				if(nd[i].pt.x > x_max) 
					x_max = nd[i].pt.x;
				if(nd[i].pt.x >0 && nd[i].pt.x <x_min)
					x_min = nd[i].pt.x;
				if(nd[i].pt.y > y_max)
					y_max = nd[i].pt.y;
				if(nd[i].pt.y>0 && nd[i].pt.y <y_min)
					y_min = nd[i].pt.y;
			}
		}else if(line[0] == '.'){
			// parse_dot by core 0
			parse_dot(line, tran);
		}
	}
	mpi_class.x_max = x_max;
	mpi_class.x_min = x_min;
	mpi_class.y_max = y_max;
	mpi_class.y_min = y_min;
	fclose(f);
	// sort resulting vector by the ckt name
	sort(ckt_layer_info);
	return 0;
}
// sort ckt according to its name and layer
// ckt name decrease, layer rising
bool Parser::sort(vector<CKT_LAYER > &a){
	CKT_LAYER tmp;
	// sort according to circuit name
	for(size_t i=0;i< a.size()-1;i++){
		int minIndex = i;
		for(size_t j = i+1;j< a.size();j++)
			if(strcmp(a[j].name, a[minIndex].name)>0)
				minIndex = j;
		if(minIndex !=i){
			tmp = a[i];
			a[i] = a[minIndex];
			a[minIndex] = tmp;
		}	
	}
	return true;
}

void Parser::set_block_geometry(float *geo, MPI_CLASS &mpi_class){
	double x, y;
	x = (double)(mpi_class.x_max-mpi_class.x_min+0.5) 
		/ mpi_class.X_BLOCKS;
	y = (double)(mpi_class.y_max-mpi_class.y_min+0.5) 
		/ mpi_class.Y_BLOCKS;
	//if( fzero(x) ) x = 1.0;
	//if( fzero(y) ) y = 1.0;
	double len_per_block_x = x;
	double len_per_block_y = y;
	double len_ovr_x = x * mpi_class.overlap_ratio;
	double len_ovr_y = y * mpi_class.overlap_ratio;
	mpi_class.len_per_block_x = x;
	mpi_class.len_per_block_y = y;
	mpi_class.len_ovr_x = len_ovr_x;
	mpi_class.len_ovr_y = len_ovr_y;
	clog<<"len_per_x, len_per_y: "<<x<<" "<<y<<endl;
	clog<<"len_ovr_x, len_ovr_y: "<<len_ovr_x
		<<" / "<<len_ovr_y<<endl;

	size_t bid = 0;
	// update block 4 corners
	for(size_t y=0;y<mpi_class.Y_BLOCKS;y++){
		for(size_t x=0;x<mpi_class.X_BLOCKS;x++){
			bid = y * mpi_class.X_BLOCKS + x;
			// lx
			geo[4*bid] = x * len_per_block_x - len_ovr_x;
			// ly
			geo[4*bid+1] = y * len_per_block_y - len_ovr_y;
			// ux
			geo[4*bid+2] = (x+1) * len_per_block_x + len_ovr_x;
			// uy
			geo[4*bid+3] = (y+1) * len_per_block_y + len_ovr_y;
		}
	}
}

// done by processor 0
void Parser::net_to_block(float *geo, MPI_CLASS &mpi_class, Tran &tran, int num_procs, int &my_id){
	static char sname[MAX_BUF];
	static char sa[MAX_BUF];
	static char sb[MAX_BUF];
	static Node nd[2];
	double value;

	int temp = 0;

	FILE *f;
	f = fopen(this->filename, "r");
	if(f==NULL) report_exit("Input file not exist!\n");

	// handles the long print sentence
	int temp_buf = 10000000;
	char line[temp_buf];
	string line_s;

	int color = 0;
	vector<FILE *> of;
	int num_blocks  = mpi_class.X_BLOCKS * mpi_class.Y_BLOCKS;
	clog<<"num_blocks. "<<num_blocks<<endl;
	InitialOF(of, num_procs, color);//num_blocks, color);
	//clog<<"after initial ofs. "<<endl;

	int count_1 = 0, count_2 = 0;
	while( fgets(line, temp_buf, f)!=NULL){
		if(line[0]=='r' || line[0] =='R' ||
		   line[0]=='v' || line[0] =='V' ||
		   line[0]=='i' || line[0]=='I' ||
		   line[0]=='c' || line[0] == 'C' ||
		   line[0]=='l' || line[0] == 'L'){
			//clog<<line<<endl;
			sscanf(line, "%s %s %s %lf", 
					sname, sa, sb, &value);
			if( sa[0] == '0' ){ nd[0].pt.set(-1,-1,-1); }
			else extract_node(sa, nd[0]);
			if( sb[0] == '0' ){ nd[1].pt.set(-1,-1,-1); }
			else	extract_node(sb, nd[1]);
		
			for(int i=0;i<num_blocks;i++){
				// at least one node is inside block
				count_1 = cpr_nd_block(nd[0], geo, i);
				count_2 = cpr_nd_block(nd[1], geo, i);
				
				if(nd[0].is_ground() && count_2 ==1)
					temp = fprintf(of[i], "%s", line);
				else if(nd[1].is_ground() && count_1 ==1)
					temp = fprintf(of[i], "%s", line);
				// write all voltage sources
				else if((count_1 + count_2 ==2)){
					temp = fprintf(of[i], "%s", line);
					//clog<<temp<<" "<<i<<" "<<line<<endl;
				}
			}
		}
		else{
			switch(line[1]){
				case 't':
				case 'w':
				case 'e':// for all procs
					for(int i=0;i<num_procs;i++){
						fprintf(of[i], "%s", line);
					}
					break;
				case 'p': // print
					write_print(tran, of, mpi_class, line);
					break;
				default:
					break;

			}
		}
	}
	free(line);
	// finally print end file symbol (not need to)
	//clog<<"finish output. "<<endl;
	fclose(f);
	//clog<<"close original file. "<<endl;
	for(int i=0;i<num_procs;i++){
		fclose(of[i]);
	}
	//clog<<"close all output file. "<<endl;
	of.clear();
}

int Parser::cpr_nd_block(Node &nd, float *geo, int &bid){
	if(nd.pt.x >= geo[4*bid] &&
	   nd.pt.x <= geo[4*bid+2] &&
	   nd.pt.y >= geo[4*bid+1] &&
	   nd.pt.y <= geo[4*bid+3]){
		return 1;
	}
	else return 0;
}

int Parser::cpr_nd_block(Node *nd, float *geo, int &bid){
	if(nd->pt.x >= geo[0] &&
	   nd->pt.x <= geo[2] &&
	   nd->pt.y >= geo[1] &&
	   nd->pt.y <= geo[3]){
		return 1;
	}
	else return 0;
}

int Parser::cpr_nd_block(Node *nd, float &lx, float &ly, float &ux, float &uy){
	if(nd->pt.x >= lx&&
	   nd->pt.x <= ux &&
	   nd->pt.y >= ly &&
	   nd->pt.y <= uy){
		return 1;
	}
	else {
		return 0;
	}
}

void Parser::add_node_inter(Node *nd_0, Node *nd_1, 
	MPI_CLASS &mpi_class, Circuit *ckt, int &my_id){
	int bx, by;
	int bid_nbr;
	float lx, ly, ux, uy;
	float lx_0, ly_0, ux_0, uy_0;
	int count_1, count_2;
	int count_10, count_20;

	by = my_id / mpi_class.X_BLOCKS;
	bx = my_id % mpi_class.X_BLOCKS;

	find_bound_line(my_id, mpi_class, lx_0, ly_0,
			ux_0, uy_0);
	
	count_10 = cpr_nd_block(nd_0, lx_0, ly_0, ux_0, uy_0);
	count_20 = cpr_nd_block(nd_1, lx_0, ly_0, ux_0, uy_0);
	
	// sw block
	if((by>=1 && bx>=1)){
		bid_nbr = my_id - mpi_class.X_BLOCKS - 1;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_sw, ckt->internal_nodelist_sw, 
		nd_0, nd_1, count_10, count_20);
	}

	// s block
	if(by>=1){
		bid_nbr = my_id - mpi_class.X_BLOCKS;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_s, ckt->internal_nodelist_s, 
		nd_0, nd_1, count_10, count_20);
	}
	
	// se block
	if((by>=1 && bx<mpi_class.X_BLOCKS-1)){
		bid_nbr = my_id - mpi_class.X_BLOCKS + 1;
	
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_se, ckt->internal_nodelist_se, 
		nd_0, nd_1, count_10, count_20);
	}

	// w block
	if(bx>=1){
		bid_nbr = my_id - 1;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_w, ckt->internal_nodelist_w, 
		nd_0, nd_1, count_10, count_20);	
	}

	// e block
	if((bx<mpi_class.X_BLOCKS-1)){
		bid_nbr = my_id + 1;
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_e, ckt->internal_nodelist_e, 
		nd_0, nd_1, count_10, count_20);
	}

	// nw block
	if((by<mpi_class.Y_BLOCKS-1 && bx>=1)){
		bid_nbr = my_id + mpi_class.X_BLOCKS -1;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_nw, ckt->internal_nodelist_nw, 
		nd_0, nd_1, count_10, count_20);
	}

	// n block
	if((by<mpi_class.Y_BLOCKS-1)){
		bid_nbr = my_id + mpi_class.X_BLOCKS;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_n, ckt->internal_nodelist_n, 
		nd_0, nd_1, count_10, count_20);
	}

	// ne block
	if((by<mpi_class.Y_BLOCKS-1 && bx<mpi_class.X_BLOCKS-1)){
		bid_nbr = my_id + mpi_class.X_BLOCKS + 1;
		
		insert_node_dir(bid_nbr, mpi_class, 
		ckt->bd_nodelist_ne, ckt->internal_nodelist_ne, 
		nd_0, nd_1, count_10, count_20);
	}
}

void Parser::InitialOF(vector<FILE *> & of, int &num_blocks, int &color){
	char  buff[100];
	FILE *f;
	for(int i=0;i<num_blocks;i++){
		sprintf(buff, "./INPUT_FILE/netlist_%d%d.txt", color, i);
		f = fopen(buff, "w");
		of.push_back(f);
	}
}

bool Parser::Is_Top_Layer_Net(Node &p, Node &q){
	bool flag = false;
	for(int i=0;i<(*p_ckts).size();i++){
		Circuit *ckt = (*p_ckts)[i];
		int layers = ckt->layers.size();
		//for(int j=0;j<layers;j++)
			//clog<<"layers: "<<ckt->layers[j]<<endl;
		// both are top layer node
		if (p.pt.z == ckt->layers[layers-1] &&
		    q.pt.z == ckt->layers[layers-1]){
			flag = true;
			break;
		}
	}
	return flag;
}

void Parser::insert_node_list(Node *nd_0, Node *nd_1, int &count_10, 
		int &count_20, int &count_1, int &count_2, NodePtrVector &list, bool &flag){
	// add into bd boundary list
	if(count_10 ==1 && count_20 ==0 && count_2 ==1){
		list.push_back(nd_1);
		if(flag == true)
			nd_1->internal_bd = 1;
	}
	else if(count_10 ==0 && count_20 ==1 && count_1 ==1){
		list.push_back(nd_0);
		if(flag == true)
			nd_0->internal_bd = 1;
	}
}

void Parser::find_bound_line(int &bid_nbr, MPI_CLASS &mpi_class, float &lx, float &ly,  float &ux, float &uy){
	lx = mpi_class.geo[4*bid_nbr];
	ux = mpi_class.geo[4*bid_nbr+2];
	ly = mpi_class.geo[4*bid_nbr+1];
	uy = mpi_class.geo[4*bid_nbr+3];
}

void Parser::insert_node_dir(int &bid_nbr, MPI_CLASS &mpi_class, NodePtrVector &bd_list, NodePtrVector &inter_list, Node *nd_0, Node *nd_1, int &count_10, int &count_20){
	float lx, ly, ux, uy;
	int count_1, count_2;
	bool flag = true;
	
	find_bound_line(bid_nbr, mpi_class, lx, ly,
			ux, uy);

	count_1 = cpr_nd_block(nd_0, lx, ly, ux, uy);
	count_2 = cpr_nd_block(nd_1, lx, ly, ux, uy);

	insert_node_list(nd_0, nd_1, count_1, 
			count_2, count_10, count_20, 
			inter_list, flag);
	flag = false;

	insert_node_list(nd_0, nd_1, count_10, 
			count_20, count_1, count_2, 
			bd_list, flag);
}

//every core parses in data simultaneously
void Parser::block_parse_dots(char *line, Tran &tran, int &my_id){
	char *chs;
	char *saveptr;
	char sname[MAX_BUF];
	Node_TR_PRINT item;
	// clear tran.nodes, especially for core 0
	const char *sep = "= v() \n";
	switch(line[1]){
		case 't': // transient steps
			sscanf(line, "%s %lf %lf", sname, 
				&tran.step_t, &tran.tot_t);
			tran.isTran = 1; // do transient ana;
			//clog<<"step: "<<tran.step_t<<" tot: "<<tran.tot_t<<endl;
			break;
		case 'w': // output length
			chs = strtok_r(line, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			tran.length = atoi(chs);
			//clog<<"out len: "<<tran.length<<endl;
			break;
		case 'p': // print
			chs = strtok_r(line, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			while(chs != NULL){
				chs = strtok_r(NULL, sep, &saveptr);
				if(chs == NULL) break;
				item.name = chs;
				tran.nodes.push_back(item);
				// disribute nodes into cores
			};
			break;
		default: 
			break;
	}
}

void Parser::write_print(Tran &tran, vector<FILE *> &of, MPI_CLASS &mpi_class, char *line){
// then write tran nodes into files
	char *chs;
	char *saveptr;
	char *saveptr_1;
	const char *sep = " _";
	const char *sep_1 = " v()V";
	string ndname;
	vector<string> name_vec;

	chs = strtok_r(line, sep_1, &saveptr_1);

	int num_blocks  = mpi_class.X_BLOCKS * mpi_class.Y_BLOCKS;
	// first print .tran to each file
	for(int j=0;j<num_blocks;j++)
		fprintf(of[j], "%s ", chs);
	chs = strtok_r(NULL, sep_1, &saveptr_1);
	for(int j=0;j<num_blocks;j++)
		fprintf(of[j], "%s ", chs);

	chs = strtok_r(NULL, sep_1, &saveptr_1);
	while(chs != NULL){
		// clog<<"chs: "<<chs<<endl;
		if(string(chs) != "\n")
			name_vec.push_back(string(chs));
		chs = strtok_r(NULL, sep_1, &saveptr_1);
	}
	// clog<<"name_vec, size(): "<<name_vec.size()<<endl;
	// then print the nodes
	for(size_t i=0;i<name_vec.size();i++){
		string tr_nd_name = name_vec[i];
		//clog<<"name: "<<tr_nd_name<<endl;
		char *p = &tr_nd_name[0];
		chs = strtok_r(p, sep, &saveptr);
		ndname = string(chs);
		// skip the X or Y part
		if(ndname[0]=='X' || ndname[0]=='Y') 
			chs = strtok_r(NULL, sep, &saveptr);
		// skip the nz
		chs = strtok_r(NULL, sep, &saveptr);
		// extract x coordinate
		double x = atoi(chs);
		chs = strtok_r(NULL, sep, &saveptr);
		// extract y coordinate
		double y = atoi(chs);

		// judge which blocks this (x,y) belongs to
		for(int j=0;j<num_blocks;j++){
			if(x >= mpi_class.geo[4*j] && x<= mpi_class.geo[4*j+2]){
				if(y >= mpi_class.geo[4*j+1]&& y <= mpi_class.geo[4*j+3]){
					//clog<<"name again: "<<name_vec[i]<<endl;
					//clog<<"belongs to block: "<<j<<" "<<x<<" "<<mpi_class.geo[4*j]<<" "<<mpi_class.geo[4*j+2]<<" "<<y<<" "<<mpi_class.geo[4*j+1]<<" "<<mpi_class.geo[4*j+3]<<endl;
					fprintf(of[j], "v(%s) ", name_vec[i].c_str());	
				}
			}
		}
	}
	for(int j=0;j<num_blocks;j++){
		fprintf(of[j], "\n");
	}
	name_vec.clear();
}
