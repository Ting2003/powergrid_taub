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
	
	long z, y, x;
	bool flag = false;
	char * chs;
	char * saveptr;
	nd.name.assign(str);

	char * l = str;
	const char * sep = "_n";
	chs = strtok_r(l, sep, &saveptr); // initialize
	if( chs[0] == 'X' ){
		flag = true;
		chs = strtok_r(NULL, sep, &saveptr);
	}
	z = atol(chs);
	chs = strtok_r(NULL, sep, &saveptr);
	x = atol(chs);
	chs = strtok_r(NULL, sep, &saveptr);
	y = atol(chs);

	nd.pt.set(x,y,z);
	nd.flag = flag;
}

// given a line, extract net and node information
void Parser::insert_net_node(char * line){
	static char sname[MAX_BUF];
	static char sa[MAX_BUF];
	static char sb[MAX_BUF];
	static Node nd[2];
	Node * nd_ptr[2];	// this will be set to the two nodes found
	double value;
	sscanf(line, "%s %s %s %lf", sname, sa, sb, &value);

	if( sa[0] == '0' ) { nd[0].pt.set(-1,-1,-1); }
	else extract_node(sa, nd[0]);

	if( sb[0] == '0' ) { nd[1].pt.set(-1,-1,-1); }
	else extract_node(sb, nd[1]);

	// insert these two node into Circuit according to the node's layer types
	// Note: 1. these two nodes may exist already, need to check
	//       2. these two nodes must be in the same circuit (network), 
	//       (except 0), so their layer_type must be the same

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
			ckt->add_node(nd_ptr[i]);
			if( nd_ptr[i]->isX() )	     // determine circuit type
				ckt->set_type(WB);

			// find the coordinate max and min
			size_t x = nd[i].pt.x;
			size_t y = nd[i].pt.y;
			if( x < ckt->x_min ) ckt->x_min = x;
			if( y < ckt->y_min ) ckt->y_min = y;
			if( x > ckt->x_max ) ckt->x_max = x;
			if( y > ckt->y_max ) ckt->y_max = y;
		}
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
	default:
		report_exit("Invalid net type!\n");
		break;
	}

	// create a Net
	Net * net = new Net(net_type, value, nd_ptr[0], nd_ptr[1]);

	// trick: when the value of a resistor via is below a threshold,
	// treat it as a 0-voltage via
	if( Circuit::MODE == (int)IT ) {
		try_change_via(net);
	}

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

	if( a->get_layer() == b->get_layer() && a->pt != b->pt ){
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
		a->flag = true;		// set a to be X node
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
void Parser::parse(int &my_id, char * filename, MPI_CLASS &mpi_class){	
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
		extract_layer(my_id, ckt_name_info);
	
	int ckt_name_info_size = ckt_name_info.size();
	
	MPI_Bcast(&ckt_name_info_size, 1, MPI_INT, 0,
			MPI_COMM_WORLD);
	if(my_id!=0) ckt_name_info.resize(ckt_name_info_size);

	MPI_Bcast(&ckt_name_info[0], ckt_name_info_size, 
			MPI_Vector, 0, MPI_COMM_WORLD);

	// first time parse:
	create_circuits(ckt_name_info);

	build_block_geo(my_id);

	second_parse(my_id, mpi_class);
}

void Parser::build_block_geo(int &my_id){
	float * geo;
	float * block_geo;
	// block info in cpu 0
	geo = new float[X_BLOCKS *Y_BLOCKS *4];
	// block info in other cpus
	block_geo = new float [4];

	// update block geometry
	if(my_id==0) set_block_geometry(geo);
	MPI_Scatter(geo, 4, MPI_FLOAT, block_geo, 4, MPI_FLOAT, 
			0, MPI_COMM_WORLD);
	
	if(my_id==0){
		net_to_block(geo);
	}
}

void Parser::second_parse(int &my_id, MPI_CLASS &mpi_class){
	char  buff[100];
	FILE *f = NULL;

	int color = 0;
	int block_size = mpi_class.block_size;
	if(block_size==0) return;
	for(int i=0;i<block_size;i++){
		sprintf(buff, "netlist_%d%d.txt", color, my_id);
		f = fopen(buff, "r");
		if(f==NULL) report_exit("Input file not exist!\n");
	}
	
	char line[MAX_BUF];
	char type;

	while( fgets(line, MAX_BUF,f)!=NULL){
		type = line[0];
		switch(type){
			case 'r': // resistor
			case 'R':
			case 'v': // VDD
			case 'V':
			case 'i': // current
			case 'I':
				insert_net_node(line);
				break;
			case '.': // command
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

int Parser::extract_layer(int &my_id, 
		vector<CKT_LAYER > &ckt_layer_info){
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
		}
	}
	fclose(f);
	// sort resulting vector by the ckt name
	sort(ckt_layer_info);
	return 0;
}

bool Parser::sort(vector<CKT_LAYER > &a){
	CKT_LAYER tmp;
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

void Parser::set_block_geometry(float *geo){
	double x, y;
	x = (double)(x_max-x_min+0.5) /X_BLOCKS;
	y = (double)(y_max-y_min+0.5) /Y_BLOCKS;
	//if( fzero(x) ) x = 1.0;
	//if( fzero(y) ) y = 1.0;
	double len_per_block_x = x;
	double len_per_block_y = y;
	double len_ovr_x = x * overlap_ratio;
	double len_ovr_y = y * overlap_ratio;
	clog<<"len_x, len_y: "<<x<<" / "<<y<<endl;
	clog<<"len_ovr_x, len_ovr_y: "<<len_ovr_x
		<<" / "<<len_ovr_y<<endl;

	size_t bid = 0;
	// update block 4 corners
	for(size_t y=0;y<Y_BLOCKS;y++){
		for(size_t x=0;x<X_BLOCKS;x++){
			bid = y * X_BLOCKS + x;
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

void Parser::net_to_block(float *geo){
	static char sname[MAX_BUF];
	static char sa[MAX_BUF];
	static char sb[MAX_BUF];
	static Node nd[2];
	double value;

	FILE *f;
	f = fopen(this->filename, "r");
	if(f==NULL) report_exit("Input file not exist!\n");

	char line[MAX_BUF];

	int color = 0;
	vector<FILE *> of;
	int num_blocks  = X_BLOCKS * Y_BLOCKS;
	InitialOF(of, num_blocks, color);

	int count_1 = 0, count_2 = 0;
	while( fgets(line, MAX_BUF,f)!=NULL){
		if(line[0]=='r' || line[0] =='R' ||
		   line[0]=='v' || line[0]=='V' ||
		   line[0]=='i' || line[0]=='I'){
			sscanf(line, "%s %s %s %lf", 
					sname, sa, sb, &value);
			if( sa[0] == '0' ){ nd[0].pt.set(-1,-1,-1); }
			else extract_node(sa, nd[0]);
			if( sb[0] == '0' ){ nd[1].pt.set(-1,-1,-1); }
			else	extract_node(sb, nd[1]);
			for(int i=0;i<X_BLOCKS * Y_BLOCKS;i++){
				count_1 = cpr_nd_block(nd[0], geo, i);
				count_2 = cpr_nd_block(nd[1], geo, i);

				if(count_1 + count_2 >=1){
					fprintf(of[i], "%s", line);
				}
			}
		}
	}
	fclose(f);
	for(int i=0;i<num_blocks;i++)
		fclose(of[i]);
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

void Parser::InitialOF(vector<FILE *> & of, int &num_blocks, int &color){
	char  buff[100];
	FILE *f;
	for(int i=0;i<num_blocks;i++){
		sprintf(buff, "netlist_%d%d.txt", color, i);
		f = fopen(buff, "w");
		of.push_back(f);
	}
}

void Parser::InitialIF(vector<FILE *> & ifs, int &my_id, int &block_size, int &color){
	char  buff[100];
	FILE *f;

	if(block_size==0) return;
	for(int i=0;i<block_size;i++){
		sprintf(buff, "netlist_%d%d.txt", color, my_id);
		f = fopen(buff, "r");
		ifs.push_back(f);
	}
}
