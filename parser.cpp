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
int Parser::create_circuits(){
	FILE * fp;		// used for popen/pclose
	fp = fopen("layer.txt", "r");
	char name[MAX_BUF]="";
			
	int layer, n_circuit=0;

	string prev_ckt_name("");
	string name_string;
	Circuit * p_last_circuit=NULL;
	// now read filename.info to create circuits (they are SORTED)
	while( fscanf(fp, "%s %d", name, &layer) != EOF ){
		//cout<<name_string<<":"<<layer<<endl;
		// compare with previous circuit name 
		name_string.assign(name);
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
	
	fclose(fp);

	// now we know the correct number of layers
	layer_in_ckt.resize(this->n_layer);
	Circuit::layer_dir.resize(this->n_layer);

	return n_circuit;
}

// parse the file
// Note: the file will be parsed twice
// the first time is to find the layer information
// and the second time is to create nodes
void Parser::parse(int &my_id, char * filename, vector<char> &grid_info){
	if(my_id==0){
		this->filename = filename;
		//count = extract_layer();
	}
	else	this->filename = "temp.txt";

	// for y200, this part cost about 5.6s
	clock_t t1, t2;
	t1 = clock();
	store_in_vector(my_id, grid_info);
	t2 = clock();
	//if(my_id==0) clog<<" send cost: "<<
	//	1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;

	return;
	FILE * f;
	f = fopen(filename, "r");
	if( f == NULL ) 
		report_exit("Input file not exist!\n");
	// first time parse:
	create_circuits();

	// second time parser:
	char line[MAX_BUF];
	string l;
	while( fgets(line, MAX_BUF, f) != NULL ){
		char type = line[0];
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

int Parser::store_in_vector(int &my_id, vector<char> &grid_info){
	FILE *f;
	f = fopen(filename, "r");
	if(f==NULL) report_exit("Input file not exist!\n");
	char line[MAX_BUF];
	
	while(fgets(line, MAX_BUF, f)!=NULL){
		for(int i=0;line[i]!='\n';i++){
			grid_info.push_back(line[i]);
		}
		grid_info.push_back('\n');
	}
	fclose(f);

	long long grid_size = grid_info.size();

	MPI_Bcast(&grid_size, 1, MPI_LONG_LONG, 0, 
		MPI_COMM_WORLD);

	if(my_id!=0) grid_info.resize(grid_size);

	MPI_Bcast(&grid_info[0], grid_size, MPI_CHAR, 0, 
		MPI_COMM_WORLD);

	//clog<<grid_info.size()<<endl;
	//if(my_id==11)
	//for(size_t i=grid_info.size()-100;i<grid_info.size();i++)
		//clog<<grid_info[i];
	return 0;
}

int Parser::extract_layer(){
	FILE *f;
	f = fopen(filename, "r");
	if(f==NULL)	report_exit("Input file not exist!\n");
	char line[MAX_BUF];
	char word[MAX_BUF];
	string word_s;
	char name[10];
	int count=0;

	FILE *fp;
	fp = fopen("layer.txt", "w");
	if(fp == NULL) report_exit("Output file not exist!\n");
	while(fgets(line, MAX_BUF, f)!=NULL){
		count++;
		if(line[0]=='*'){
			// copy the entire line into stringstream
			stringstream ss;
			ss<< line;	
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
						if(vdd_count==1)
							fprintf(fp, "%s ", name);
						vdd_count++;
					}
				}
				// extract layer number
				if(word_count==4){
					fprintf(fp, "%s \n", word);
				}
				word_count++;
			}
		}
	}
	fclose(f);
	fclose(fp);
	return count;
}
