// ----------------------------------------------------------------//
// Filename : parser.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of Parser class
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#ifndef __PARSER_H__
#define __PARSER_H__

#include <vector>
#include "global.h"
#include "circuit.h"
#include "mpi_class.h"
#include <algorithm>
#include <map>
#include <string>
using std::vector;

struct CKT_LAYER{
	char name[10];
	int layer;
};
// given an input file, parse it and store corresponding result into Circuit
class Parser{
public:
	// supply Circuit objects
	Parser(vector<Circuit*> * ckts);
	~Parser();

	// parser a input file and construct the circuit
	void parse(int &my_id, char * filename, MPI_CLASS &mpi_class);

	int get_num_layers() const;

	// functions for block
	void set_block_geometry(float *geo, MPI_CLASS &mpi_class);
	int cpr_nd_block(Node &nd, float *geo, int &bid);
	
	int cpr_nd_block(Node *nd, float *geo, int &bid);

	void net_to_block(map<string, string> &vdd_map, float *geo, MPI_CLASS &mpi_class);

	void build_block_geo(int &my_id, MPI_CLASS &mpi_class);

	void set_vdd_map(map<string, string> &vdd_map);

	void second_parse(int &my_id, MPI_CLASS &mpi_class);
	
	void InitialOF(vector<FILE *> & of, int &num_blocks, int &color);

	void InitialIF(vector<FILE *> & ifs, int &my_id, int &block_size, int &color);
	
private:
	int create_circuits(vector<CKT_LAYER> &ckt_name_info);		// parse the file and create circuits

	int extract_layer(int &my_id, vector<CKT_LAYER >&ckt_layer_info, MPI_CLASS &mpi_class);
	bool sort(vector <CKT_LAYER> &a);
	
	bool Is_Top_Layer_Net(Node &p, Node &q);
	void try_change_via(Net *);

	//void insert_net_node(string line);
	void insert_net_node(char * line, int &count, MPI_CLASS &mpi_class);
	void extract_node(char * str, Node & nd);
	void update_node(Net * net);

	char * filename;		  // input file name
	int n_layer;			  // total number of layers
	vector<Circuit*> * p_ckts;	  // pointer to circkt list
	vector<int> layer_in_ckt;	  // which circuit a layer belong
};

// Trick: try to modify the net
inline void Parser::try_change_via(Net * net){
	// is it a via?
	if( net->ab[0]->get_layer() == net->ab[1]->get_layer() )
		return;

	// make it a zero voltage via
	if( net->type == RESISTOR && net->value < 1e-4 ){
		net->type = VOLTAGE;
		net->value = 0.0;
	}
}

#endif
