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
#include "transient.h"
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
	void parse(int &my_id, char * filename, MPI_CLASS &mpi_class, Tran &tran);

	int get_num_layers() const;

	// functions for block
	void set_block_geometry(float *geo, MPI_CLASS &mpi_class);
	int cpr_nd_block(Node &nd, float *geo, int &bid);
	
	int cpr_nd_block(Node *nd, float *geo, int &bid);

	int cpr_nd_block_x(Node *nd, float *geo, int &bid);

	int cpr_nd_block_y(Node *nd, float *geo, int &bid);

	int cpr_nd_block(Node *nd, float &lx, float &ly, 
		float &ux, float &uy);

	void add_node_inter(Node *nd_0, Node *nd_1, 
		MPI_CLASS &mpi_class, Circuit *ckt, int &bid);
	
	void insert_node_list(Node *nd_0, Node *nd_1, 
		int &count_10, int &count_20, int &count_1, 
		int &count_2, NodePtrVector &list, bool &flag);
	
	void find_bound_line(int &bid_nbr, MPI_CLASS &mpi_class, 
		float &lx, float &ly,  float &ux, float &uy);
	
	void insert_node_dir(int &bid_nbr, MPI_CLASS &mpi_class, 
		NodePtrVector &bd_list, NodePtrVector &inter_list, 
		Node *nd_0, Node *nd_1, int &count_10, int &count_20);
	
	void net_to_block(float *geo, MPI_CLASS &mpi_class, Tran &tran);

	void build_block_geo(int &my_id, MPI_CLASS &mpi_class, Tran &tran);

	void set_vdd_map(map<string, string> &vdd_map);

	void second_parse(int &my_id, MPI_CLASS &mpi_class, Tran &tran);
	void parse_dot(char *line, Tran &tran);
	void block_parse_dots(char *line, Tran &tran);
	
	void InitialOF(vector<FILE *> & of, int &num_blocks, int &color);

	void InitialIF(vector<FILE *> & ifs, int &my_id, int &block_size, int &color);
	
private:
	int create_circuits(vector<CKT_LAYER> &ckt_name_info);		// parse the file and create circuits

	int extract_layer(int &my_id, vector<CKT_LAYER >&ckt_layer_info, MPI_CLASS &mpi_class, Tran &tran);
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
