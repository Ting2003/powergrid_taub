#ifndef __NODE_H__
#define __NODE_H__

#include <string>
#include <algorithm>
#include "global.h"
#include "point.h"
#include "net.h"
#include <vector>
#include <iostream>
//#include "block.h"
using namespace std;

class Net;
class Circuit;
// a Node in the network
class Node{
public:
	// member functions
	Node();
	Node(string name, Point _pt, bool flag=false, double v=0.0);
	Node(const Node & nd);
	Node & operator = (const Node & nd);
	void set_nbr(DIRECTION dir, Net * name);
	Net * get_nbr_net(DIRECTION dir) const;

	// Ting: get_nbr_node
	Node * get_nbr_node(Node *node, DIRECTION dir) const;
	Node get_nbr_node(Node node, DIRECTION dir) const;
	Node * get_nbr_node(DIRECTION dir) const;

	int get_layer() const;

	bool isX() const;
	bool is_ground() const;

	double get_value() const;
	void set_value(double v);

	vector<size_t> & get_block_id() const;

	bool inside_block(size_t block_id) const;

	bool is_mergeable() const;

	friend ostream & operator << (ostream & os, const Node & node);
	friend class Circuit;
	friend class Block;
	friend class Parser;

	////////////////////////////////////////////////////////////
	// member variables
	string name;		// node name
	Point pt;		// coordinate
	// only 2 possible cases:
	// {TOP, BOTTOM, EAST, WEST}
	// {TOP, BOTTOM, NORTH, SOUTH}
	Net * nbr[6];		// neighboring nets

	size_t rid;		// id in rep_list

private:
	double value;		// voltage
	bool flag;		// mark if the node is an X
	Node * rep;		// representative, if somewhere is short-circuit
	vector<size_t> blocklist;	// belongs to which block
	vector<size_t> id_in_block;	// local index inside block	

	Node * end[4];		// south - north (west-east) ends
	double eqvr[4];		// equivalent resisotrs
};      	

inline bool Node::isX() const{return flag;}

//inline bool Node::is_ground() const{return name == "0";}
// use a tricky way to speed up
inline bool Node::is_ground() const{return pt.x<0;}

inline int Node::get_layer() const{ return pt.z; }

inline double Node::get_value() const{return value;}

inline void Node::set_value(double v){value = v;}

inline void Node::set_nbr(DIRECTION dir, Net * net){ nbr[dir] = net; }

inline vector<size_t> & Node::get_block_id() const{return rep->blocklist;}

inline bool Node::inside_block(size_t block_id) const{
	vector<size_t>::const_iterator it;
	it = find(blocklist.begin(), blocklist.end(), block_id);
	if( it == blocklist.end() ) return false;
	return true;
}

inline Node * Node::get_nbr_node(DIRECTION dir) const{
	if( nbr[dir] == NULL ) return NULL;
	Node * nbr_node = nbr[dir]->ab[0];
	return (nbr_node != this? nbr_node: nbr[dir]->ab[1]);
}

inline Net * Node::get_nbr_net(DIRECTION dir) const{
	return nbr[dir];
}

inline bool Node::is_mergeable() const{
	bool flag;
	flag = nbr[TOP] == NULL && nbr[BOTTOM] == NULL &&
	     ((nbr[EAST]  != NULL && nbr[WEST]  != NULL &&
	       nbr[NORTH] == NULL && nbr[SOUTH] == NULL)
	    ||(nbr[NORTH] != NULL && nbr[SOUTH] != NULL &&
	       nbr[EAST]  == NULL && nbr[WEST] == NULL));

	return flag;
}

#endif
