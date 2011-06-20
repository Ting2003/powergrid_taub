#ifndef __NET_H__
#define __NET_H__

#include <string>
#include "global.h"
#include "point.h"
using namespace std;

class Node;

class Net{
public:
	//Net(NET_TYPE type, string name, double value, Node * a, Node * b);
	Net(NET_TYPE type, double value, Node * a, Node * b);
	friend ostream & operator << (ostream & os, const Net & net);

	NET_TYPE type;	// identify the type
	//string name;
	double value;
	Node* ab[2];	// two connection node
};

#endif
