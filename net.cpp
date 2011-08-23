#include <iostream>
#include <iomanip>
#include "node.h"
#include "net.h"
using namespace std;

// a, b are pointers to two nodes
//Net::Net(NET_TYPE t, string n, double v, Node * a, Node * b): 
//type(t), name(n), value(v){
Net::Net(NET_TYPE t, double v, Node * a, Node * b):
	type(t),value(v){
	ab[0]=a;
	ab[1]=b;
	flag_bd = false;
}

ostream & operator << (ostream & os, const Net & net){
	os//<<net.name
		<<"("
	    <<net.ab[0]->name<<","
	    <<net.ab[1]->name<<")="
 	    <<scientific<<net.value;
	return os;
}
