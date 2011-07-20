#ifndef __MAIN_H__
#define __MAIN_H__

#include <iostream>
#include <sstream>
#include <unistd.h>
#include "global.h"
#include "util.h"
#include "point.h"
#include "vec.h"
#include "triplet.h"
#include "parser.h"
#include "circuit.h"
//#include "hash_mat.h"
using namespace std;

// assgign procs for circuit
void Assign_Task(int &num_tasks, int & num_procs, int *start_task,int *end_task, int *tasks_n);

void Assign_color(int &my_id, int &color, int **ranks, int &n,
	int *procs_n, int *start_proc);


#endif
