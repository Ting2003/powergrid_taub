#ifndef __MPI_CLASS_H_
#define __MPI_CLASS_H_

#include <iostream>
#include "circuit.h"
#include "global.h"
#include "mpi.h"

using namespace std;

const int MAX_ITERATION = 1000;

class MPI_CLASS{
public:
	// member function
	MPI_CLASS();
	~MPI_CLASS();
	// member
	int *start_task;
	int *end_task;
	int *tasks_n;

	int block_size;

	// function
	void MPI_Assign_Task(int & num_procs);
};

#endif
