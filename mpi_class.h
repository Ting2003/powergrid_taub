#ifndef __MPI_CLASS_H_
#define __MPI_CLASS_H_

#include <iostream>

using namespace std;

class MPI_CLASS{
public:
	// member function
	MPI_CLASS();
	~MPI_CLASS();
	// member
	int X_BLOCKS;
	int Y_BLOCKS;
	long x_max;
	long y_max;
	long x_min;
	long y_min;
	float overlap_ratio;

	int *start_task;
	int *end_task;
	int *tasks_n;

	float *geo;
	float *block_geo;

	int block_size;

	// function
	void MPI_Assign_Task(int & num_procs);
};

#endif
