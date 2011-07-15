#ifndef __MPI_CLASS_H_
#define __MPI_CLASS_H_

#include <iostream>
#include "circuit.h"
#include "mpi.h"

using namespace std;

const int MAX_ITERATION = 1000;

class MPI_CLASS{
public:
	// member function
	MPI_CLASS();
	~MPI_CLASS();
	void solve_mpi_IT(int &my_id, int &num_procs);
	void MPI_Assign_Task(int &num_tasks, int &num_procs);
	void block_mpi_setup();
	// updates nodes value in each iteration
	void solve_iteration_1();

	double solve_iteration_2();
	void solve_mpi_iteration(int &my_id, int&num_procs);

	double solve_iteration(int &my_id, int&num_procs);
	
	// member
	Circuit *ckt;
	// start_task and end_task stores begin
	// and end block number for each processor
	int *start_task;
	int *end_task;
	// tasks_n stores # of blocks for each processor
	int *tasks_n;
	// b_new_info is the global solution array
	float *b_new_info;
	// L_h is the global factorized array
	float *L_h;
	// # of nz in L of each processor
	int *L_nz_d;
	// base for L_h for each processor
	int *base_nz_d;
	// # of n in each processor
	int *L_n_d;
	// base for rhs for each processor
	int *base_n_d;
	// # of nz in each processor
	int *send_nz;
	// # of n in each processor
	int *send_n;
	// L for each processor
	float *L_d;
	// solution array for each processor
	float *b_x_d;
	// # of nz in each block within a processor
	int *L_nz_dd;
	// # of n in each block within a processor
	int *L_n_dd;
	
	// L_nz and L_n are # of nz in L and # of n in rhs
	// L_nz here is the non-zero of triplet including coord
	int L_nz; int L_n;
	// block_size is the number of tasks for each cpu
	int block_size;
	
	// first, scatter send_nz into processors to initialize
	int num_blocks;

	size_t total_n, total_nz;
};

#endif
