#include "mpi_class.h"

MPI_CLASS::MPI_CLASS(){
	X_BLOCKS = 2; // # of blocks along x axis
	Y_BLOCKS = 2; // # of blocks along y axis
	x_max = 0;
	y_max = 0;
	x_min = 0;
	y_min = 0;
	overlap_ratio = 0.5;

	start_task = NULL;
	end_task = NULL;
	tasks_n = NULL;

	geo = NULL;
	block_geo = NULL;

	block_size = 0;
}

MPI_CLASS::~MPI_CLASS(){
	delete [] start_task;
	delete [] end_task;
	delete [] tasks_n;
	delete [] geo;
	delete [] block_geo;
}

// assgign tasks
void MPI_CLASS::MPI_Assign_Task(int & num_procs){
	int num_tasks = X_BLOCKS * Y_BLOCKS;
	size_t base = 0;
	for(int i=0;i<num_procs;i++){
		tasks_n[i] = num_tasks / (num_procs);
		if(num_tasks % (num_procs) != 0) {
			if(i < num_tasks % (num_procs))
				tasks_n[i] += 1;
		}
		start_task[i] = base;
		end_task[i] = base + tasks_n[i];
		base += tasks_n[i];
	}
}
