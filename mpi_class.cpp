#include "mpi_class.h"

MPI_CLASS::MPI_CLASS(){
	NUM_NET_TYPE =3;
	X_BLOCKS = 10; // # of blocks along x axis
	Y_BLOCKS = 8; // # of blocks along y axis
	num_blocks = 0;
	cktlist_size = 0;
	x_max = 0;
	y_max = 0;
	x_min = 0;
	y_min = 0;
	overlap_ratio = 0.2;

	len_per_block_x = 0;
	len_per_block_y = 0;
	len_ovr_x = 0;
	len_ovr_y = 0;

	start_task = NULL;
	end_task = NULL;
	tasks_n = NULL;

	geo = NULL;
	block_geo = NULL;
	block_geo_origin = NULL;

	block_size = 0;

	start_proc = NULL;
	end_proc = NULL;
	procs_n = NULL;
	color = 0;
	new_size = 0;
	ranks = NULL;
	layer_color.clear();
}

MPI_CLASS::~MPI_CLASS(){
	delete [] start_task;
	delete [] end_task;
	delete [] tasks_n;
	delete [] geo;
	delete [] block_geo;
	delete [] block_geo_origin;
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

// set 4 internal boundary line
void MPI_CLASS::set_geo_origin(MPI_CLASS &mpi_class){
	// lx
	block_geo_origin[0] = block_geo[0] + len_ovr_x;
	// ly
	block_geo_origin[1] = block_geo[1] + len_ovr_x;
	// ux
	block_geo_origin[2] = block_geo[2] - len_ovr_y;
	// uy
	block_geo_origin[3] = block_geo[3] - len_ovr_y;
}

void MPI_CLASS::Assign_Task(int &num_tasks, int & num_procs){
	size_t base = 0;
	for(int i=0;i<num_procs;i++){
		procs_n[i] = num_tasks / (num_procs);
		if(num_tasks % (num_procs) != 0) {
			if(i < num_tasks % (num_procs))
				procs_n[i] += 1;
		}
		start_proc[i] = base;
		end_proc[i] = base + procs_n[i];
		base += procs_n[i];
	}
}

// n is cktlist_size
void MPI_CLASS::Assign_color(int &my_id, int &n){
	ranks = new int *[n];
	for(int i=0;i<n;i++){
		int procs_num = procs_n[i];
		ranks[i] = new int[procs_num];
		//if(my_id==0)
			//clog<<"ranks range for "<<i<<" th ckt: "<<endl;
		for(int j=0;j<procs_num;j++){
			ranks[i][j] = start_proc[i]+j;
			if(my_id==ranks[i][j]){
				color = i;
			}
			//if(my_id==0)
				//clog<<ranks[i][j]<<" ";
		}
		//if(my_id==0) clog<<endl;
	}
}

void MPI_CLASS::Assign_color_ckt(int &my_id, int &num_procs){
	// start and end proc stores proc number of each cktlist
	start_proc = new int [cktlist_size];
	end_proc = new int [cktlist_size];
	procs_n = new int [cktlist_size];

	Assign_Task(num_procs, cktlist_size);
	ranks = new int *[cktlist_size];
	
	Assign_color(my_id, cktlist_size);

}

void MPI_CLASS::get_char(){
	char c;
	c=getchar();
}
