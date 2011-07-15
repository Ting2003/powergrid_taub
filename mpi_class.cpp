#include "mpi_class.h"

MPI_CLASS::MPI_CLASS():
	ckt(NULL),
	start_task(NULL),
	end_task(NULL),
	tasks_n(NULL),
	b_new_info(NULL),
	L_h(NULL),
	L_nz_d(NULL),
	base_nz_d(NULL),
	L_n_d(NULL),
	base_n_d(NULL),
	send_nz(NULL),
	send_n(NULL),
	L_d(NULL),
	b_x_d(NULL),
	L_nz_dd(NULL),
	L_n_dd(NULL){
	total_n=0;
	total_nz=0;
}

MPI_CLASS::~MPI_CLASS(){
	delete [] start_task;
	delete [] end_task;
	delete [] tasks_n;
	delete [] b_new_info;
	delete [] L_h;
	delete [] L_nz_d;
	delete [] base_nz_d;
	delete [] L_n_d;
	delete [] base_n_d;
	delete [] send_nz;
	delete [] send_n;
	delete [] L_d;
	delete [] b_x_d;
	delete [] L_nz_dd;
	delete [] L_n_dd;
}

void MPI_CLASS::solve_mpi_IT(int &my_id, int &num_procs){
	// shoule be included in solve_mpi_IT	
	int num_tasks = ckt->block_info.size();

	start_task = new int [num_procs];
	end_task = new int [num_procs];
	tasks_n = new int [num_procs];	
	send_nz = new int [num_procs];
	send_n = new int [num_procs];
	base_nz_d = new int [num_procs];
	base_n_d = new int [num_procs];

	for(int i=0;i<num_procs;i++){
		start_task[i]=0;
		end_task[i]=0;
		tasks_n[i]=0;
		send_nz[i]=0;
		send_n[i]=0;
		base_nz_d[i]=0;
		base_n_d[i]=0;
	}
	
	if(my_id==0){
		MPI_Assign_Task(num_tasks, num_procs);
		// calculate index arrays
		block_mpi_setup();
	}
	
	// first, scatter send_nz into processors to initialize
	num_blocks = ckt->block_info.size();

	// bcast total block numbers to others
	MPI_Bcast(&num_blocks, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// scatter total block numbers for each processor
	MPI_Scatter(tasks_n, 1, MPI_INT, &block_size, 
			1, MPI_INT, 0, MPI_COMM_WORLD);
	
	// scatter total nz for each processor
	MPI_Scatter(send_nz, 1,  MPI_INT, &L_nz, 1, 
			MPI_INT, 0, MPI_COMM_WORLD);
		
	// scatter total n for each processor
	MPI_Scatter(send_n, 1, MPI_INT, &L_n, 1, 
			MPI_INT, 0, MPI_COMM_WORLD);

	L_d = new float [L_nz];
	b_x_d = new float [L_n];

	// scatter L into corresponding processors
	MPI_Scatterv(L_h, send_nz, 
		base_nz_d, MPI_FLOAT, 
		L_d, L_nz, MPI_FLOAT, 0, 
		MPI_COMM_WORLD);

	L_nz_dd = new int [block_size];
	L_n_dd = new int [block_size];
	
	// scatter nz of L of each block within each processor
	MPI_Scatterv(L_nz_d, tasks_n, 
		start_task, MPI_INT, 
		L_nz_dd, block_size, MPI_INT, 
		0, MPI_COMM_WORLD);
	// scatter n of rhs of each block within each processor
	MPI_Scatterv(L_n_d, tasks_n, 
		start_task, MPI_INT, 
		L_n_dd, block_size, MPI_INT, 
		0, MPI_COMM_WORLD);
}

void MPI_CLASS::solve_mpi_iteration(int &my_id, int&num_procs){
	// 0 rank cpu scatter rhs to all other processors
	MPI_Scatterv(b_new_info, send_n, base_n_d, 
		MPI_FLOAT, b_x_d, L_n, MPI_FLOAT, 0, 
		MPI_COMM_WORLD);
	
	//MPI_Barrier(MPI_COMM_WORLD);

	// there are in total 'block_size' blocks 
	// in each processor
	int base_nz = 0;
	int base_n = 0;
	for(int i=0;i<block_size;i++){
		Algebra::solve_CK_for_back_sub(my_id,
				L_nz_dd[i], L_d, L_n_dd[i], 
				b_x_d, base_nz, base_n);
		/*if(my_id==0 && i ==0){
		clock_t t1, t2;
		t1 = clock();
		for(int j=0;j<1000;j++)
			Algebra::solve_CK_for_back_sub(my_id,
				L_nz_dd[i], L_d, L_n_dd[i], 
				b_x_d, base_nz, base_n);

		t2 = clock();
		clog<<"time for bid_0 is: (1000 iter: ): "<<1.0*(t2-t1) / CLOCKS_PER_SEC<<endl;
		}*/

		base_nz += L_nz_dd[i];
		base_n += L_n_dd[i];
	}
	//MPI_Barrier(MPI_COMM_WORLD);

	// 0 rank cpu will gather all the solution from b_x_d
	// to b_new_info
	MPI_Gatherv(b_x_d, L_n, MPI_FLOAT, b_new_info, send_n,
		base_n_d, MPI_FLOAT, 0, MPI_COMM_WORLD);
}

double MPI_CLASS::solve_iteration(int &my_id, int&num_procs){
	int iter=0;
	double diff = 0.0;
	bool successful = false;

	clock_t t1, t2;
	t1= clock();
	while(iter<MAX_ITERATION){
		if(my_id==0){
			solve_iteration_1();
		}

		solve_mpi_iteration(my_id, num_procs);

		if(my_id==0){
			diff = solve_iteration_2();
		}
		iter++;
		MPI_Bcast(&diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//if(my_id==0)
			//clog<<"iter, diff: "<<iter<<" "<<diff<<endl;
		if(diff < Circuit::EPSILON){
			successful=true;
			break;
		}
	}

	t2 = clock();
	if(my_id==0) clog<<"iter time is: "<<1.0*(t2-t1) / CLOCKS_PER_SEC<<endl;
	if(my_id==0){
		clog<<"# iter: "<<iter<<endl;
		ckt->get_voltages_from_block_LU_sol();
		ckt->get_vol_mergelist();
		for(size_t i=0;i<ckt->block_info.size();i++){
			if(ckt->block_info[i].count > 0)
				ckt->block_info[i].free_block_cholmod(ckt->cm);
		}
		cholmod_finish(ckt->cm);
	}
	return successful;

	return diff;
}

// assgign tasks
void MPI_CLASS::MPI_Assign_Task(int &num_tasks,int & num_procs){
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

void MPI_CLASS::block_mpi_setup(){
	for(size_t i=0;i<ckt->block_info.size();i++){	
		total_n += ckt->block_info[i].count;
		total_nz += ckt->block_info[i].L_h_nz;
	}

	b_new_info = new float [total_n];
	for(size_t i=0;i<total_n;i++){
		b_new_info[i]=0;
	}

	// total_nz is the number of the triplet, 
	// length of L_h should be 3*total_nz;
	L_h = new float [3*total_nz];
	size_t base =0;
	for(size_t i=0;i<ckt->block_info.size();i++){
		for(size_t j=0;j<3*ckt->block_info[i].L_h_nz;j++){
			L_h[base+j] = ckt->block_info[i].L_h[j];
		}
		base += 3*ckt->block_info[i].L_h_nz;
	}
	L_nz_d = new int [ckt->block_info.size()];
	L_n_d = new int [ckt->block_info.size()];

	size_t base_nz = 0;
	size_t base_n = 0;
	size_t j=0;
	size_t k=0;
	for(size_t i=0;i<ckt->block_info.size();i++){
		L_nz_d[i] = 3*ckt->block_info[i].L_h_nz;
		L_n_d[i] = ckt->block_info[i].count;
		if((int)i == start_task[j] && start_task[j]<end_task[j]){
			base_nz_d[j] = base_nz;
			base_n_d[j] = base_n;
			j++;
		}
		if((int) i == end_task[k] && start_task[k]<end_task[k]){
			send_nz[k] = base_nz- base_nz_d[k];
			send_n[k] = base_n - base_n_d[k];
			k++;
		}
		base_nz += 3*ckt->block_info[i].L_h_nz;
		base_n += ckt->block_info[i].count;
	}
	// define the last send_nz and send_n
	send_nz[k] = 3*total_nz - base_nz_d[k];
	send_n[k] = total_n - base_n_d[k];
}

// solve blocks with mpi: multi-core
// One iteration during solving the circuit, for any block B:
// 1. update the righthand-side of the matrix of B
// 2. solve the matrix
// 3. update node voltages
// 4. track the maximum error of solution
void MPI_CLASS::solve_iteration_1(){	
	// 0 rank cpu update all the rhs
	size_t base = 0;
	for(size_t i=0;i<ckt->block_info.size();i++){
		Block &block = ckt->block_info[i];
		if(block.count ==0) continue;
		block.update_rhs();
		// backup the old voltage value	
		for(size_t j=0; j<block.count;j++){
			b_new_info[base+j] = block.bnewp[j];
			block.x_old[j] = block.x_new[j];
		}
		base += block.count;
	}
}

double MPI_CLASS::solve_iteration_2(){
	float diff = .0, max_diff = .0;
	// modify node voltage with OMEGA and old voltage value
	size_t base = 0;
	for(size_t i=0;i<ckt->block_info.size();i++){
		Block &block = ckt->block_info[i];
		// copy block.x_new from large b_new_info
		for(size_t j=0;j<block.count;j++){
			block.x_new[j] = b_new_info[base+j];
		}
		base += block.count;
		diff = ckt->modify_voltage(block, block.x_old);
		if( max_diff < diff ) max_diff = diff;
	}

	return max_diff;
}
