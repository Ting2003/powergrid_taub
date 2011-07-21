#include "main.h"
#include "cholmod.h"
#include "mpi.h"

const char * usage="Usage: %s [-eorbILifl] benchmark\n\
    -e EPSILON\n\
    -o OMEGA\n\
    -r overlap ratio\n\
    -b max block nodes\n\
    -I block iterative (default)\n\
    -L direct LU\n\
    -i input file\n\
    -f output file\n\
    -l log file (default to screen)\n"
;

const char * usage2="Usage: %s -i input -f output\n";

int main(int argc, char * argv[]){
	int my_id;
	int num_procs;

	//double mpi_t1, mpi_t2;
	//mpi_t1 = MPI_Wtime();
	
	int c;
	int mode=0;
	double epsilon, omega, overlap_ratio;
	size_t max_block_nodes;
	//char * logfile="/dev/null";
	char * logfile=NULL;
	char * input=NULL, * output=NULL;
	bool input_flag = false, output_flag = false;
	Circuit::get_parameters(epsilon, omega, overlap_ratio, 
			max_block_nodes, mode);

	while( ( c = getopt(argc, argv, "i:f:e:o:r:b:l:LI")) != -1 ){
		switch(c){
		case 'e':
			epsilon = atof(optarg);
			break;
		case 'o':
			omega = atof(optarg);
			break;
		case 'r':
			overlap_ratio = atof(optarg);
			break;
		case 'b':
			max_block_nodes = atof(optarg);
			break;
		case 'L':
			mode = 1;
			break;
		case 'I':
			mode = 0;
			break;
		case 'l':
			logfile = optarg;
			break;
		case 'i':
			input = optarg;
			input_flag = true;
			break;
		case 'f':
			output = optarg;
			output_flag = true;
			break;
		case '?':
		default:
			fprintf(stderr, usage2, argv[0]);
			exit(EXIT_FAILURE);
		}
	}
	//if( argc == optind ) report_exit(usage2);
	if( !input_flag || ! output_flag ){
		fprintf(stderr, usage2, argv[0]);
		exit(EXIT_FAILURE);
	}
	open_logfile(logfile);

	// mpi variable for grouping
	MPI_Comm comm, new_comm;
	MPI_Group orig_group, new_group;
	int new_rank, rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
	//if(my_id==0) clog<<"num_procs: "<<num_procs<<endl;	

	if( freopen(output, "w", stdout) == NULL )
		report_exit("Ouptut file error\n");

	Circuit::set_parameters(epsilon, omega, overlap_ratio, 
			max_block_nodes, mode);
	// start to parfile
	vector<Circuit *> cktlist;
	Parser parser(&cktlist);
	// cktlist size are found now
	parser.parse_1(my_id, input);

	int cktlist_size = cktlist.size();
	// start and end proc stores proc number of each cktlist
	int *start_proc, *end_proc, *procs_n;
	// color as a index for grouping processors
	int color, new_size;
	start_proc = new int [cktlist_size];
	end_proc = new int [cktlist_size];
	procs_n = new int [cktlist_size];

	Assign_Task(num_procs, cktlist_size, start_proc,
		end_proc, procs_n);
	
	int **ranks;
	ranks = new int *[cktlist_size];
	
	Assign_color(my_id, color, ranks, cktlist_size, 
		procs_n, start_proc);

	MPI_Comm_split(MPI_COMM_WORLD, color, my_id, 
		&new_comm);
	MPI_Comm_group(new_comm, &new_group);
	
	MPI_Group_size(new_group, &new_size);
	MPI_Group_rank(new_group, &new_rank);
	//if(my_id==0) clog<<"new size is: "<<new_size<<endl;
	//clog<<"old and new rank: "<<my_id<<" "<<new_rank<<endl;		
	
	clock_t t1,t2;
	t1=clock();
	parser.parse_2(new_rank, color);
	MPI_Barrier(MPI_COMM_WORLD);
	// after parsing, this mem can be released
	t2=clock();
	//if(new_rank==0) clog<<"Parse time="<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;

	double mpi_t11, mpi_t12;
	mpi_t11 = MPI_Wtime();
	
	// now each group of processor only deals with its color ckt
	int i = color;
	Circuit * ckt = cktlist[i];	

	if(new_rank ==0)
		clog<<"Solving "<<ckt->get_name()<<endl;
	ckt->solve(new_rank, new_size, new_comm);

	//Write_send(my_id, color, new_rank, cktlist_size, ckt);

	if(new_rank ==0)
		ckt->print(color);

	free(ckt);
	
	mpi_t12 = MPI_Wtime();
	
	// output a single ground node
	if(new_rank==0)
		clog<<"solve using: "<<1.0*(mpi_t12-mpi_t11)<<endl;

	if(my_id==0){
		close_logfile();
	}

	//MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Group_free(&new_group);
	//MPI_Group_rank(orig_group, &rank);
	//clog<< "rank after union: "<<rank<<endl;

	MPI_Group_free(&orig_group);

	MPI_Finalize();
	//close_logfile();
	//cout<<"close logfile. "<<endl;
	return 0;
}

// assgign procs for circuit
void Assign_Task(int &num_tasks, int & num_procs, int *start_task,int *end_task, int *tasks_n){
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

// n is cktlist_size
void Assign_color(int &my_id, int &color, int **ranks, int &n,
	int *procs_n, int *start_proc){
	ranks = new int *[n];
	for(int i=0;i<n;i++){
		int procs_num = procs_n[i];
		ranks[i] = new int[procs_num];
		if(my_id==0)
			clog<<"ranks range for "<<i<<" th ckt: "<<endl;
		for(int j=0;j<procs_num;j++){
			ranks[i][j] = start_proc[i]+j;
			if(my_id==ranks[i][j])
				color = i;
			if(my_id==0)
				clog<<ranks[i][j]<<" ";
		}
		if(my_id==0) clog<<endl;
	}
}

// serially output into logfile
void Write_send(int &my_id, int &color, int &new_rank, int &cktlist_size, Circuit *ckt){
	int wf = 0;
	int send_id=0;
	for(int i=0;i<cktlist_size;i++){
		if(color == wf && new_rank ==0){
			send_id = my_id;
			clog<<"print id is: "<<my_id<<endl;
			//ckt->print();
		}
		wf ++;
		MPI_Bcast(&wf, 1, MPI_INT, my_id, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if(my_id==0)
		printf("G %.5e\n", 0.0);
}
