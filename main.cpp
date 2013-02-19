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
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	if(my_id==0) clog<<"num_procs: "<<num_procs<<endl;	

	if( freopen(output, "w", stdout) == NULL )
		report_exit("Ouptut file error\n");

	Circuit::set_parameters(epsilon, omega, overlap_ratio, 
			max_block_nodes, mode);
	// start to parfile
	vector<Circuit *> cktlist;
	MPI_CLASS mpi_class;
	// allocate class Tran for all cores
	Tran tran;
	if(my_id==0){
		mpi_class.start_task = new int [num_procs];
		mpi_class.end_task = new int [num_procs];
		mpi_class.tasks_n = new int [num_procs];
		mpi_class.MPI_Assign_Task(num_procs);
	}
	MPI_Scatter(mpi_class.tasks_n, 1, MPI_INT, 
		&mpi_class.block_size, 
		1, MPI_INT, 0, MPI_COMM_WORLD);
		
	Parser parser(&cktlist);
	clock_t t1,t2;
	t1=clock();
	parser.parse(my_id, input, mpi_class, tran, num_procs);
	MPI_Barrier(MPI_COMM_WORLD);

	// after parsing, this mem can be released
	t2=clock();
	if(my_id==0) clog<<"Parse time="<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	
	double mpi_t11, mpi_t12;
	mpi_t11 = MPI_Wtime();
	
	for(size_t i=1;i<cktlist.size();i++){
		Circuit * ckt = cktlist[i];
		if(my_id==0){
			clog<<"<======== solving: "<<ckt->get_name()<<" =========>"<<my_id<<endl;
		}
		ckt->solve(my_id, num_procs, mpi_class, 
				tran);	
		free(ckt);

		// need to print out nodes for all cores
		/*if(my_id==0)
	 		tran.print_tr_nodes();*/
		//clog<<"before barrier: "<<my_id<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	
	mpi_t12 = MPI_Wtime();
	
	// output a single ground node
	if(my_id==0){
		printf("G  %.5e\n", 0.0);
		clog<<"solve using: "<<1.0*(mpi_t12-mpi_t11)<<endl;
		close_logfile();
		//clog<<"after close file/ "<<endl;
	}

	// use too much memory, need to check
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_id==0) clog<<"after mpi barrier. "<<endl;
	MPI_Finalize();
	clog<<"after finalize. "<<endl;
	//close_logfile();
	//cout<<"close logfile. "<<endl;
	return 0;
}
