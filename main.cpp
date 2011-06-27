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
	double mpi_t1, mpi_t2;

	mpi_t1 = MPI_Wtime();
	
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
	if( freopen(output, "w", stdout) == NULL )
		report_exit("Ouptut file error\n");

	Circuit::set_parameters(epsilon, omega, overlap_ratio, 
			max_block_nodes, mode);
	// start to parfile
	vector<Circuit *> cktlist;
	Parser parser(&cktlist);
	clock_t t1,t2;
	t1=clock();
	parser.parse(input);
	t2=clock();
	//if(my_id==0) clog<<"Parse time="<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	//if( cktlist.size()>0 ) cktlist[0]->check_sys();
	
	// do the job
	//clog<<"number of layers: "<<Circuit::get_total_num_layer()<<endl;
	//if( mode == 0 ) clog<<"Solve using block-iterative."<<endl;
	//else clog<<"Solve using direct LU."<<endl;
	
	//clog<<"start mpi init process. "<<endl;
	//return 0;
	int my_id;
	int num_procs;

	int new_rank;
	MPI_Group orig_group, new_group;
	MPI_Comm new_comm;
	MPI_Comm comm;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

	// divide tasks into 2 distinct groups based upon rank
	size_t base = 0;
	int ranks1[12]={0,1,2,3,4,5, 6, 7,8,9,10,11};
	int ranks2[12]={12,13,14,15,16,17,18,19,20,21, 22, 23};
	/*int ** ranks;
	ranks = new int *[cktlist.size()];
	for(size_t i=0;i<cktlist.size();i++)
		ranks[i]= new int [num_procs/cktlist.size()];
	//clog<<"start to assign. "<<endl;
	for(size_t i=0;i<cktlist.size();i++){
		for(size_t j=0;j<num_procs/cktlist.size();j++){
			ranks[i][j] = j+base;
			//if(my_id==0)
			//clog<<"j, ranks: "<<j<<" "<<ranks[j]<<endl;
		}
		base += num_procs/cktlist.size();
	}*/
	if(my_id <12)
		MPI_Group_incl(orig_group, num_procs/cktlist.size(), ranks1, &new_group);
	else
		MPI_Group_incl(orig_group, num_procs/cktlist.size(), ranks2, &new_group);
	
	MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
	MPI_Group_rank(new_group, &new_rank);
	//clog<<"old_rank: "<<my_id<<" new_rank: "<<new_rank<<endl;
	
	t1 = clock();
	double mpi_t11, mpi_t12;
	mpi_t11 = MPI_Wtime();

	size_t i=0;
	if(my_id<12) i = 0;
	else	i=1;
	Circuit * ckt = cktlist[i];
	if(new_rank ==0)
		clog<<"Solving "<<ckt->get_name()<<endl;
	int new_num_procs = num_procs / cktlist.size();
	ckt->solve(new_rank, new_num_procs, new_comm);
	// DEBUG: output each circuit to separate file
	//char ofname[MAX_BUF];
	//sprintf(ofname,"%s.%s",filename,ckt->get_name().c_str());
	//freopen(ofname,"w", stdout);
	if(my_id ==0)
		cktlist[i]->print();
	//clog<<(*ckt)<<endl;
	clog<<endl;
	// after that, this circuit can be released
	delete ckt;
	t2 = clock();
	mpi_t12 = MPI_Wtime();
	// sync all the processors

	if(my_id ==0)
	clog<<"solve using: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<" "<<1.0*(mpi_t12-mpi_t11)<<endl;
	
	//fclose(stdout);
	// output a single ground node
	printf("G  %.5e\n", 0.0);

	close_logfile();

	mpi_t2 = MPI_Wtime();
	if(my_id ==0)	
	clog<<"mpi time for my_id is: "<<my_id<<" "<<1.0*(mpi_t2-mpi_t1)<<endl;
	MPI_Finalize();
	
	return 0;
}
