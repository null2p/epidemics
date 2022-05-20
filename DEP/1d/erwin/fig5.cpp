#include <iostream>
#include <iomanip>
#include <mpi.h>
#include "../DEP_algorithm.h"
#include "../DEP_spatiotemporal.h"

using namespace std;

class Params{
	public:
		int L,z,T_max, init_total_number, final_total_number, steps;
		std::string init_dist;

		Params(int L, int z, int T_max, int init_total_number, int final_total_number, int steps, std::string init_dist): L(L), z(z), T_max(T_max), init_total_number(init_total_number), final_total_number(final_total_number), steps(steps), init_dist(init_dist)
		{
		}
};
void WriteFile(int *mpi_A, int *mpi_B,int parallel,Params params){
	std::ofstream Afile, Bfile;
	//Afile.open("./fig5/A_subdiff.txt"); Bfile.open("./fig5/B_subdiff.txt");
	Afile.open("./pattern/rho676/A.txt"); Bfile.open("./pattern/rho676/B.txt");
	int L = params.L, T_max = params.T_max;

	for(int t =0; t<T_max; t++){
		for(int i=0;i<L;i++){
			Afile << (double)mpi_A[t*L + i]/parallel<<"\t";
			Bfile << (double)mpi_B[t*L + i]/parallel<<"\t";
		}
		Afile<<"\n";	
		Bfile<<"\n";
	}

	Afile.close(); Bfile.close();

}

int main(int argc, char* argv[]){

	MPI::Init();
	unsigned WORKER = static_cast<unsigned>(MPI::COMM_WORLD.Get_rank());
    	unsigned NUM_WORKERS = static_cast<unsigned>(MPI::COMM_WORLD.Get_size());

	int parallelization = 1, run_by_onecore = parallelization/NUM_WORKERS;
	int L = 1024, z =2,T_max = 10000;//1e+4
	int init_total_number= 6.765*L, final_total_number = 0.4*L, steps = 80;
	std::string init = "single";

	Params params(L, z, T_max, init_total_number, final_total_number, steps,init);

	int *alltime_A = new int[L*T_max]{};
	int *alltime_B = new int[L*T_max]{};
	int *mpi_A = new int[L*T_max]{};
	int *mpi_B = new int[L*T_max]{};
	int total_number = params.init_total_number;

	for(int i =0; i< run_by_onecore ; i++){
		SpatioTemporalMode pattern(params.L,total_number,params.z,params.T_max,params.init_dist);
		pattern.SpatioTemporalPattern(T_max, alltime_A, alltime_B);
	}

	MPI::COMM_WORLD.Reduce(alltime_A, mpi_A ,L*T_max,MPI::INT,MPI::SUM,0);
	MPI::COMM_WORLD.Reduce(alltime_B, mpi_B ,L*T_max,MPI::INT,MPI::SUM,0);

	if(WORKER==0) WriteFile(mpi_A,mpi_B,parallelization,params);

	MPI::Finalize();
	return 0;
}


