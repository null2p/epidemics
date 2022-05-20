#include <iostream>
#include <iomanip>
#include <mpi.h>
#include "../DEP_algorithm.h"

using namespace std;

class Params{
	public:
		int L,z,T_max, init_total_number, final_total_number, steps;
		std::string init_dist;
		Params(int L, int z, int T_max, int init_total_number, int final_total_number, int steps,std::string init_dist): L(L), z(z), T_max(T_max), init_total_number(init_total_number), final_total_number(final_total_number), steps(steps), init_dist(init_dist)
		{
		}
};

int main(int argc, char* argv[]){

	MPI::Init();
	unsigned WORKER = static_cast<unsigned>(MPI::COMM_WORLD.Get_rank());
    	unsigned NUM_WORKERS = static_cast<unsigned>(MPI::COMM_WORLD.Get_size());

	int parallelization = 10020, run_by_onecore = parallelization/NUM_WORKERS;
	int L = atoi(argv[1]), z = 2,T_max = 100000;//1e+7
	int init_total_number= 6.765*L, final_total_number = 0.4*L, steps = 80;
	std::string init = "homogeneous";
	Params params(L, z, T_max, init_total_number, final_total_number, steps, init);

	//int max_steps = params.final_total_number - params.init_total_number;
	//if(params.steps > max_steps) params.steps = max_steps;
	int total_number = params.init_total_number;
	int seed = WORKER;

	//double* density_B = new double[params.T_max]{};
	//double* P_survival = new double[params.T_max]{};
	//double* mpi_density_B = new double[params.T_max]{};
	//double* mpi_P_survival = new double[params.T_max]{};

	int measure_T_max = ((int)std::log10(T_max))*10 - 10+1;

	double* density_B = new double[measure_T_max]{};
	double* P_survival = new double[measure_T_max]{};
	double* mpi_density_B = new double[measure_T_max]{};
	double* mpi_P_survival = new double[measure_T_max]{};


	if(WORKER==0)  cout<<"#t\tdensity_B\tP_survival"<<endl;

	for(int i=0;i<run_by_onecore;i++){
		seed = WORKER + i*NUM_WORKERS;
		DEPmodel1d density_production(params.L,total_number,params.z,params.T_max, params.init_dist, seed);
		//density_production.DanielMaiaAlgorithm(T_max);
		//density_production.GillespieAlgorithm(T_max);
		density_production.NextReactionMethod(T_max);
		for(int t=0;t<measure_T_max;t++)  {
			density_B[t] += density_production.density_B[t]/run_by_onecore;
			P_survival[t] += (double)density_production.survival[t]/run_by_onecore;
			if(density_production.density_B[t] == 0) break;
		}
	}
	MPI::COMM_WORLD.Reduce(density_B, mpi_density_B ,measure_T_max,MPI::DOUBLE,MPI::SUM,0);
	MPI::COMM_WORLD.Reduce(P_survival, mpi_P_survival ,measure_T_max,MPI::DOUBLE,MPI::SUM,0);

	if(WORKER==0){
		for(int t=0;t<measure_T_max;t++){
			double time = std::pow(10,1.0 + t/10.0);
			cout<<time<<"\t"<<mpi_density_B[t]/NUM_WORKERS<<"\t"<<mpi_P_survival[t]/NUM_WORKERS<<endl;
			if(mpi_density_B[t]+mpi_P_survival[t]==0) break;
		}
	}

	MPI::Finalize();
	return 0;

}


