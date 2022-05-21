#include <cstdio>
//#include <memory.h>
//#include <cmath>
#include <random>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mpi.h>
#include <cblas.h>
#include <cstring>
#include <string>

#include "Random123/philox.h"
#include "Random123/examples/uniform.hpp"

typedef r123::Philox4x64_R<7> RNG;

using namespace std;

class Params
{
    public:
	unsigned L, parallelization;
        double p_lambda;
        float theta, R_btoa, newR;
        unsigned T, midT;
        float D_b;
        char * seed_type;
};

void Initialization(unsigned L, unsigned &L_A, unsigned &L_B, vector<int> &A_arg, vector<int> &A, vector<int> &B,vector<int> &B_arg, vector<double> &B_time, char * seed_type);
void openFile(ofstream & producefile, unsigned NUM_WORKERS, Params & params);
void ABdynamics(Params& params, unsigned &L_A, unsigned &L_B, RNG::key_type &k, RNG::ctr_type &c, vector<int> &A, vector<int> &B,vector<int> &A_arg, vector<int> &B_arg, vector<double> &B_time, double &t);
int ReturnClusterSize(vector<int> A,unsigned L, unsigned N_G, unsigned index, unsigned * visit);
void surface_roughness(unsigned L,const vector<int> &A,const vector<int> &B, const double rho_A, const double rho_AB, double &rough_A, double &rough_B, double &rough_AB);
void RadiusSiteDensity(Params &params,unsigned& L_A, unsigned &L_B,vector<int>&A, vector<int> &B, int time,double* radius_A, double* radius_B, double* density_Asite, double*density_Bsite, double* N_Acluster, double* N_Bcluster);

int main(int argc, char* argv[]) {

	MPI::Init();
	unsigned WORKER = static_cast<unsigned>(MPI::COMM_WORLD.Get_rank());
	unsigned NUM_WORKERS = static_cast<unsigned>(MPI::COMM_WORLD.Get_size());
	// NUM_WORKERS == ENSENBLE SIZE


	Params params;
	params.L = atoi(argv[1]);
	params.p_lambda = atof(argv[2]);
	params.theta = atof(argv[3]), params.R_btoa = atof(argv[4]);
	params.T = atoi(argv[5]);
	params.D_b = atof(argv[6]);
	params.seed_type = argv[7];
	params.parallelization = atoi(argv[8]);

	int run_by_onecore = params.parallelization/NUM_WORKERS;

	unsigned L_A = 0, L_B = 0;
	vector<int> A_arg;
	vector<int> A;
	vector<int> B_arg;
	vector<int> B(params.L,0);
	vector<double> B_time;
	vector<int> zero_arg;
	double* rho_A , *rho_A_active, *rho_B ,*rho_B_active, *survP , *radius_A, *radius_B, *N_Acluster, *N_Bcluster; 

	int measure_T = 10 + ((std::log(params.T)/std::log(10))-1)*30 ;
	rho_A = new double[measure_T]{};
	rho_A_active = new double[measure_T]{};
	rho_B = new double[measure_T]{};
	rho_B_active = new double[measure_T]{};
	survP = new double[measure_T]{};
	radius_A = new double[measure_T]{};
	radius_B = new double[measure_T]{};
	N_Acluster = new double[measure_T]{};
	N_Bcluster = new double[measure_T]{};
	//rhogh_A = new double[measure_T]{};
	//rhogh_B = new double[measure_T]{};
	//rhogh_AB = new double[measure_T]{};

	ofstream producefile, sp_corrfile, tp_corrfile;
	if(WORKER == 0)  openFile(producefile, NUM_WORKERS, params);

	//cout<<"measure_T : "<<measure_T<<endl;

	for(int i = 0; i<run_by_onecore;i++){
		
		Initialization(params.L, L_A, L_B, A_arg, A, B, B_arg, B_time, params.seed_type);

		//rho_A[0] += L_A/(double)params.L; rho_A_active[0] += L_A/(double)params.L; rho_AB[0] = (L_A+L_B)/(double)params.L; rho_AB_active[0] = (L_A+L_B)/(double)params.L;	survP[0] = 1;       

		double t=0;
		unsigned steps = 0;
		int measure_time = 1, measure_count =0;
		while(t<params.T){
			//if(WORKER ==0){ 
			//    cout<<"STEPS : "<<steps <<", t : "<<t<<"========================="<<100*t/(double)params.T<<"%\r";
			//}
			RNG::key_type k = {{10*WORKER,steps*2}};
			RNG::ctr_type c = {{0,10,WORKER+i*NUM_WORKERS,1000}};

			ABdynamics(params, L_A, L_B,k,c, A, B,A_arg, B_arg, B_time, t);

			if(L_A+L_B == 0)  break;
			//if(L_A+L_B == 0){cout<<"\ndead\n"; break;}

			//surface_roughness(params.L,A,B, rho_A,rho_AB,rough_A,rough_B,rough_AB);
			if(t >= measure_time){
			
				//if(WORKER==0)  cout<<"\t"<<"time : "<<t<<", measure_time : "<<measure_time<<" "<<", count : "<<measure_count<<endl;
				RadiusSiteDensity(params,L_A,L_B,A,B,measure_count,radius_A,radius_B,rho_A_active,rho_B_active,N_Acluster,N_Bcluster);

				rho_A[measure_count] += L_A/(double)params.L;
				rho_B[measure_count] += L_B/(double)params.L;
				survP[measure_count] += 1;       

				measure_count++;
				
				if(measure_time<=10)  measure_time = measure_count+1;
				else measure_time = std::pow(10.0,(double)(measure_count+22)/30.0);

			}
			steps ++; 
		
		}
	}

	double* global_rho_A , *global_rho_A_active, *global_rho_B ,*global_rho_B_active, *global_survP , *global_radius_A, *global_radius_B , *global_N_Acluster, *global_N_Bcluster; 

	global_rho_A = new double[measure_T]{};
	global_rho_A_active = new double[measure_T]{};
	global_rho_B = new double[measure_T]{};
	global_rho_B_active = new double[measure_T]{};
	global_survP = new double[measure_T]{};
	global_radius_A = new double[measure_T]{};
	global_radius_B = new double[measure_T]{};
	global_N_Acluster = new double[measure_T]{};
	global_N_Bcluster = new double[measure_T]{};

	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Reduce(rho_A,global_rho_A,measure_T,MPI::DOUBLE,MPI::SUM,0);
	MPI::COMM_WORLD.Reduce(rho_B,global_rho_B,measure_T,MPI::DOUBLE,MPI::SUM,0);

	MPI::COMM_WORLD.Reduce(rho_A_active,global_rho_A_active,measure_T,MPI::DOUBLE,MPI::SUM,0);
	MPI::COMM_WORLD.Reduce(rho_B_active,global_rho_B_active,measure_T,MPI::DOUBLE,MPI::SUM,0);

	MPI::COMM_WORLD.Reduce(radius_A,global_radius_A,measure_T,MPI::DOUBLE,MPI::SUM,0);
	MPI::COMM_WORLD.Reduce(radius_B,global_radius_B,measure_T,MPI::DOUBLE,MPI::SUM,0);
	MPI::COMM_WORLD.Reduce(N_Acluster,global_N_Acluster,measure_T,MPI::DOUBLE,MPI::SUM,0);
	MPI::COMM_WORLD.Reduce(N_Bcluster,global_N_Bcluster,measure_T,MPI::DOUBLE,MPI::SUM,0);
	MPI::COMM_WORLD.Reduce(survP,global_survP,measure_T,MPI::DOUBLE,MPI::SUM,0);
	//MPI::COMM_WORLD.Reduce(&rough_B,&global_rough_B,measure_T,MPI::DOUBLE,MPI::SUM,0);
	//MPI::COMM_WORLD.Reduce(&rough_AB,&global_rough_AB,measure_T,MPI::DOUBLE,MPI::SUM,0);
	//MPI::COMM_WORLD.Reduce(&N_gA,&global_N_gA,measure_T,MPI::UNSIGNED,MPI::SUM,0);
	//MPI::COMM_WORLD.Reduce(&N_gB,&global_N_gB,measure_T,MPI::UNSIGNED,MPI::SUM,0);


	if(WORKER == 0 ){
		int time = 0, parallel = params.parallelization;
		for(int t= 0;t<measure_T;t++){
			if(t<=10) time = t+1;
			else time = std::pow(10.0,(t+22)/30.0);

			if(time > params.T) break;
			if(global_rho_A[t] + global_rho_B[t] == 0) break;
			//cout<<time<<"\t";
			producefile << time <<" "<< global_rho_A[t]/(double)parallel<<" "<< global_rho_A_active[t]/(double)parallel<<" "<<global_rho_B[t]/(double)parallel<<" "<<global_rho_B_active[t]/(double)parallel<<" "<<global_radius_A[t]/(double)parallel<<" "<<global_radius_B[t]/(double)parallel<<" "<<global_N_Acluster[t]/(double)parallel<<" "<<global_N_Bcluster[t]/parallel<<" "<<global_survP[t]/parallel<<endl;
		}

	}

	producefile.close();


	MPI::Finalize();
	return 0;
}


void Initialization(unsigned L, unsigned &L_A, unsigned &L_B, vector<int> &A_arg, vector<int> &A, vector<int> &B,vector<int>& B_arg, vector<double>& B_time, char * seed_type)
{
    A_arg.clear();
    A.clear();
    B.clear();
    B_arg.clear();
    B_time.clear();
    L_A = 0; L_B = 0;

    //initialization
    if(!strcmp(seed_type,"single")){
	for(unsigned i=0;i<L;i++){
       	    if(i==L/2) { // single seed
                A_arg.push_back(i);
                A.push_back(1);
                L_A++;
            }
            else {
                //zero_arg.push_back(i);
                A.push_back(0);
            }
            B.push_back(0);
         }
    }
    else if(!strcmp(seed_type,"half")){
        for(unsigned i=0;i<L;i++){
            if(i%2==0) { // half seed
                A_arg.push_back(i);
                A.push_back(1);
                L_A++;
            }
            else {
                //zero_arg.push_back(i);
                A.push_back(0);
            }
            B.push_back(0);
        }
    }
    //rho_A = L_A/(double)L;
    //rho_AB = (L_A+L_B)/(double)L;
    //survP = 0;
}

void openFile(ofstream & producefile, unsigned NUM_WORKERS, Params & params){

    std::stringstream pdname, sp_corrname, tp_corrname;
    pdname<<"./rho_scaling/seed_"<<params.seed_type<<"/theta"<<params.theta<<"/L"<<params.L<<"T"<<params.T<<"ens"<<params.parallelization<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";

    producefile.open(pdname.str().c_str());
    producefile <<"#t    <rhoA[t]>    <siteA_active[t]>    <rhoB[t]>    <siteB_active[t]>   <radius_A[t]>    <radius_B[t]>  <N_A[t]>    <N_B[t]>   <survP[t]>...  contact process (L="<<params.L<<", T="<<params.T<<", ENS="<<NUM_WORKERS<<", p_lambda = "<<params.p_lambda<<", theta = "<<params.theta<<", R_btoa = "<<params.R_btoa<< endl;


}

void ABdynamics(Params& params, unsigned &L_A, unsigned &L_B, RNG::key_type &k, RNG::ctr_type &c, vector<int> &A, vector<int> &B,vector<int> &A_arg, vector<int> &B_arg, vector<double> &B_time, double &t)
{
	unsigned L  = params.L;
	double p_lambda = params.p_lambda;
	float theta = params.theta;
	float D_b = params.D_b;
	float R_btoa = params.R_btoa;

        RNG rng;
        RNG::ctr_type r;            

        for(unsigned x=0;x<L_A+L_B;x++){
            r = rng(c,k);
            //double p = mt_uniform(generator);
            double p = r123::u01fixedpt<double>(r.v[0]);

            unsigned index = p*(L_A+L_B);
            //if(p==1) index = L_A+L_B-1;
            //Dynamics of A
            if(index < L_A){
                //double p_A = mt_uniform(generator);
                double p_A = r123::u01fixedpt<double>(r.v[1]);
                //A ->B
                if(p_A < p_lambda){
                    B_arg.push_back(A_arg.at(index));

                    //double tau = pow((1-mt_uniform(generator)),-1.0/theta) + t;
                    double tau = pow( (1-r123::u01fixedpt<double>(r.v[2]) ),-1.0/theta) + t;
                    B_time.push_back(tau);

                    A.at(A_arg.at(index))--;  B.at(A_arg.at(index))++;
                    A_arg.erase(A_arg.begin()+index);

                    L_B++;  L_A--;
                }
                else{
                    // 0A0 -> 0AA or AA0
                    // q set the direction to branch
                    //double q = mt_uniform(generator);
                    double q = r123::u01fixedpt<double>(r.v[2]);
		    unsigned neighbor_idx = 0;
                    //right
                    if( q < D_b )  neighbor_idx = (A_arg.at(index)+L+1)%L;
                    //left
                    else  neighbor_idx = (A_arg.at(index)+L-1)%L;

                    if(A.at(neighbor_idx)==0){
                            A_arg.push_back(neighbor_idx);
                            A.at(neighbor_idx) ++;  L_A ++;

                            //double tmp_R_sq = neighbor_idx - (float)L/2;
                            //if (tmp_R_sq > R_sq) R_sq = tmp_R_sq;
                    }
                }
            }
            //B dynamics
            else{
                index = index - L_A ;
                if(B_time.at(index)<t){
                    //float r = mt_uniform(generator);
                    double r_db = r123::u01fixedpt<double>(r.v[3]);
                    // B -> A
                    if(r_db<R_btoa){
                        A_arg.push_back(B_arg.at(index));
                        A.at(B_arg.at(index))++;
                        L_A++;
                    }
                    // if B -> 0, nothing happens but B is gone
                    B.at(B_arg.at(index))--;

                    B_arg.erase(B_arg.begin()+index);
                    B_time.erase(B_time.begin()+index);
                    L_B --;
                }
            }
            if(L_A+L_B!=0)t = t + 1.0/(float)(L_A+L_B);
            ++c[0];
        }   


}



int ReturnClusterSize(vector<int> A,unsigned L, unsigned N_G, unsigned index, unsigned * visit){

    int rightstop=0, leftstop=0;
    for(unsigned j = 0; j < L-1; j++){
        int right = (index + j)%L;
        int left = (index - j - 1 + L)%L;
        if(A[index]*A[right]==1&visit[right]==0&rightstop==0) {
            visit[right] = 1;
            N_G++;
        }
        else rightstop=1;

        if(A[index]*A[left]==1&visit[left]==0&leftstop==0) {
            visit[left] = 1;
            N_G++;
        }
        else leftstop=1;

        if(rightstop*leftstop)  break;
    }

    return N_G;
}

//!! the result is squared roughenss !! w^2
void surface_roughness(unsigned L,const vector<int> &A,const vector<int> &B, const double rho_A, const double rho_AB, double &rough_A, double &rough_B, double &rough_AB)
{

    vector<double> heightA(L,0);
    vector<double> heightB(L,0);
    vector<double> heightAB(L,0);

    heightA[0] = A[0] - rho_A;  heightB[0] = B[0] - rho_AB +rho_A;  heightAB[0] = A[0]+B[0] - rho_AB;
    
    double heightA_avg =0, heightB_avg =0, heightAB_avg =0;
    rough_A =0; rough_B =0; rough_AB =0;
    
    for(unsigned i = 1; i<L;i++){
        heightA[i] = heightA[i-1] + A[i] - rho_A;
        heightB[i] = heightB[i-1] + B[i] - rho_AB + rho_A;
        heightAB[i] = heightAB[i-1] + A[i] + B[i] - rho_AB;
    }

    heightA_avg = accumulate(heightA.begin(),heightA.end(),0.);
    heightB_avg = accumulate(heightB.begin(),heightB.end(),0.);
    heightAB_avg = accumulate(heightAB.begin(),heightAB.end(),0.);

    heightA_avg /= L;
    heightB_avg /= L;
    heightAB_avg /= L;

    for(unsigned i =0;i<L;i++){
        rough_A += (heightA[i]-heightA_avg)*(heightA[i]-heightA_avg);
        rough_B += (heightB[i]-heightB_avg)*(heightB[i]-heightB_avg);
        rough_AB += (heightAB[i]-heightAB_avg)*(heightAB[i]-heightAB_avg);
    }

    rough_A /= L;
    rough_B /= L;
    rough_AB /= L;

}

void RadiusSiteDensity(Params &params,unsigned& L_A, unsigned &L_B,vector<int>&A, vector<int> &B, int time,double* radius_A, double* radius_B, double* density_Asite, double*density_Bsite, double* N_Acluster, double* N_Bcluster){

        int minA_idx = params.L-1, maxA_idx = 0, minB_idx = params.L-1, maxB_idx = 0, NAsite =0, NBsite=0;
        //double density_Bsite = 0, density_Asite = 0, N_A_cluster =0, N_B_cluster = 0, radius_A = 0, radius_B = 0;

        for(int i =0; i<params.L;i++){
                if(B[i]!=0){
                        NBsite ++;
                        if(i < minB_idx) minB_idx = i;
                        if(i > maxB_idx) maxB_idx = i;
                }
                if(A[i]!=0){
                        NAsite ++;
                        if(i < minA_idx) minA_idx = i;
                        if(i > maxA_idx) maxA_idx = i;
                }

        }
        //if (density_Bsite == 0)  return ;

        if(maxA_idx >= minA_idx){
                N_Acluster[time] += (double)L_A/(maxA_idx - minA_idx+1);
                radius_A[time] += (maxA_idx-minA_idx)*(maxA_idx-minA_idx);
        }

        if(maxB_idx >= minB_idx){
                N_Bcluster[time] += (double)L_B/(maxB_idx - minB_idx+1);
                radius_B[time] += (maxB_idx-minB_idx)*(maxB_idx-minB_idx);
        }

        density_Asite[time] += NAsite/(double)params.L;
        density_Bsite[time] += NBsite/(double)params.L;

}

