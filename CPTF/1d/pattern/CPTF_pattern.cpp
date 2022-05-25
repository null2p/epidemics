#include <cstdio>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mpi.h>

#include "Random123/philox.h"
#include "Random123/examples/uniform.hpp"

#include "CPTF.h"

typedef r123::Philox4x64_R<7> RNG;

using namespace std;

class patternParams : public Params
{
    public:
	unsigned init_T,midT,final_T;
	float newR;
};

void Initialization(unsigned L, unsigned &L_A, unsigned &L_B, vector<int> &A_arg, vector<int> &A, vector<int> &B,vector<int> &B_arg, vector<double> &B_time, char * seed_type);
void openFile(ofstream& Afile, ofstream& Bfile,ofstream& Sfile, vector<int> &A, vector<int> &B,const int WORKER,patternParams& params);
//void ABdynamics(unsigned L, unsigned &L_A, unsigned &L_B,RNG::key_type &k , RNG::ctr_type &c,  vector<int> &A, vector<int> &B,vector<int> &A_arg, vector<int> &B_arg, vector<double> &B_time, double p_lambda, double &t ,float theta,float D_b, float R_btoa,double &R_sq);
void ABdynamics(Params & params, unsigned &L_A, unsigned &L_B, RNG::key_type &k, RNG::ctr_type &c, vector<int> &A, vector<int> &B, vector<int> &A_arg, vector<int> &B_arg, vector<double> &B_time, double &t, double& R_sq);
void RadiusSiteDensity(Params &params,unsigned& L_A, unsigned &L_B,vector<int>&A, vector<int> &B, int time, ofstream& Sfile);

int main(int argc, char* argv[]) {

    MPI::Init();
    unsigned WORKER = static_cast<unsigned>(MPI::COMM_WORLD.Get_rank());
    unsigned NUM_WORKERS = static_cast<unsigned>(MPI::COMM_WORLD.Get_size());
    // NUM_WORKERS == ENSENBLE SIZE

    patternParams params;
    params.L = atoi(argv[1]);
    params.p_lambda = atof(argv[2]);
    params.theta = atof(argv[3]), params.R_btoa = atof(argv[4]), params.newR = atof(argv[5]);
    params.init_T = atoi(argv[6]),params.midT = atoi(argv[7]), params.final_T = atoi(argv[8]);
    params.D_b = atof(argv[9]);
    params.seed_type = argv[10];

    //int dt, dN;
    unsigned L_A = 0,L_B=0;
    //int * n;
    //double p, q, t;

    vector<int> A_arg;
    vector<int> A;
    vector<int> B(params.L,0);
    vector<int> B_arg;
    vector<double> B_time;
    //vector<int> zero_arg;
    double rho_A = 0, rho_AB = 0, survP = 0, R_sq = 0; 

    double t=0;

    Initialization(params.L, L_A, L_B, A_arg, A, B, B_arg, B_time, params.seed_type);

    rho_A = L_A/(double)params.L;
    rho_AB = (L_A+L_B)/(double)params.L;
    survP = 0; 
 
    ofstream Afile,Bfile,Sfile;
    openFile(Afile,Bfile,Sfile,A,B,WORKER,params);

    unsigned steps = 0, measure_time = 0;
    while(steps<params.final_T){
    /*
        if(WORKER ==0){ 
           cout<<"\nSAMPLE : "<<sample <<", t : "<<t<<"========================="<<100*t/(double)final_T<<"%\r";
        }
    */
        //dN=0;

	if(steps == params.midT)  params.R_btoa = params.newR;
	
	RNG::key_type k = {{10*WORKER,steps*2}};
        RNG::ctr_type c = {{0,10,100,1000}};

        //ABdynamics(L, L_A, L_B, k, c, A, B,A_arg, B_arg, B_time, p_lambda, t, theta,D_b, R_btoa, R_sq);
        ABdynamics(params, L_A, L_B, k, c, A,B, A_arg, B_arg, B_time,t,R_sq);

        if(steps >= params.init_T){
		if( t >= measure_time ){
			RadiusSiteDensity(params, L_A, L_B,A,B, measure_time,Sfile);
			for(unsigned j=0;j<params.L;j++){
			    Afile<<A[j]<<" ";
			    Bfile<<B[j]<<" ";
			}
			Afile<<"\n"; Bfile<<"\n";
			measure_time ++;
		}
        }
	steps++;
    } 

    Afile.close(); Bfile.close(); Sfile.close();
    //producefile.close();
    //taufile.close();
    MPI::Finalize();
    return 0;
}


void openFile(ofstream& Afile, ofstream& Bfile,ofstream&Sfile, vector<int> &A, vector<int> &B,const int WORKER,patternParams& params)
{

    std::stringstream Aname, Bname, Sname;

    Aname<<"./seed_"<<params.seed_type<<"/theta"<<params.theta<<"/Aparticle_L"<<params.L<<"T"<<params.final_T<<"ens"<<WORKER<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";
    Bname<<"./seed_"<<params.seed_type<<"/theta"<<params.theta<<"/Bparticle_L"<<params.L<<"T"<<params.final_T<<"ens"<<WORKER<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";
    Sname<<"./seed_"<<params.seed_type<<"/theta"<<params.theta<<"/Scailing_L"<<params.L<<"T"<<params.final_T<<"ens"<<WORKER<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";

    Afile.open(Aname.str().c_str());  Bfile.open(Bname.str().c_str()); Sfile.open(Sname.str().c_str());
    Afile <<"#    <A[x=0]>    <A[x=1]>  ...  <A[x=L-1]> ...  contact process (L="<<params.L<<", T="<<params.final_T<<", ENS="<<WORKER<<", p_lambda = "<<params.p_lambda<<", theta = "<<params.theta<<", R_btoa = "<<params.R_btoa<< endl;
    Bfile <<"#    <B[x=0]>    <B[x=1]>  ...  <B[x=L-1]> ...  contact process (L="<<params.L<<", T="<<params.final_T<<", ENS="<<WORKER<<", p_lambda = "<<params.p_lambda<<", theta = "<<params.theta<<", R_btoa = "<<params.R_btoa<< endl;
    Sfile<<"#t\tdensity_A\tdensity_B\tP_survival\tdensity_Asite\tdensity_Bsite\tradius_A\tradius_B\tN_Acluster\tN_Bcluster"<<endl;

    for(unsigned i=0; i<params.L;i++) {
        Afile <<" "<<A[i]; Bfile<< " " <<B[i];
    }
    Afile << endl;  Bfile<<endl;

}

void RadiusSiteDensity(Params &params,unsigned& L_A, unsigned &L_B,vector<int>&A, vector<int> &B, int time, ofstream& Sfile){

	int minA_idx = params.L-1, maxA_idx = 0, minB_idx = params.L-1, maxB_idx = 0;
	double density_Bsite = 0, density_Asite = 0, N_A_cluster =0, N_B_cluster = 0, radius_A = 0, radius_B = 0;

	for(int i =0; i<params.L;i++){
		if(B[i]!=0){
			density_Bsite ++;
			if(i < minB_idx) minB_idx = i;
			if(i > maxB_idx) maxB_idx = i;
		}
		if(A[i]!=0){
			density_Asite ++;
			if(i < minA_idx) minA_idx = i;
			if(i > maxA_idx) maxA_idx = i;
		}

	}
	//if (density_Bsite == 0)  return ;

	if(maxA_idx >= minA_idx){
		N_A_cluster = (double)L_A/(maxA_idx - minA_idx+1);
		radius_A = (maxA_idx-minA_idx)*(maxA_idx-minA_idx);
	}

	if(maxB_idx >= minB_idx){
		N_B_cluster = (double)L_B/(maxB_idx - minB_idx+1);
		radius_B = (maxB_idx-minB_idx)*(maxB_idx-minB_idx);
	}

	density_Asite /= params.L;
	density_Bsite /= params.L;

	Sfile<<time<<"\t"<<L_A/(double)params.L<<"\t"<<L_B/(double)params.L<<"\t"<<"1"<<"\t"<<density_Asite<<"\t"<<density_Bsite<<"\t"<<radius_A<<"\t"<<radius_B<<"\t"<<N_A_cluster<<"\t"<<N_B_cluster<<"\n";

}


