#include <cstdio>
//#include <memory.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mpi.h>
#include <cmath>

#include "Random123/philox.h"
#include "Random123/examples/uniform.hpp"

typedef r123::Philox4x64_R<7> RNG;

using namespace std;

class Params
{
    public:
        unsigned L;
        double p_lambda;
        float theta, R_btoa;
        unsigned init_T, final_T;
        float D_b;
};

void Initialization(unsigned L, unsigned &L_A, unsigned &L_B, vector<int> &A_arg, vector<int> &A, vector<int> &B,vector<int> &B_arg, vector<double> &B_time);
void openFile(ofstream& Afile, ofstream& Bfile, ofstream& Sfile, vector<int> &A, vector<int> &B,const int WORKER,Params& params);
void ABdynamics(Params & params, unsigned WORKER, unsigned steps, unsigned &L_A, unsigned &L_B, RNG::key_type &k, RNG::ctr_type &c, vector<int> &A, vector<int> &B, vector<int> &A_arg, vector<int> &B_arg, vector<double> &B_time, double &t);
void RadiusSiteDensity(Params &params,unsigned& L_A, unsigned &L_B,vector<int>&A, vector<int> &B, int time, ofstream& Sfile);

int main(int argc, char* argv[]) {

    MPI_Init(&argc,&argv);
    int WORKER,NUM_WORKERS;
    MPI_Comm_rank(MPI_COMM_WORLD, &WORKER);
    MPI_Comm_size(MPI_COMM_WORLD, &NUM_WORKERS);

    // NUM_WORKERS == ENSENBLE SIZE

    Params params;
    params.L = atoi(argv[1]);
    params.p_lambda = atof(argv[2]);
    params.theta = atof(argv[3]), params.R_btoa = atof(argv[4]);
    params.init_T = atoi(argv[5]), params.final_T = atoi(argv[6]);
    params.D_b = atof(argv[7]);

    unsigned argv_seed = atoi(argv[8]);
    unsigned sample,L_A = 0,L_B=0;
    int * n, * surv_seed, *global_surv_seed;
    surv_seed = (int *) calloc(NUM_WORKERS, sizeof(int)); 
    global_surv_seed = (int *) calloc(NUM_WORKERS, sizeof(int)); 

    vector<int> A_arg;
    vector<int> A;
    vector<int> B_arg;
    vector<int> B(params.L, 0);
    vector<double> B_time;
    //vector<int> zero_arg;
    //double rho_A = 0, rho_AB = 0, survP = 0; 
    Initialization(params.L, L_A, L_B, A_arg, A,B, B_arg, B_time);


    /////////////////////////////////////////////////////////////
    //this 'while' will find seed making survived 'A' in surv_T// 
    unsigned seed = 100+2*WORKER;
    //std::uniform_real_distribution<double> mt_uniform(0.0,1.0);
    //std::mt19937 generator (100+2*WORKER);

    RNG::key_type k = {{2*WORKER,argv_seed}};
    RNG::ctr_type c = {{0,10,100,1000}};
    unsigned steps = 0;
    double test_t = 0;

    unsigned surv_T=100;
    while(test_t<=surv_T){
        //dN=0;
        ABdynamics(params, WORKER,steps, L_A, L_B, k, c, A,B, A_arg, B_arg, B_time,test_t);

        if(L_A == 0&steps>10) {
	    //cout<<"WORKER"<<WORKER<<" dead at "<<steps<<endl;
            break;
        }
        steps++;
    }
    ///////////////////////////////////////////////////////////////

    //make survived array and sharing for all WORKERS
    if( L_A != 0 )  surv_seed[WORKER] = 1;
    MPI_Allreduce(surv_seed,global_surv_seed,NUM_WORKERS,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    //////////////////////////////////////////////////////
    //receiving order to use only limited number of seed//  
    unsigned survived_order = 0;
    int order = 0;

    int *surv_arg;
    surv_arg = (int *)calloc(4,sizeof(int));
    
    for(int i=0;i<4;i++)  surv_arg[i] = -1;

    if(WORKER==0){
        for(unsigned i =0; i<NUM_WORKERS;i++){
            if(global_surv_seed[i]!=0) {
                surv_arg[order] = i;
                order ++;
            }
            if(order > 3) break;
        }
    }
    MPI_Bcast(surv_arg,4,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////
    //Init main function with SURVIVED & SELECTED 4 seeds//

    if(WORKER == surv_arg[0]|WORKER == surv_arg[1]|WORKER == surv_arg[2]|WORKER == surv_arg[3]){

        Initialization(params.L, L_A, L_B, A_arg, A,B, B_arg, B_time);
        cout<<"selected WORKER : "<<WORKER<<endl;

        ////////////////////////////////////////////
        //file opening for saving patterns and tau//
        ofstream Afile,Bfile,Sfile,taufile;
	openFile(Afile,Bfile, Sfile, A,B,WORKER,params);
        //std::stringstream pdname,tauname;
        //pdname<<"./patternsCPTF_ens"<<WORKER<<".txt";
        //pdname<<"./seed_single/theta"<<theta<<"/pattern_L"<<L<<"T_i"<<init_T<<"T_f"<<final_T<<"ens"<<WORKER<<"p"<<p_lambda<<"theta"<<theta<<"R"<<R_btoa<<"Db"<<D_b<<".txt";
        //tauname<<"./seed1/theta"<<theta<<"/tau_L"<<L<<"T_i"<<init_T<<"T_f"<<final_T<<"ens"<<WORKER<<"p"<<p_lambda<<"theta"<<theta<<"R"<<R_btoa<<"Db"<<D_b<<".txt";
        //producefile.open(pdname.str().c_str());
        //taufile.open(tauname.str().c_str());
        /////////////////////////////////////////////

        RNG::ctr_type reseted_c = {{0,10,100,1000}};
        double t= 0 ; steps = 0 ;
	unsigned measure_time = 0;

        while(steps < params.final_T){
            /*
            if(WORKER ==0){ 
               cout<<"\nSAMPLE : "<<sample <<", t : "<<t<<"========================="<<100*t/(double)final_T<<"%\r";
            }
            */
            ABdynamics(params, WORKER,steps, L_A, L_B, k, reseted_c, A,B, A_arg, B_arg, B_time,t);

	    if(steps >= params.init_T){
/*
		n = (int *) calloc(L,sizeof(int));
		for(unsigned j=0;j<B_arg.size();j++) n[B_arg.at(j)] =-1;

		//if else is divided for decreasing of visiting time
		if(A_arg.size() < L){
		    for(unsigned j=0;j<A_arg.size();j++){
			//dN+=n[j];
			unsigned A_idx = A_arg.at(j);
			if(A.at(A_idx) > 1) n[A_idx] = 2;
			else if(A.at(A_idx)==1) n[A_idx] = 1;
		    }
		}
		else{
		    for(unsigned j=0;j<L;j++){
			//dN+=n[j];
			if(A.at(j) > 1) n[j] = 2;
			else if(A.at(j)==1) n[j] = 1;
		    }
		}

		for(unsigned j=0;j<L;j++) producefile << n[j] << " ";
		producefile << "\n";
		free(n);
*/

		    if(t >= measure_time){
			RadiusSiteDensity(params,L_A,L_B,A,B,measure_time,Sfile);
			for(unsigned j=0;j< params.L;j++){
			    Afile << A[j] << " ";
			    Bfile << B[j] << " ";
			}
			Afile<<"\n"; Bfile<<"\n";
			//taufile << "\n"; be careful here, it could be hard to read file. # of tau is frequently different 
		  	measure_time++;
		    }
	    } 
            steps++;
        }
 
        Afile.close();  Bfile.close();
        //producefile.close();
        //taufile.close();
    }
    ///////////////////////////////////////////////////////////

    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}


void Initialization(unsigned L, unsigned &L_A, unsigned &L_B, vector<int> &A_arg, vector<int> &A, vector<int> &B, vector<int>& B_arg, vector<double>& B_time)
{
    A_arg.clear();
    A.clear();
    B.clear();

    B_arg.clear();
    B_time.clear();
    L_A = 0; L_B = 0;

    //initialization 
    for(unsigned i=0;i<L;i++){
        //if((i==1*L/5)||(i==2*L/5)||(i==3*L/5)||(i==4*L/5)) { // single seed
	if(i==L/2){
	//if(i%2==0) { // half seed
            A_arg.push_back(i);
            A.push_back(1); 
            L_A++;
            
	    B.push_back(0); 
        }
        else {
            //zero_arg.push_back(i);
            A.push_back(0);
	    
	    B.push_back(0); 
        }
    }
    //rho_A = L_A/(double)L;
    //rho_AB = (L_A+L_B)/(double)L;
    //survP = 0;       
}

void openFile(ofstream& Afile, ofstream& Bfile,ofstream& Sfile, vector<int> &A, vector<int> &B,const int WORKER,Params& params)
{

    std::stringstream Aname, Bname, Sname;

    Aname<<"./seed_single/theta"<<params.theta<<"/Aparticle_L"<<params.L<<"T"<<params.final_T<<"ens"<<WORKER<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";
    Bname<<"./seed_single/theta"<<params.theta<<"/Bparticle_L"<<params.L<<"T"<<params.final_T<<"ens"<<WORKER<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";
    Sname<<"./seed_single/theta"<<params.theta<<"/Scailing_L"<<params.L<<"T"<<params.final_T<<"ens"<<WORKER<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";

    //Aname<<"./seed_four/theta"<<params.theta<<"/Aparticle_L"<<params.L<<"T"<<params.final_T<<"ens"<<WORKER<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";
    //Bname<<"./seed_four/theta"<<params.theta<<"/Bparticle_L"<<params.L<<"T"<<params.final_T<<"ens"<<WORKER<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";
    //Sname<<"./seed_four/theta"<<params.theta<<"/Scailing_L"<<params.L<<"T"<<params.final_T<<"ens"<<WORKER<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";

    Afile.open(Aname.str().c_str());  Bfile.open(Bname.str().c_str());  Sfile.open(Sname.str().c_str());
    Afile <<"#    <A[x=0]>    <A[x=1]>  ...  <A[x=L-1]> ...  contact process (L="<<params.L<<", T="<<params.final_T<<", ENS="<<WORKER<<", p_lambda = "<<params.p_lambda<<", theta = "<<params.theta<<", R_btoa = "<<params.R_btoa<< endl;
    Bfile <<"#    <B[x=0]>    <B[x=1]>  ...  <B[x=L-1]> ...  contact process (L="<<params.L<<", T="<<params.final_T<<", ENS="<<WORKER<<", p_lambda = "<<params.p_lambda<<", theta = "<<params.theta<<", R_btoa = "<<params.R_btoa<< endl;
    Sfile<<"#t\tdensity_A\tdensity_B\tP_survival\tdensity_Asite\tdensity_Bsite\tradius_A\tradius_B\tN_Acluster\tN_Bcluster"<<endl;

    for(unsigned i=0; i<params.L;i++) {
        Afile <<" "<<A[i]; Bfile<< " " <<B[i];
    }
    Afile << endl;  Bfile<<endl;

}



void ABdynamics(Params & params, unsigned WORKER, unsigned steps, unsigned &L_A, unsigned &L_B, RNG::key_type &k, RNG::ctr_type &c, vector<int> &A, vector<int> &B, vector<int> &A_arg, vector<int> &B_arg, vector<double> &B_time, double &t)
{
	unsigned L = params.L;
        double p_lambda = params.p_lambda;
        float theta = params.theta;
        float D_b = params.D_b;
        float R_btoa = params.R_btoa;

	RNG rng;
	RNG::ctr_type r;

  //Dynamics of A and B
        for(unsigned x=0;x<L_A+L_B;x++){
	    r = rng(c,k);
            double p = r123::u01fixedpt<double>(r.v[0]);
            unsigned index = p*(L_A+L_B);
	    //Dynamics of A
            if(index < L_A){

                double p_A = r123::u01fixedpt<double>(r.v[1]);
                //A ->B
                if(p_A < p_lambda){
                    B_arg.push_back(A_arg.at(index));

                    double tau = pow((1-r123::u01fixedpt<double>(r.v[2])),-1.0/theta);
                    //taufile<<tau<<" ";
                    tau += t;
                    B_time.push_back(tau);

                    A.at(A_arg.at(index))--;  B.at(A_arg.at(index))++;
                    A_arg.erase(A_arg.begin()+index);

                    L_B++; L_A--;
                }
                else{
                    // 0A0 -> 0AA or AA0
                    // q set the direction to branch
                    double q = r123::u01fixedpt<double>(r.v[2]);
                    //right
                    if( q < D_b ) {
                        unsigned y=(A_arg.at(index)+L+1)%L;
                        if(A.at(y)==0){
                            A_arg.push_back(y);
                            A.at(y) ++;  L_A ++;
                        }
                    }
                    else{
                        //left
                        unsigned y=(A_arg.at(index)+L-1)%L;
                        if(A.at(y)==0){
                            A_arg.push_back(y);
                            A.at(y) ++;  L_A ++;
                        }
                    }
                } 
            }
            //B dynamics
            else{
                index = index - L_A ;
                if(B_time.at(index)<t){
                    double r_db = r123::u01fixedpt<double>(r.v[3]);
                    // B -> A
                    if(r_db < R_btoa){
                        A_arg.push_back(B_arg.at(index));
                        A.at(B_arg.at(index))++;
                        L_A++;
                    }
                    // if B -> 0, nothing happens but B is gone including right above B -> A
		    B.at(B_arg.at(index))--;

                    B_time.erase(B_time.begin()+index);
                    B_arg.erase(B_arg.begin()+index);
                    L_B --;
                }
            }
            if(L_A+L_B>0) t = t + 1.0/(L_A+L_B);
            ++c[0]; 
	} 

 
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

