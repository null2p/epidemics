#include <cstdio>
//#include <memory.h>
#include <random>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mpi.h>

using namespace std;

void Initialization(unsigned L, unsigned &L_A, unsigned &L_B, vector<int> &A_arg, vector<int> &A, vector<int> &B_arg, vector<double> &B_time);
void ABdynamics(unsigned L, unsigned &L_A, unsigned &L_B,std::mt19937 &generator,  std::uniform_real_distribution<double> &mt_uniform,  vector<int> &A, vector<int> &A_arg, vector<int> &B_arg, vector<double> &B_time, double p_lambda, double &t ,float theta,float D_b, float R_btoa);

int main(int argc, char* argv[]) {

    MPI::Init();
    unsigned WORKER = static_cast<unsigned>(MPI::COMM_WORLD.Get_rank());
    unsigned NUM_WORKERS = static_cast<unsigned>(MPI::COMM_WORLD.Get_size());
    // NUM_WORKERS == ENSENBLE SIZE

    unsigned x, y, dN, sample,L_A = 0,L_B=0, L = atoi(argv[1]);
    unsigned init_T = atoi(argv[5]), final_T = atoi(argv[6]),surv_T=100;
    int * n, * surv_seed, *global_surv_seed;
    surv_seed = (int *) calloc(NUM_WORKERS, sizeof(int)); 
    global_surv_seed = (int *) calloc(NUM_WORKERS, sizeof(int)); 
    float theta = atof(argv[3]), R_btoa = atof(argv[4]), D_b=atof(argv[7]);
    double dT, p, q,t,p_lambda = atof(argv[2]);

    vector<int> A_arg;
    vector<int> A;
    vector<int> B_arg;
    vector<double> B_time;
    //vector<int> zero_arg;
    //double rho_A = 0, rho_AB = 0, survP = 0; 
    Initialization(L, L_A, L_B, A_arg, A, B_arg, B_time);


    /////////////////////////////////////////////////////////////
    //this 'while' will find seed making survived 'A' in surv_T// 
    unsigned seed = 100+2*WORKER;
    std::uniform_real_distribution<double> mt_uniform(0.0,1.0);
    std::mt19937 generator (100+2*WORKER);

    unsigned steps = 0;
    
    double test_t = 0;
    while(test_t<=surv_T){
        //dN=0;
        ABdynamics(L, L_A, L_B, generator, mt_uniform, A, A_arg, B_arg, B_time, p_lambda, test_t, theta,D_b, R_btoa);

        if(L_A == 0&steps>10) {
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

        Initialization(L, L_A, L_B, A_arg, A, B_arg, B_time);
        cout<<"selected WORKER : "<<WORKER<<endl;

        ////////////////////////////////////////////
        //file opening for saving patterns and tau//
        ofstream producefile,taufile;
        std::stringstream pdname,tauname;
        //pdname<<"./patternsCPTF_ens"<<WORKER<<".txt";
        pdname<<"./seed1/theta"<<theta<<"/pattern_L"<<L<<"T_i"<<init_T<<"T_f"<<final_T<<"ens"<<WORKER<<"p"<<p_lambda<<"theta"<<theta<<"R"<<R_btoa<<"Db"<<D_b<<".txt";
        //tauname<<"./seed1/theta"<<theta<<"/tau_L"<<L<<"T_i"<<init_T<<"T_f"<<final_T<<"ens"<<WORKER<<"p"<<p_lambda<<"theta"<<theta<<"R"<<R_btoa<<"Db"<<D_b<<".txt";
        producefile.open(pdname.str().c_str());
        //taufile.open(tauname.str().c_str());
        /////////////////////////////////////////////

        std::mt19937 reseted_gen (100+2*WORKER);
        t= 0 ; steps = 0 ;

        while(steps<final_T){
            /*
            if(WORKER ==0){ 
               cout<<"\nSAMPLE : "<<sample <<", t : "<<t<<"========================="<<100*t/(double)final_T<<"%\r";
            }
            */
            ABdynamics(L, L_A, L_B, reseted_gen, mt_uniform, A, A_arg, B_arg, B_time, p_lambda, t, theta,D_b, R_btoa);

            if(steps > init_T){
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
		free(n);
		producefile << "\n";
		//taufile << "\n"; be careful here, it could be hard to read file. # of tau is frequently different 
            } 

            steps++;
        }
 
        producefile.close();
        //taufile.close();
    }
    ///////////////////////////////////////////////////////////

    MPI_Barrier(MPI_COMM_WORLD);
    MPI::Finalize();

    return 0;
}


void Initialization(unsigned L, unsigned &L_A, unsigned &L_B, vector<int> &A_arg, vector<int> &A, vector<int>& B_arg, vector<double>& B_time){

    A_arg.clear();
    A.clear();
    B_arg.clear();
    B_time.clear();
    L_A = 0; L_B = 0;

    //initialization 
    for(unsigned i=0;i<L;i++){
        if(i==L/2) { // single seed
        //if(i%2==0) { // half seed
            A_arg.push_back(i);
            A.push_back(1); 
            L_A++;
        }
        else {
            //zero_arg.push_back(i);
            A.push_back(0);
        }
    }
    //rho_A = L_A/(double)L;
    //rho_AB = (L_A+L_B)/(double)L;
    //survP = 0;       
}


void ABdynamics(unsigned L, unsigned &L_A, unsigned &L_B, std::mt19937 &generator, std::uniform_real_distribution<double> &mt_uniform, vector<int> &A, vector<int> &A_arg, vector<int> &B_arg, vector<double> &B_time, double p_lambda, double &t, float theta,float D_b, float R_btoa)
{
        //Dynamics of A and B
        for(unsigned x=0;x<L_A+L_B;x++){
            float p = mt_uniform(generator);
            unsigned index = p*(L_A+L_B);
            //Dynamics of A
            if(index < L_A){
                double p_A = mt_uniform(generator);
                //A ->B
                if(p_A < p_lambda){
                    B_arg.push_back(A_arg.at(index));

                    double tau = pow((1-mt_uniform(generator)),-1.0/theta);
                    //taufile<<tau<<" ";
                    tau += t;
                    B_time.push_back(tau);

                    A.at(A_arg.at(index))--;
                    A_arg.erase(A_arg.begin()+index);

                    L_B++; L_A--;
                }
                else{
                    // 0A0 -> 0AA or AA0
                    // q set the direction to branch
                    float q = mt_uniform(generator);
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
                    float r = mt_uniform(generator);
                    // B -> A
                    if(r<R_btoa){
                        A_arg.push_back(B_arg.at(index));
                        A.at(B_arg.at(index))++;
                        L_A++;
                    }
                    // if B -> 0, nothing happens but B is gone
                    B_time.erase(B_time.begin()+index);
                    B_arg.erase(B_arg.begin()+index);
                    L_B --;
                }
            }
            if(L_A+L_B>0) t = t + 1.0/(L_A+L_B);
        } 

}

