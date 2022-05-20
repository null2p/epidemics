#include <cstdio>
//#include <memory.h>
#include <random>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mpi.h>
//#define L 512
#define T 100000
//#define p_lambda 0.4959
//#define theta 1

using namespace std;

void Diffusion(unsigned WORKER, float D_A, float D_B,unsigned L,vector<int> &A, vector<int> &B,unsigned steps);
void Reaction(unsigned WORKER,unsigned L,vector<int> &A, vector<int> &B,unsigned steps);
void Reactivation(unsigned WORKER, unsigned L, vector<int> &A, vector<int> &B,unsigned steps);
int main(int argc, char* argv[]) {

    MPI::Init();
    unsigned WORKER = static_cast<unsigned>(MPI::COMM_WORLD.Get_rank());
    unsigned NUM_WORKERS = static_cast<unsigned>(MPI::COMM_WORLD.Get_size());
    // NUM_WORKERS == ENSENBLE SIZE

    unsigned i, j, x, y, dt, dN, sample, L = atoi(argv[1]),N_r=0;
    int * n;
    double p, q,t,D_A = 0.5, D_B = 0.5;
    vector<int> A(L);
    vector<int> B(L);
    //vector<int> A_arg;
    //vector<int> B_arg;
    //vector<double> B_time;
    //vector<int> zero_arg;

    double rho_tot = atof(argv[2]), rho_A = 0, rho_AB = 0, survP = 0; 

    t=0;

    //A_arg.clear();
    //A.clear();
    //B_arg.clear();
    //B_time.clear();

    //initialization, half of them are infected 'A' and others are susceptible 'B' 
    std::mt19937 generator (100+2*WORKER);
    std::uniform_real_distribution<double> mt_uniform(0.0,1.0);
    for(i=0;i<L*rho_tot;i++){
        unsigned randindex = L*mt_uniform(generator);
        if(i%2==0) {
            //A_arg.push_back(randindex);
            A.at(randindex)++; 
            //L_A++;
        }
        else {
            //B_arg.push_back(randindex);
            B.at(randindex)++; 
            //L_B++;
        }
    }

    unsigned A_sum = std::accumulate(A.begin(),A.end(),0);
    unsigned B_sum = std::accumulate(B.begin(),B.end(),0);
    rho_A = A_sum/(double)(A_sum+B_sum);
    survP = 0; 


    /////////////////////////////
    //file open for saving data//
    ofstream producefile;
    std::stringstream pdname;
    pdname<<"./rho_scaling/L"<<L<<"T"<<T<<"ens"<<NUM_WORKERS<<"rho"<<rho_tot<<"D_A"<<D_A<<"D_B"<<D_B<<".txt";
    if(WORKER == 0) {
        producefile.open(pdname.str().c_str());
        producefile <<"#t    U[t]    P[t]    Delta[t]    U_active[t]    P_active[t]    Delta_active[t] ...  Diffusive Epidemic Process (L="<<L<<", T="<<T<<", ENS="<<NUM_WORKERS<< endl;
        //producefile << t<<" "<< rho_A<<" "<< rho_AB<<" "<<survP<<endl;
    }
    /////////////////////////////

    unsigned steps = 0;
    while(steps<T){
        /*
        if(WORKER ==0){ 
            cout<<"SAMPLE : "<<sample <<", t : "<<t<<"========================="<<100*t/(double)T<<"%\r";
        }
        */
        //Dynamics of A and B
        //cout<<"\n------------------------------------\nstep = "<<steps<<"\nA : ";
        //for(i=0;i<L;i++) cout<<A.at(i)<<" ";
        //cout<<endl;
        //cout<<"B : ";
        //for(i=0;i<L;i++) cout<<B.at(i)<<" ";
        //cout<<endl;
        Diffusion(WORKER, D_A, D_B,L, A, B,steps);
        Reaction(WORKER, L, A, B,steps);
        /*   
        cout<<"Diffusion Reaction ends \nA : ";
        for(i=0;i<L;i++) cout<<A.at(i)<<" ";
        cout<<endl;
        cout<<"B : ";
        for(i=0;i<L;i++) cout<<B.at(i)<<" ";
        cout<<endl;
        */
        double rho_A_active = 0;
        for(i = 0;i<L;i++)  if(A.at(i)!=0) rho_A_active++;
        A_sum = std::accumulate(A.begin(),A.end(),0);
        B_sum = std::accumulate(B.begin(),B.end(),0);
        rho_A_active = rho_A_active/(double)L;
        rho_A = A_sum/(double)(A_sum+B_sum);

        if(rho_A != 0){
            survP ++;
            //cout<<"\nA_sum = "<<A_sum<<endl;
        }
        else{
            //cout<<"\nReactivation activated\n";
            A_sum = std::accumulate(A.begin(),A.end(),0);
            B_sum = std::accumulate(B.begin(),B.end(),0);
            //if (B_sum == 0) break;
            Reactivation(WORKER,L, A, B,steps);
            rho_A = A_sum/(double)(A_sum+B_sum);
            for(i = 0;i<L;i++)  if(A.at(i)!=0) rho_A_active++;
            rho_A_active = rho_A_active/(double)L;
        }
        double global_rho_A = 0, global_survP = 0, global_rho2_A = 0, global_rho3_A = 0, global_rho4_A = 0;
        double global_rho_A_active = 0, global_rho2_A_active = 0, global_rho3_A_active = 0, global_rho4_A_active = 0;
        MPI::COMM_WORLD.Barrier();
        MPI::COMM_WORLD.Reduce(&rho_A,&global_rho_A,1,MPI::DOUBLE,MPI::SUM,0);
        MPI::COMM_WORLD.Reduce(&rho_A_active,&global_rho_A_active,1,MPI::DOUBLE,MPI::SUM,0);
        rho_A = rho_A*rho_A;
        rho_A_active = rho_A_active*rho_A_active;
        MPI::COMM_WORLD.Reduce(&rho_A,&global_rho2_A,1,MPI::DOUBLE,MPI::SUM,0);
        MPI::COMM_WORLD.Reduce(&rho_A_active,&global_rho2_A_active,1,MPI::DOUBLE,MPI::SUM,0);
        rho_A = rho_A*rho_A;
        rho_A_active = rho_A_active*rho_A_active;
        MPI::COMM_WORLD.Reduce(&rho_A,&global_rho3_A,1,MPI::DOUBLE,MPI::SUM,0);
        MPI::COMM_WORLD.Reduce(&rho_A_active,&global_rho3_A_active,1,MPI::DOUBLE,MPI::SUM,0);
        rho_A = rho_A*rho_A;
        rho_A_active = rho_A_active*rho_A_active;
        MPI::COMM_WORLD.Reduce(&rho_A,&global_rho4_A,1,MPI::DOUBLE,MPI::SUM,0);
        MPI::COMM_WORLD.Reduce(&rho_A_active,&global_rho4_A_active,1,MPI::DOUBLE,MPI::SUM,0);
        MPI::COMM_WORLD.Reduce(&survP,&global_survP,1,MPI::DOUBLE,MPI::SUM,0);
        if(WORKER == 0 ){
            double U = (NUM_WORKERS*global_rho2_A*global_rho3_A-global_rho_A*global_rho2_A*global_rho2_A)/(NUM_WORKERS*global_rho_A*global_rho4_A-global_rho_A*global_rho2_A*global_rho2_A);
            double delta = L*(global_rho2_A/(double)NUM_WORKERS-global_rho_A*global_rho_A/(double)(NUM_WORKERS*NUM_WORKERS));
            double U_active = (NUM_WORKERS*global_rho2_A_active*global_rho3_A_active-global_rho_A_active*global_rho2_A_active*global_rho2_A_active)/(NUM_WORKERS*global_rho_A_active*global_rho4_A_active-global_rho_A_active*global_rho2_A_active*global_rho2_A_active);
            double delta_active = L*(global_rho2_A_active/(double)NUM_WORKERS-global_rho_A_active*global_rho_A_active/(double)(NUM_WORKERS*NUM_WORKERS));
            producefile << steps<<" "<< U <<" "<<global_rho_A/(float)NUM_WORKERS<<" "<<delta<<" "<<U_active<<" "<<global_rho_A_active/(float)NUM_WORKERS<<" "<<delta_active<<endl;
        }
        survP = 0;
        steps ++; 
    } 
    //Be careful that <Rho> is averaged on the same 'steps' variable for WORKERS. The 't' would be different by the MPI WORKERS. 

    producefile.close();
    MPI::Finalize();
    return 0;
}


void Diffusion(unsigned WORKER, float D_A, float D_B,unsigned L,vector<int> &A, vector<int> &B,unsigned steps)
{
    std::mt19937 generator (200+20*WORKER+2*steps);
    std::uniform_real_distribution<double> mt_uniform(0.0,1.0);

    vector<int> A_tmp(L);
    vector<int> B_tmp(L);
    A_tmp = A; B_tmp = B;

    //Every A, B in the node can hope to its neighbors
    for(unsigned x=0;x<L;x++){
        if(A.at(x)!=0){
            unsigned A_size = A.at(x);

            for(unsigned i=0;i<A_size;i++){
                float p_move = mt_uniform(generator);
                if(D_A > p_move){
                    double p_direction = mt_uniform(generator);
                    if(p_direction>0.5){
                        //right
                        A_tmp.at(x)--;
                        unsigned right = (x + L + 1)%L;
                        A_tmp.at(right)++;
                    }
                    else{
                        //left
                        A_tmp.at(x)--;
                        unsigned left = (x + L - 1)%L;
                        A_tmp.at(left)++;
                    }
                }
            }

        }
        if(B.at(x)!=0){
            unsigned B_size = B.at(x);

            for(unsigned i=0;i<B_size;i++){
                float p_move = mt_uniform(generator);
                if(D_B > p_move){
                    double p_direction = mt_uniform(generator);
                    if(p_direction>0.5){
                        //right
                        B_tmp.at(x)--;
                        unsigned right = (x + L + 1)%L;
                        B_tmp.at(right)++;
                    }
                    else{
                        //left
                        B_tmp.at(x)--;
                        unsigned left = (x + L - 1)%L;
                        B_tmp.at(left)++;
                    }
                }
            }

        }
    }
    A = A_tmp; B = B_tmp;
}

void Reaction(unsigned WORKER,unsigned L,vector<int> &A, vector<int> &B,unsigned steps)
{
    std::mt19937 generator (10+2*WORKER+4*steps);
    std::uniform_real_distribution<double> mt_uniform(0.0,1.0);

    float mu_c = 2,mu_r=1;

    for(unsigned x =0; x<L;x++){

        if(A.at(x)!=0){
            unsigned A_size = A.at(x), B_size = B.at(x);
            double t_max = -log(1-mt_uniform(generator))/(double)(A_size+B_size), t_q = 0, Q = mu_c*A_size*B_size/(double)(A_size+B_size)+A_size*mu_r,kappa = A_size*B_size*mu_c/(double)((A_size+B_size)*Q);

            while(1){ 

                if(kappa > mt_uniform(generator)){
                    A_size++ ; B_size-- ;
                }
                else{
                    A_size-- ; B_size++ ;
                }

                Q = mu_c*A_size*B_size/(double)(A_size+B_size)+A_size*mu_r;
                kappa = A_size*B_size*mu_c/(double)((A_size+B_size)*Q);
                t_q += -log(1-mt_uniform(generator))/Q;
                //cout <<"x :"<<x<<", A.at(x) : "<<A_size<<", B.at(x) : "<<B_size<<", t_max : "<<t_max<<", Q:"<<Q<<", kappa : "<<kappa<<", t_q :"<<t_q<<endl;

                if((t_q>t_max)|(A_size==0)) break;
            }

            A.at(x) = A_size;  B.at(x) = B_size;
        }

    }

}

void Reactivation(unsigned WORKER,unsigned L, vector<int> &A, vector<int> &B,unsigned steps)
{ 
    std::mt19937 generator (1000+2*WORKER+6*steps);
    std::uniform_real_distribution<double> mt_uniform(0.0,1.0);
   
    vector<unsigned> B_arg;
    for(unsigned i = 0 ; i<L; i++){
        if(B.at(i)!=0) B_arg.push_back(i);
    }
    unsigned idx_rand = B_arg.at(B_arg.size()*(mt_uniform(generator)));
    unsigned tmp_number = B.at(idx_rand);
    A.at(idx_rand) = tmp_number;
    B.at(idx_rand) = 0;
}
