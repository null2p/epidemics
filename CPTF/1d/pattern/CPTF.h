#ifndef CPTF_H
#define CPTF_H

#include <cstdio>
#include <random>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
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
	unsigned L;
        double p_lambda;
        float theta, R_btoa;
        unsigned T;
        float D_b;
        char * seed_type;
};

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
                B.push_back(0);
            }
            else {
                //zero_arg.push_back(i);
                A.push_back(0);
                B.push_back(0);
            }
         }
    }
    else if(!strcmp(seed_type,"half")){
        for(unsigned i=0;i<L;i++){
            if(i%2==0) { // half seed
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
    }
    //rho_A = L_A/(double)L;
    //rho_AB = (L_A+L_B)/(double)L;
    //survP = 0;
}

void open_file(ofstream& Afile, ofstream& Bfile, vector<int> &A, vector<int> &B,const unsigned NUM_WORKERS,Params& params)
{
    std::stringstream Aname, Bname;
    Aname<<"./static_pattern/seed_"<<params.seed_type<<"/theta"<<params.theta<<"/Aparticle_L"<<params.L<<"T"<<params.T<<"ens"<<NUM_WORKERS<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";
    Bname<<"./static_pattern/seed_"<<params.seed_type<<"/theta"<<params.theta<<"/Bparticle_L"<<params.L<<"T"<<params.T<<"ens"<<NUM_WORKERS<<"p"<<params.p_lambda<<"theta"<<params.theta<<"R"<<params.R_btoa<<"Db"<<params.D_b<<".txt";

    Afile.open(Aname.str().c_str());  Bfile.open(Bname.str().c_str());
    Afile <<"#t    <A[x=0]>    <A[x=1]>  ...  <A[x=L-1]> ...  contact process (L="<<params.L<<", T="<<params.T<<", ENS="<<NUM_WORKERS<<", p_lambda = "<<params.p_lambda<<", theta = "<<params.theta<<", R_btoa = "<<params.R_btoa<< endl;
    Bfile <<"#t    <B[x=0]>    <B[x=1]>  ...  <B[x=L-1]> ...  contact process (L="<<params.L<<", T="<<params.T<<", ENS="<<NUM_WORKERS<<", p_lambda = "<<params.p_lambda<<", theta = "<<params.theta<<", R_btoa = "<<params.R_btoa<< endl;
    Afile << 0;  Bfile << 0; // the first time = 0;
    for(unsigned i=0; i<params.L;i++) {
        Afile <<" "<<A[i]; Bfile<< " " <<B[i];
    }
    Afile << endl;  Bfile<<endl;

}

void ABdynamics(Params& params, unsigned &L_A, unsigned &L_B, RNG::key_type &k, RNG::ctr_type &c, vector<int> &A, vector<int> &B,vector<int> &A_arg, vector<int> &B_arg, vector<double> &B_time, double &t,double& R_sq)
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
            double p = r123::u01fixedpt<double>(r.v[0]);

            unsigned index = p*(L_A+L_B);
            //if(p==1) index = L_A+L_B-1;
            //Dynamics of A
            if(index < L_A){
                double p_A = r123::u01fixedpt<double>(r.v[1]);
                //A ->B
                if(p_A < p_lambda){
                    B_arg.push_back(A_arg.at(index));

                    double tau = pow( (1-r123::u01fixedpt<double>(r.v[2]) ),-1.0/theta) + t;
                    B_time.push_back(tau);

                    A.at(A_arg.at(index))--;  B.at(A_arg.at(index))++;
                    A_arg.erase(A_arg.begin()+index);

                    L_B++;  L_A--;
                }
                else{
                    // 0A0 -> 0AA or AA0
                    // q set the direction to branch
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



#endif
