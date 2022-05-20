#pragma once

#include<random>
#include<iostream>
#include<iomanip>
#include<deque>
#include<string>
#include<fstream>
#include <Random123/philox.h>
#include <Random123/uniform.hpp>

#include"priority_queue.cpp"

//using namespace std;

typedef r123::Philox4x64 RNG;

class DEPmodel
{
	public:
		const int L, z, total_number;
		int L_A = 0, L_B = 0, L_AB = 0;
		int *nn ;
		int *particle_A, *particle_B, *particle_AB;
		std::string initial_dist; //homogeneous or single
		int seed;
		const double discrete_time = 1.0;
		//const double diffusion_A = 2/3.0, diffusion_B = 2/3.0, infection_k1 = 1, decay_k2 = 0.03278;
		//const double diffusion_A = 0.5, diffusion_B = 0.25, infection_k1 = 0.5, decay_k2 = 0.5;
		//const double diffusion_A = 1, diffusion_B = 0.5, infection_k1 = 1.0/6, decay_k2 = 5.0/6;
		//const double diffusion_A = 1, diffusion_B = 0.5, infection_k1 = 0.2, decay_k2 = 1-std::exp(-1);
		//const double diffusion_A = 2, diffusion_B = 1, infection_k1 = 0.2, decay_k2 = 1;
		const double diffusion_A = 0.5*discrete_time, diffusion_B = 0.25*discrete_time, infection_k1 = 0.5*discrete_time, decay_k2 = 0.5*discrete_time;


		DEPmodel(int L, int total_number, int z, std::string initial_dist, int seed);

		void ShowVariables();
		void Printnn();
		void PrintParticle(int * particle);
		~DEPmodel(){
			delete []nn;
			delete []particle_A;
			delete []particle_B;
			delete []particle_AB;
		}
};

class DEPmodel1d : public DEPmodel
{
	public :
		int survival_T;
		int* survival;	
		double* density_B;
		double* density_Bsite, *radius_B, *radius2_B;

		DEPmodel1d(int L , int total_number,int z, int T_max, std::string initial_dist, int seed);

		//PRE 2000 61(6), 6330 (simple algorithm)
		//void Resurect(int * particle_A, int * particle_B, int& L_A, int& L_B){
		void Resurrect(RNG::ctr_type &c, RNG::key_type &k, RNG::ctr_type & r);
		void SimpleDiffusion(int* particle, double diffusion_rate, RNG::ctr_type &c, RNG::key_type &k, RNG::ctr_type & r );
		void SimpleInfection( RNG::ctr_type &c, RNG::key_type &k, RNG::ctr_type & r);
		void SimpleHealing( RNG::ctr_type &c, RNG::key_type &k, RNG::ctr_type & r);
		void RadiusSiteDensity(int time);
		//PRE 2001 63, 066118, U.L. Fulco
		void FulcoAlgorithm(int T_max);
		void GillespieAlgorithm(int T_max);
		//void GillespieReaction();
		void NextReactionMethod(int T_max);
		void MakeDependencyGraph(int * dependency_idx);
		void UpdatePropensityTau(int* dependency_idx, double * propensity,IdxPriorityQueue& tau_q, int reaction_number,int diffusion_direction, double time, unsigned count);
		void UpdateTau(int* dependency_idx, double* old_propensity, double * propensity, IdxPriorityQueue& tau_q, int reaction_number,int diffusion_direction, double time, unsigned count);

		//void SpatioTemporalPattern(int T_max,int *alltime_A, int * alltime_B);
		void DanielMaiaAlgorithm(int T_max);

		~DEPmodel1d(){
			delete []density_B;
			delete []density_Bsite;
			delete []survival;
			delete []radius_B;
			delete []radius2_B;
		}
};


template <typename T>
double avg_sample (T * sample, int T_max);
template <typename T>
void binning (T * sample, int T_max);
