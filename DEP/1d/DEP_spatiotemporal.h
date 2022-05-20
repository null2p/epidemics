#pragma once

#include "./DEP_algorithm.h"

//#include<random>
//#include<iostream>
//#include<iomanip>
//#include<deque>
//#include<string>
//#include<fstream>

class SpatioTemporalMode : public DEPmodel
{
	public :
		int * survival;

		SpatioTemporalMode(int L , int total_number,int z, int T_max, std::string initial_dist);
		void SimpleDiffusion(int* particle, double diffusion_rate);
		void GillespieReaction();
		void SpatioTemporalPattern(int T_max,int *alltime_A, int * alltime_B);
		void NextReactionMethod(int T_max,int* alltime_A, int* alltime_B);
		void MakeDependencyGraph(int * dependency_idx);
		void UpdatePropensityTau(int* dependency_idx, double * propensity,IdxPriorityQueue& tau_q, int reaction_number , int diffusion_direction, double time);
		void UpdateTau(int* dependency_idx, double* old_propensity, double * propensity, IdxPriorityQueue& tau_q, int reaction_number,int diffusion_direction, double time);


		~SpatioTemporalMode(){
			delete []survival;
		}
};

