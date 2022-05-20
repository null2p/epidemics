#include"DEP_spatiotemporal.h"

/*
SpatioTemporalMode::SpatioTemporalMode(int L , int total_number,int z, int T_max, std::string initial_dist) : DEPmodel(L,total_number, z, initial_dist)
{
	particle_A = new int[L]{};
	particle_B = new int[L]{};
	particle_AB = new int[L]{};
	//nn = new int[z*L]{};

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> rng;
	L_A = total_number - 10;  L_B = 10;
	//L_A = total_number/2;  L_B = total_number - L_A;
	//L_A = 0; L_B = total_number;

	if(initial_dist == "homogeneous"){
		std::uniform_int_distribution<> rng_L(0,L-1);
		for(int i =0; i<L_A;i++){
			int idx = rng_L(mt);
			//if (idx%2 ==0)  idx = nn[2*idx];
			particle_A[idx] ++;
		}
		for(int i =0; i<L_B;i++){
			int idx = rng_L(mt);
			//if (idx%2 !=0)  idx = nn[2*idx];
			particle_B[idx] ++;
		}
	}
	else if(initial_dist == "single"){
	
		L_A = total_number - 10;  L_B = 10;

		for(int idx = 0 ; idx < L; idx++){
		
			particle_A[idx] = total_number/L;
			//if(idx == L/2) particle_A[idx] = 0;

		}
		
		int residue = total_number - L*(total_number/L);

		std::uniform_int_distribution<> rng_L(0,L-1);
		for(int i =0; i<residue;i++) particle_A[rng_L(mt)] ++;

		particle_A[L/2] -= L_B;
		if(particle_A[L/2]<0){
			residue = -particle_A[L/2];
			particle_A[L/2] = 0;
			for(int i =0; i<residue; ){
				int idx = rng_L(mt);
				if(idx == L/2)  continue;
				particle_A[idx] --;
				i++;
			}
		}
		particle_B[L/2] = L_B;
	}
	for(int i =0; i<L;i++)  particle_AB[i] = particle_A[i]*particle_B[i];

}
*/

SpatioTemporalMode::SpatioTemporalMode(int L , int total_number,int z, int T_max, std::string initial_dist) : DEPmodel(L,total_number, z, initial_dist)
{
        particle_A = new int[L]{};
        particle_B = new int[L]{};
        particle_AB = new int[L]{};
        //nn = new int[z*L]{};

        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> rng;

        if(initial_dist == "homogeneous"){
                L_A = 0; L_B = total_number;

                for (int idx =0; idx<L ; idx++) particle_B[idx] = L_B / L;

                int residue = total_number - L*(L_B/L);
                std::uniform_int_distribution<> rng_L(0,L-1);

                for (int i =0; i<residue ; i++){
                        int idx = rng_L(mt);
                        particle_B[idx] ++;
                }
        }
        else if(initial_dist == "single"){
                L_A = total_number - 10 ; L_B = 10;

                for (int idx = 0 ; idx < L ; idx ++) particle_A[idx] = total_number / L;

                int residue = total_number - L*(total_number/L);
                std::uniform_int_distribution<> rng_L(0,L-1);

                for (int i =0; i<residue ; i++){
                        int idx = rng_L(mt);
                        particle_A[idx] ++;
                }

                particle_B[L/2] = L_B;
                residue = particle_B[L/2] - particle_A[L/2];
                particle_A[L/2] = 0;

                for(int i =0; i<residue; i++){
                        int idx = -1;
                        while(1){
                                idx = rng_L(mt);
                                if(particle_A[idx]!=0) break;
                        }
                        particle_A[idx] --;
                }
                if(residue < 0){
                        for(int i =0; i< -residue; i++){
                                int idx = -1;
                                while(1){
                                        idx = rng_L(mt);
                                        if(particle_A[idx]!=0) break;
                                }
                                particle_A[idx] ++;
                        }
                }
        }
        //for(int i =0; i<L_B;i++)  particle_B[rng_L(mt)] ++;
        for(int i =0; i<L;i++)  particle_AB[i] = particle_A[i]*particle_B[i];

}


void SpatioTemporalMode::SimpleDiffusion(int* particle, double diffusion_rate){
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> rng;
	int *tmp_particle = new int[L]{};
	//O(particle_number), each particle is asked to diffuse or not
	for(int i=0;i<L;i++){
		if(particle[i] != 0){
			for(int j=0;j<particle[i];j++){
				double rng_diff = rng(mt);
				if(diffusion_rate/2 > rng_diff){
					tmp_particle[i] --; tmp_particle[nn[2*i]]++;
				}
				else if(diffusion_rate > rng_diff){
					tmp_particle[i] --; tmp_particle[nn[2*i+1]]++;
				}
			}
			// fixed number by the diff rate is used. not a random number. caution : small number of particle does not diffuse at all.
			//int diffuseNumber = particle[i]*diffusion_rate;
			//particle[nn[2*i]] += diffuseNumber; particle[nn[2*i+1]] += diffuseNumber; particle[i] -= 2*diffuseNumber;
		}
	}

	for(int i =0;i<L;i++)  particle[i] += tmp_particle[i];
	delete []tmp_particle;

}




void SpatioTemporalMode::GillespieReaction(){
	
	std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> rng(0.0,1.0);


	for(int i=0;i<L;i++){
		if(particle_B[i] > 0){
			double t_max = -std::log(1-rng(mt))/(particle_A[i]+particle_B[i]);
			int tmp_B = 0, tmp_A = 0;

			while(t_max>0){
				double infection_p = infection_k1*particle_B[i]*particle_A[i];
				double R_gill = infection_p + decay_k2*particle_B[i];
				infection_p /= R_gill;
				double t_update = -std::log(1-rng(mt))/R_gill;

				if(infection_p > rng(mt)){
					tmp_B++; tmp_A--;
				}
				else{
					tmp_B--; tmp_A++;
				}

				t_max -= t_update;
				if(particle_B[i]==-tmp_B) break;
			}
			particle_B[i] += tmp_B; particle_A[i] += tmp_A;
			L_B += tmp_B; L_A += tmp_A;
		}
	}

}

void SpatioTemporalMode::SpatioTemporalPattern(int T_max,int *alltime_A, int * alltime_B){



	NextReactionMethod(T_max,alltime_A, alltime_B);

}

void SpatioTemporalMode::NextReactionMethod(int T_max, int *alltime_A, int*alltime_B){

        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> rng(0.0,1.0);

	survival = new int[T_max]{};

        int* dependency_idx;
        dependency_idx = new int[4*L*9]{};
        MakeDependencyGraph(dependency_idx);

        double *propensity, *tau;
        propensity = new double[4*L]{};
        tau = new double[4*L]{};

        for(int i=0;i<L;i++){

                propensity[4*i] = infection_k1*particle_A[i]*particle_B[i]; //infection
                propensity[4*i+1] = decay_k2*particle_B[i];                 //decay
                propensity[4*i+2] = diffusion_A*particle_A[i];              //diffusion of A
                propensity[4*i+3] = diffusion_B*particle_B[i];              //diffusion of B

                if(propensity[4*i] == 0)  tau[4*i] = 2147483647;
                else tau[4*i] = -std::log(1-rng(mt))/propensity[4*i];

                if(propensity[4*i+1] == 0)  tau[4*i+1] = 2147483647;
                else tau[4*i+1] = -std::log(1-rng(mt))/propensity[4*i+1];

                if(propensity[4*i+2] == 0)  tau[4*i+2] = 2147483647;
                else tau[4*i+2] = -std::log(1-rng(mt))/propensity[4*i+2];

                if(propensity[4*i+3] == 0)  tau[4*i+3] = 2147483647;
                else tau[4*i+3] = -std::log(1-rng(mt))/propensity[4*i+3];

        }

        IdxPriorityQueue tau_q(4*L,tau);
        delete []tau;


	int measure_time = 1;
	int measure_count = 0;
        for(double t=0;t<T_max;){

                int reaction_number = tau_q.heap_element[0].reaction;
                double tau = tau_q.heap_element[0].tau;

                int channel = reaction_number%4;
                int node = reaction_number/4;
                int diffusion_direction = -1;
                if(channel == 0){
                        particle_A[node]--;
                        particle_B[node]++;
                        L_A--; L_B++;
                }
                else if(channel == 1){
                        particle_A[node] ++;
                        particle_B[node] --;
                        L_A++; L_B--;
                }
                else if(channel == 2){
                        particle_A[node]--;
                        if(rng(mt)>0.5){
                                particle_A[nn[2*node]]++;
                                diffusion_direction = 0;
                        }
                        else{
                                particle_A[nn[2*node+1]]++;
                                diffusion_direction = 1;
                        }
                }
                else if(channel == 3){
                        particle_B[node]--;
                        if(rng(mt)>0.5){
                                particle_B[nn[2*node]]++;
                                diffusion_direction = 2;
                                //A
                                //left_B=1;
                        }
                        else{
                                particle_B[nn[2*node+1]]++;
                                diffusion_direction = 3;
                                //right_B=1;
                        }
                }

                t = tau;

                if(L_B == 0) break;

                UpdatePropensityTau(dependency_idx, propensity,tau_q, reaction_number,diffusion_direction,t);

                //UpdateTau(dependency_idx, propensity, tau_q, reaction_number, t);


                if(t >= measure_time){
                        survival[measure_count] = 1;
			
			for(int i =0; i<L;i++){
		       		alltime_A[L*measure_count + i] += particle_A[i];
		       		alltime_B[L*measure_count + i] += particle_B[i];
			}
                        measure_time ++;
                        measure_count++;
                }
//std::cout<<t<<"\t"<<(double)L_B/L<<"\n";

        }

        delete []dependency_idx;
        delete []propensity;

}

void SpatioTemporalMode::MakeDependencyGraph(int * dependency_idx){


        for(int i = 0 ; i < 4*L ; i++ ){

                int channel = i%4;
                int node = i/4;

                if(channel == 0){
                        //INFECTION
                        dependency_idx[9*i + 0] = i;   //infection
                        dependency_idx[9*i + 1] = i+1; //decay
                        dependency_idx[9*i + 2] = i+2; //diffusion_A
                        dependency_idx[9*i + 3] = i+3; //diffusion_B
                        dependency_idx[9*i + 4] = -1;  //nothing
                        dependency_idx[9*i + 5] = -1;  //nothing
                }
                else if(channel == 1){
                        //DECAY
                        dependency_idx[9*i + 0] = i;   //decay
                        dependency_idx[9*i + 1] = i-1; //infection
                        dependency_idx[9*i + 2] = i+1; //diffusion_A
                        dependency_idx[9*i + 3] = i+2; //diffusion_B
                        dependency_idx[9*i + 4] = -1;  //nothing
                        dependency_idx[9*i + 5] = -1;  //nothing

                }
                else if(channel == 2){
                        //DIFFUSION A
                        dependency_idx[9*i + 0] = i;                      //diffusion_A
                        dependency_idx[9*i + 1] = i-2;                    //infection
                        dependency_idx[9*i + 2] = 4*nn[2*node] + 0;       //infection left
                        dependency_idx[9*i + 3] = 4*nn[2*node] + channel; //diffusion_A, left
                        dependency_idx[9*i + 4] = 4*nn[2*node+1] + 0;     //infection right
                        dependency_idx[9*i + 5] = 4*nn[2*node+1]+channel; //diffusion_A, right
                        dependency_idx[9*i + 6] = -1;  //nothing

                }else if(channel == 3){
                        //DIFFUSION B
                        dependency_idx[9*i + 0] = i;                      //diffusion_B
                        dependency_idx[9*i + 1] = i-3;                    //infection
                        dependency_idx[9*i + 2] = i-2;                    //decay
                        dependency_idx[9*i + 3] = 4*nn[2*node] + 0;       //infection left
                        dependency_idx[9*i + 4] = 4*nn[2*node] + 1;       //decay left
                        dependency_idx[9*i + 5] = 4*nn[2*node]+channel;   //diffusion_B, left
                        dependency_idx[9*i + 6] = 4*nn[2*node+1] + 0;     //infection right
                        dependency_idx[9*i + 7] = 4*nn[2*node+1] + 1;     //decay right
                        dependency_idx[9*i + 8] = 4*nn[2*node+1]+channel; //diffusion_B, right
                }


        }

}


void SpatioTemporalMode::UpdatePropensityTau(int* dependency_idx, double * propensity,IdxPriorityQueue& tau_q, int reaction_number , int diffusion_direction, double time){


        int channel = reaction_number%4;
        int node = reaction_number/4;

        double * old_propensity;
        old_propensity = new double[9]{};

        if(channel == 0){
                for(int i=0;i<4;i++) old_propensity[i] = propensity[dependency_idx[9*reaction_number+i]];

                propensity[dependency_idx[9*reaction_number+0]] = particle_A[node]*particle_B[node]*infection_k1;
                propensity[dependency_idx[9*reaction_number+1]] = particle_B[node]*decay_k2;
                propensity[dependency_idx[9*reaction_number+2]] = particle_A[node]*diffusion_A;
                propensity[dependency_idx[9*reaction_number+3]] = particle_B[node]*diffusion_B;

                UpdateTau(dependency_idx,old_propensity, propensity, tau_q, reaction_number,diffusion_direction, time);
        }
        else if(channel == 1){
                for(int i=0;i<4;i++) old_propensity[i] = propensity[dependency_idx[9*reaction_number+i]];

                propensity[dependency_idx[9*reaction_number+0]] = particle_B[node]*decay_k2;
                propensity[dependency_idx[9*reaction_number+1]] = particle_A[node]*particle_B[node]*infection_k1;
                propensity[dependency_idx[9*reaction_number+2]] = particle_A[node]*diffusion_A;
                propensity[dependency_idx[9*reaction_number+3]] = particle_B[node]*diffusion_B;

                UpdateTau(dependency_idx,old_propensity, propensity, tau_q, reaction_number,diffusion_direction, time);
        }
        else if(channel == 2){
                for(int i=0;i<6;i++) old_propensity[i] = propensity[dependency_idx[9*reaction_number+i]];

                propensity[dependency_idx[9*reaction_number+0]] = particle_A[node]*diffusion_A;
                propensity[dependency_idx[9*reaction_number+1]] = particle_A[node]*particle_B[node]*infection_k1;

                if(diffusion_direction == 0){
                        propensity[dependency_idx[9*reaction_number+2]] = particle_A[nn[2*node]]*particle_B[nn[2*node]]*infection_k1;
                        propensity[dependency_idx[9*reaction_number+3]] = particle_A[nn[2*node]]*diffusion_A;
                }
                else if(diffusion_direction == 1){
                        propensity[dependency_idx[9*reaction_number+4]] = particle_A[nn[2*node+1]]*particle_B[nn[2*node+1]]*infection_k1;
                        propensity[dependency_idx[9*reaction_number+5]] = particle_A[nn[2*node+1]]*diffusion_A;
                }
                UpdateTau(dependency_idx,old_propensity, propensity, tau_q, reaction_number,diffusion_direction, time);
        }
        else if(channel == 3){
                for(int i=0;i<9;i++) old_propensity[i] = propensity[dependency_idx[9*reaction_number+i]];

                propensity[dependency_idx[9*reaction_number+0]] = particle_B[node]*diffusion_B;
                propensity[dependency_idx[9*reaction_number+1]] = particle_A[node]*particle_B[node]*infection_k1;
                propensity[dependency_idx[9*reaction_number+2]] = particle_B[node]*decay_k2;

                if(diffusion_direction ==2){
                        propensity[dependency_idx[9*reaction_number+3]] = particle_A[nn[2*node]]*particle_B[nn[2*node]]*infection_k1;
                        propensity[dependency_idx[9*reaction_number+4]] = particle_B[nn[2*node]]*decay_k2;
                        propensity[dependency_idx[9*reaction_number+5]] = particle_B[nn[2*node]]*diffusion_B;
                }
                else if (diffusion_direction ==3){
                        propensity[dependency_idx[9*reaction_number+6]] = particle_A[nn[2*node+1]]*particle_B[nn[2*node+1]]*infection_k1;
                        propensity[dependency_idx[9*reaction_number+7]] = particle_B[nn[2*node+1]]*decay_k2;
                        propensity[dependency_idx[9*reaction_number+8]] = particle_B[nn[2*node+1]]*diffusion_B;
                }
                UpdateTau(dependency_idx,old_propensity, propensity, tau_q, reaction_number,diffusion_direction, time);
        }

        delete []old_propensity;

}

void SpatioTemporalMode::UpdateTau(int* dependency_idx, double* old_propensity, double * propensity, IdxPriorityQueue& tau_q, int reaction_number,int diffusion_direction, double time){

        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> rng(0.0,1.0);

        int channel = reaction_number%4;

        if(channel < 2){
                for(int i = 0; i < 9; i++){

                        int new_reaction_number = dependency_idx[9*reaction_number+i];
                        //since No Reaction is ... dependency_idx[] = -1;
                        if(new_reaction_number == -1) break;

                        int heap_idx = tau_q.heap_index[new_reaction_number];
                        double new_tau = 2147483647;

                        if(propensity[new_reaction_number]!=0){
                                if(i==0)  new_tau = time - std::log(1-rng(mt))/propensity[new_reaction_number];
                                else{
                                        if(old_propensity[i]!=0) new_tau = time + old_propensity[i]*(tau_q.heap_element[heap_idx].tau - time)/propensity[new_reaction_number];
                                        else new_tau = time - std::log(1-rng(mt))/propensity[new_reaction_number];

                                }
                        }
                        tau_q.update(heap_idx, new_tau);
                }
        }
        else if(channel == 2){
                int * diff_index;
                if(diffusion_direction == 0) diff_index = new int[4]{0,1,2,3};
                else if(diffusion_direction == 1) diff_index = new int[4]{0,1,4,5};
                for(int i = 0; i < 4; i++){
                        int new_reaction_number = dependency_idx[9*reaction_number+diff_index[i]];
                        //since No Reaction is ... dependency_idx[] = -1;

                        int heap_idx = tau_q.heap_index[new_reaction_number];
                        double new_tau = 2147483647;

                        if(propensity[new_reaction_number]!=0){
                                if(i==0)  new_tau = time - std::log(1-rng(mt))/propensity[new_reaction_number];
                                else{
                                        if(old_propensity[diff_index[i]]!=0) new_tau = time + old_propensity[diff_index[i]]*(tau_q.heap_element[heap_idx].tau - time)/propensity[new_reaction_number];
                                        else new_tau = time - std::log(1-rng(mt))/propensity[new_reaction_number];

                                }
                        }
                        tau_q.update(heap_idx, new_tau);
                }
                delete []diff_index;

        }
        else if(channel == 3){
                int * diff_index;
                if(diffusion_direction == 2) diff_index = new int[6]{0,1,2,3,4,5};
                else if(diffusion_direction == 3) diff_index = new int[6]{0,1,2,6,7,8};
                for(int i = 0; i < 6; i++){
                        int new_reaction_number = dependency_idx[9*reaction_number+diff_index[i]];
                        //since No Reaction is ... dependency_idx[] = -1;

                        int heap_idx = tau_q.heap_index[new_reaction_number];
                        double new_tau = 2147483647;

                        if(propensity[new_reaction_number]!=0){
                                if(i==0)  new_tau = time - std::log(1-rng(mt))/propensity[new_reaction_number];
                                else{
                                        if(old_propensity[diff_index[i]]!=0) new_tau = time + old_propensity[diff_index[i]]*(tau_q.heap_element[heap_idx].tau - time)/propensity[new_reaction_number];
                                        else new_tau = time - std::log(1-rng(mt))/propensity[new_reaction_number];

                                }
                        }
                        tau_q.update(heap_idx, new_tau);
                }
                delete []diff_index;

        }


}
