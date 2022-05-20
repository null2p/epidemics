#include"DEP_algorithm.h"

//typedef r123::Philox4x64 RNG;

DEPmodel::DEPmodel(int L, int total_number, int z, std::string initial_dist, int seed) : L(L), total_number(total_number), z(z), initial_dist(initial_dist), seed(seed)
//DEPmodel::DEPmodel(int L, int total_number, int z)
{
	//for(int i =0;i<L*L;i++)  particle_A[i] = ( rng(mt) > 0.5? 1:-1 );
	//for(int i =0;i<L*L;i++)  particle_A[i] = 1;
	if(z==2){
		nn = new int[2*L]{};
		for(int i=0; i<L; i++){
			int left = (i-1+L)%L, right = (i+1)%L;
			nn[z*i] = left; nn[z*i+1] = right;
		}
	}
	else if(z==4){
		int N = L*L;
		nn = new int[4*L*L]{};
		for(int i = 0; i<N;i++){
			int left = (i-1)%L + (i/L)*L, right = (i+1)%L + (i/L)*L, down=(i+L)%N, up=(i-L)%N;
			if(up<0) up+=N;
			nn[z*i] = left; nn[z*i+1] = right; nn[z*i+2] = up; nn[z*i+3] = down;
		}
		nn[0] = L-1;
	}
	else if(z==6){
		int N = L*L;
		nn = new int[6*L*L]{};
		for(int i = 0; i<N;i++){
			int left = (i-1)%L + (i/L)*L, right = (i+1)%L + (i/L)*L, down=(i+L)%N, up=(i-L)%N, rup = (right -L)%N, ldown = (left+L)%N;
			nn[z*i] = left; nn[z*i+1] = right; nn[z*i+2] = up; nn[z*i+3] = down; nn[z*i+4] = rup;
		}
	}
}

void DEPmodel::ShowVariables(){
	std::cout<<"\n ------------\n L = "<<L<<"\n z = "<<z<<"\n Total number = "<<total_number<<"\n ------------\n";
}

void DEPmodel::Printnn(){
	if(z==2){
		for(int i =0 ; i < 2*L; i++){
			std::cout<<nn[i]<<" ";
		}
		std::cout<<std::endl;
	}
	else if(z==4){
		for(int i =0 ; i < 4*L*L; i++){
			std::cout<<nn[i]<<" ";
			if(i%(4*L) == 4*L-1) std::cout<<std::endl;
		}
	}
}

void DEPmodel::PrintParticle(int * particle){
	int total = 0;
	for(int i =0 ; i < L; i++){
		std::cout<<std::setw(2)<<particle[i]<<" ";
		total += particle[i];
		//if(particle[i]==1)  std::cout<<setw(2)<<"#"<<" ";
		//if(particle[i]==-1)  std::cout<<setw(2)<<"0"<<" ";
		//if(i%L == L-1) std::cout<<endl;
	}
	std::cout<<"Total number : "<<total<<std::endl;
	std::cout<<"Density : "<<(double)total/total_number<<std::endl;
}



DEPmodel1d::DEPmodel1d(int L , int total_number,int z, int T_max, std::string initial_dist, int seed) : DEPmodel(L,total_number, z, initial_dist, seed)
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
/*
void DEPmodel1d::Resurrect(RNG::ctr_type &c, RNG::key_type &k, RNG::ctr_type & r){

	RNG rng;
        c[0]++;  r = rng(c,k);

	std::deque<int> dq_candidate ;
	for(int i =0; i<L; i++){
		if(particle_A[i] > 0)  dq_candidate.push_back(i);
	}	

	int node_resurect = dq_candidate[dq_candidate.size()*r123::u01fixedpt<double>(r.v[0])];
	L_B = 1; L_A -= 1;
	particle_B[node_resurect] = 1; particle_A[node_resurect] -= 1;

	//L_B = particle_A[node_resurect]; L_A -= particle_A[node_resurect];
	//particle_B[node_resurect] = L_B; particle_A[node_resurect] = 0;

}
*/


void DEPmodel1d::Resurrect(RNG::ctr_type &c, RNG::key_type &k, RNG::ctr_type & r){

	RNG rng;
        c[0]++;  r = rng(c,k);

	int node_resurrect = -1;

	for(int i = 0; ;i++){
		node_resurrect = L*r123::u01fixedpt<double>(r.v[i%4]);
		if(particle_A[i] != 0) break;
	}

	L_B = 1; L_A -=1;
	particle_B[node_resurrect] = 1; particle_A[node_resurrect] -= 1;

}

/*
void DEPmodel1d::SimpleDiffusion(int* particle, double diffusion_rate,int t){
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> rng;
	int *tmp_particle = new int[L]{};
	int min_i = L/2, max_i = L/2;
	//O(particle_number), each particle is asked to diffuse or not
	for(int i=0;i<L;i++){
		if(particle[i] != 0){
			if(i < min_i) min_i = i;
			if(i > max_i) max_i = i;
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

	radius_B[t] = ( (L/2 - min_i)*(L/2 - min_i) + (max_i - L/2)*(max_i - L/2) )/2 ;

	for(int i =0;i<L;i++)  particle[i] += tmp_particle[i];
	delete []tmp_particle;

}
*/

void DEPmodel1d::SimpleDiffusion(int* particle, double diffusion_rate, RNG::ctr_type &c, RNG::key_type &k, RNG::ctr_type & r){

	RNG rng;
	//r = rng(c,k);

	int *tmp_particle = new int[L]{};
	//O(particle_number), each particle is asked to diffuse or not
	for(int i=0;i<L;i++){
		if(particle[i] != 0){
			for(int j=0;j<particle[i];j++){
	
				if(j%4==0){
					c[0]++;  r = rng(c,k);
				}

				double rng_diff = r123::u01fixedpt<double>(r.v[j%4]);
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


void DEPmodel1d::SimpleInfection( RNG::ctr_type &c, RNG::key_type &k, RNG::ctr_type & r){

	RNG rng;
	//r =rng(c,k);
	
	for(int i = 0; i<L;i++){
		if(particle_B[i]!=0){
//				if(particle_A[i]!=0 & (infection_k1 > rng(mt))){
//					L_A -= particle_A[i]; L_B += particle_A[i];
//					particle_B[i] +=particle_A[i];
//					particle_A[i] = 0;
//				}
			//fulco, each A has same infection rate k1 with at least one B particle.
//				if(particle_A[i]!=0){
//					int tmp_infection = 0;
//					for(int j=0;j<particle_A[i];j++){
//						if( infection_k1 > rng(mt)) tmp_infection++;
//					}
//					particle_B[i] += tmp_infection; particle_A[i] -= tmp_infection;
//					L_B += tmp_infection; L_A -= tmp_infection;
//					
//
//				}
			//double no_infection_rate = std::pow(1-infection_k1,particle_B[i]);
			if(particle_A[i]!=0){
				int tmp_infection = 0;
				for(int j=0; j< particle_A[i] ; j++){
					if(j%4==0){
						c[0]++; r = rng(c,k);
					}
					if( infection_k1 > r123::u01fixedpt<double>(r.v[j%4]) ) tmp_infection++;
					//if( infection_k1 > rng(mt)) infection_tmp[i]++;
					//if(1-no_infection_rate > rng(mt))  tmp_infection ++;
				}
				particle_B[i] += tmp_infection; particle_A[i] -= tmp_infection;
				L_B += tmp_infection; L_A -= tmp_infection;
			}

		}
	}


}

void DEPmodel1d::SimpleHealing( RNG::ctr_type &c, RNG::key_type &k, RNG::ctr_type & r){

	RNG rng;
	//r = rng(c,k);

	for(int i=0;i<L;i++){

		if(particle_B[i]!=0){
			int tmp_decay = 0;
			for(int j=0;j<particle_B[i];j++){
				if(j%4==0){
					c[0]++; r = rng(c,k);
				}

				if(decay_k2 > r123::u01fixedpt<double>(r.v[j%4]) ){
					tmp_decay ++;
				}
			}
			particle_B[i] -= tmp_decay; particle_A[i] += tmp_decay;
			L_B -= tmp_decay; L_A += tmp_decay;
		}
	}
}

void DEPmodel1d::RadiusSiteDensity(int time){

	int min_idx = L-1, max_idx = 0, tmp_radius =0;
	//int tmp_radius =0;
	for(int i =0; i<L;i++){
		if(particle_B[i]!=0){
		
			density_Bsite[time] ++;
			radius_B[time] += (i-L/2)*(i-L/2);

			if(i < min_idx) min_idx = i;
			if(i > max_idx) max_idx = i;

		}
	}
	if (density_Bsite[time] == 0)  return ;

	radius_B[time] /= density_Bsite[time];
	radius2_B[time] = (max_idx - min_idx)*(max_idx-min_idx);

	density_Bsite[time] /= L;

}


//void DEPmodel1d::RadiusSiteDensity(int time){
//
//	int min_idx = L-1, max_idx = 0, tmp_radius =0;
//	for(int i =0; i<L;i++){
//		if(particle_B[i]!=0){
//			if(tmp_radius == 0) tmp_radius ++;
//
//			if(particle_B[nn[2*i+1]]!=0) tmp_radius ++;
//			else if (particle_B[nn[2*i+1]]==0){
//				radius_B[time] += tmp_radius*tmp_radius;
//				if(tmp_radius != 0) N_Bcluster[time]++;
//				
//				tmp_radius = 0;
//			}
//
//
//			density_Bsite[time] ++;
//			//if(i < min_idx) min_idx = i;
//			//if(i > max_idx) max_idx = i;
//		}
//	}
//	if (density_Bsite[time] == 0)  return ;
//
//	//radius_B[time] = (max_idx-min_idx)*(max_idx-min_idx);
//	if(N_Bcluster[time]!=0)  N_BinCluster[time] = L_B/(double)N_Bcluster[time];
//	density_Bsite[time] /= L;
//
//}

//PRE 2001 63, 066118, U.L. Fulco
void DEPmodel1d::FulcoAlgorithm(int T_max){


	int measure_T_max = ((int)std::log(T_max)/std::log(10)-1)*30 + 1 + 10;
	//int measure_T_max = T_max;

	density_B = new double[measure_T_max]{};
	density_Bsite = new double[measure_T_max]{};
        survival = new int[measure_T_max]{};
	radius_B = new double[measure_T_max]{};
	radius2_B = new double[measure_T_max]{};

	RNG rng;
        RNG::ctr_type c={{1,0,0,0}};
        RNG::ukey_type uk={{1,2}};
        uk[0] = seed; // some user_supplied_seed
        RNG::key_type k=uk;
        RNG::ctr_type r;

	survival_T = T_max;
	double time_exp = 1.0 + 1/30.0;
	int measure_time = 1;
	int measure_count = 0;
	unsigned count = 0 ;



	for(double t=0;t<T_max;){
	//for(int t=0;t<T_max;t++){

		//diffusion of A and B
		SimpleDiffusion(particle_A,diffusion_A, c, k, r);
		SimpleDiffusion(particle_B,diffusion_B, c, k, r);
		SimpleInfection(c,k,r);
		SimpleHealing(c,k,r);

		if (L_B == 0)  break;

		if(t >= measure_time){

			//std::cout<<"t : "<<t<<", measure_time : "<<measure_time<<", measure_count : "<<measure_count<<std::endl;

			//RadiusSiteDensity(measure_count);
			density_B[measure_count] = (double)L_B/L;
                	survival[measure_count] = 1;
			//measure_time ++;
			if(t<10)  measure_time ++;
			else{
				time_exp += 1/30.0;
				measure_time = std::pow(10.0,time_exp);	
			}	
			measure_count++;

		}

		t= t+discrete_time;
	

		/*	
		if(L_B == 0){
			survival_T = t;
			Resurrect(c,k,r);
			//break;
		}


		density_B[t] = (double)L_B/L;
		*/
	}
}

void DEPmodel1d::GillespieAlgorithm(int T_max){

        density_B = new double[T_max]{};
        survival = new int[T_max]{};
	density_Bsite = new double[T_max]{};
	radius_B = new double[T_max]{};

	RNG rng;
        RNG::ctr_type c={{1,0,0,0}};
        RNG::ukey_type uk={{1,2}};
        uk[0] = seed; // some user_supplied_seed
        RNG::key_type k=uk;
        RNG::ctr_type r;

				

        for(int t=0;t<T_max;t++){

                //diffusion of A and B
                SimpleDiffusion(particle_A,diffusion_A,c,k,r);
                SimpleDiffusion(particle_B,diffusion_B,c,k,r);

                for(int i=0;i<L;i++){
                        if(particle_B[i] > 0){
				if(i%4==0){
					c[0]++; r = rng(c,k);
				}

                                double t_max = -std::log(r123::u01fixedpt<double>(r.v[i%4]))/(particle_A[i]+particle_B[i]);
                                int tmp_B = 0, tmp_A = 0, steps = 0;

                                while(t_max>0){
					c[0]++; r = rng(c,k);

                                        double R_gill = infection_k1*particle_B[i]*particle_A[i] + decay_k2*particle_B[i];
                                        double infection_p = infection_k1*particle_B[i]*particle_A[i];
                                        infection_p /= R_gill;
                                        double t_update = -std::log(r123::u01fixedpt<double>(r.v[0]) )/R_gill;

                                        if(infection_p >r123::u01fixedpt<double>(r.v[1]) ){
                                                tmp_B ++; tmp_A--;
                                        }
                                        else{
                                                tmp_B--; tmp_A++;
                                        }

                                        t_max -= t_update;
                                        if(particle_B[i]==-tmp_B) break;

					steps ++;
                                }
                                particle_B[i] += tmp_B; particle_A[i] += tmp_A;
                                L_B += tmp_B; L_A += tmp_A;
                        }
                }
		

                if(L_B == 0) break;
		
		//RadiusSiteDensity(t);
                density_B[t] = (double)L_B/L;
                survival[t] = 1;

        }

}

void DEPmodel1d::NextReactionMethod(int T_max){

	//std::random_device rd;
        //std::mt19937 mt(rd());
        //std::uniform_real_distribution<double> rng(0.0,1.0);

	RNG rng;
        RNG::ctr_type c={{1,0,0,0}};
        RNG::ukey_type uk={{1,2}};
        uk[0] = seed; // some user_supplied_seed
        RNG::key_type k=uk;
        RNG::ctr_type r;

	//int measure_T_max = 10000;
	//int measure_T_max = ((int)std::log10(T_max))*10 - 10+1;
	int measure_T_max = ((int)std::log10(T_max)-1)*30 + 1 + 10;

        density_B = new double[measure_T_max]{};
        survival = new int[measure_T_max]{};
	density_Bsite = new double[measure_T_max]{};
	radius_B = new double[measure_T_max]{};
	radius2_B = new double[measure_T_max]{};

	int* dependency_idx;
	dependency_idx = new int[4*L*9]{};
	MakeDependencyGraph(dependency_idx);

	double *propensity, *tau;
	propensity = new double[4*L]{};
	tau = new double[4*L]{};

	for(int i=0;i<L;i++){

		c[0] = i; c[1] = i+L;
		r = rng(c,k);

		propensity[4*i] = infection_k1*particle_A[i]*particle_B[i]; //infection
		propensity[4*i+1] = decay_k2*particle_B[i];		    //decay
		propensity[4*i+2] = diffusion_A*particle_A[i];		    //diffusion of A
		propensity[4*i+3] = diffusion_B*particle_B[i];		    //diffusion of B

		if(propensity[4*i] == 0)  tau[4*i] = 2147483647;
		else tau[4*i] = -std::log(r123::u01fixedpt<double>(r.v[0]))/propensity[4*i];

		if(propensity[4*i+1] == 0)  tau[4*i+1] = 2147483647;
		else tau[4*i+1] = -std::log(r123::u01fixedpt<double>(r.v[1]))/propensity[4*i+1];

		if(propensity[4*i+2] == 0)  tau[4*i+2] = 2147483647;
		else tau[4*i+2] = -std::log(r123::u01fixedpt<double>(r.v[2]))/propensity[4*i+2];

		if(propensity[4*i+3] == 0)  tau[4*i+3] = 2147483647;
		else tau[4*i+3] = -std::log(r123::u01fixedpt<double>(r.v[3]))/propensity[4*i+3];

	}

	IdxPriorityQueue tau_q(4*L,tau);
	delete []tau;

	double time_exp = 1.0 + 1/30.0;
	int measure_time = 1;
	int measure_count = 0;
	unsigned count = 0 ;

        for(double t=0;t<T_max;){

		int reaction_number = tau_q.heap_element[0].reaction;
		double tau = tau_q.heap_element[0].tau;

		int channel = reaction_number%4;
		int node = reaction_number/4;
		int diffusion_direction = -1;
//cout<<"\nreaction # : "<<reaction_number<<", tau : "<<tau<<"\n";
//cout<<"\nchannel : "<<channel<<", node : "<<node<<", A[node] : "<<particle_A[node]<<", B[node] : "<<particle_B[node]<<"\n";
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

			++c[0];
			r = rng(c,k);

			particle_A[node]--;
			if(r123::u01fixedpt<double>(r.v[0])>0.5){
				particle_A[nn[2*node]]++;
				diffusion_direction = 0;
			}
			else{
			       	particle_A[nn[2*node+1]]++;
				diffusion_direction = 1;
			}
		}
		else if(channel == 3){

			++c[1];
			r = rng(c,k);

			particle_B[node]--;
			if(r123::u01fixedpt<double>(r.v[1])>0.5){
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
		
//cout<<"\nAFTER REACTION ... A[node] : "<<particle_A[node]<<", B[node] : "<<particle_B[node]<<"\n";
		t = tau;
                
                if(L_B == 0) Resurrect(c,k,r);
                //if(L_B == 0) break;

		UpdatePropensityTau(dependency_idx, propensity,tau_q, reaction_number,diffusion_direction,t,count);

		//UpdateTau(dependency_idx, propensity, tau_q, reaction_number, t);

		if(t >= measure_time){


			//std::cout<<"t : "<<t<<", measure_time : "<<measure_time<<std::endl;

			//RadiusSiteDensity(measure_count);
			density_B[measure_count] = (double)L_B/L;
                	survival[measure_count] = 1;
			//measure_time ++;
			if(t<10)  measure_time ++;
			else{
				time_exp += 1/30.0;
				measure_time = std::pow(10.0,time_exp);	
			}	
			measure_count++;

		}
		count ++;
//std::cout<<t<<"\t"<<(double)L_B/L<<"\n";

        }

	delete []dependency_idx;
	delete []propensity;

}

void DEPmodel1d::MakeDependencyGraph(int * dependency_idx){


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
			dependency_idx[9*i + 0] = i;			  //diffusion_A
			dependency_idx[9*i + 1] = i-2;			  //infection
			dependency_idx[9*i + 2] = 4*nn[2*node] + 0;       //infection left
			dependency_idx[9*i + 3] = 4*nn[2*node] + channel; //diffusion_A, left
			dependency_idx[9*i + 4] = 4*nn[2*node+1] + 0;     //infection right
			dependency_idx[9*i + 5] = 4*nn[2*node+1]+channel; //diffusion_A, right
			dependency_idx[9*i + 6] = -1;  //nothing

		}else if(channel == 3){
			//DIFFUSION B
			dependency_idx[9*i + 0] = i;			  //diffusion_B
			dependency_idx[9*i + 1] = i-3;			  //infection
			dependency_idx[9*i + 2] = i-2;			  //decay
			dependency_idx[9*i + 3] = 4*nn[2*node] + 0;       //infection left
			dependency_idx[9*i + 4] = 4*nn[2*node] + 1;       //decay left
			dependency_idx[9*i + 5] = 4*nn[2*node]+channel;   //diffusion_B, left
			dependency_idx[9*i + 6] = 4*nn[2*node+1] + 0;     //infection right
			dependency_idx[9*i + 7] = 4*nn[2*node+1] + 1;     //decay right
			dependency_idx[9*i + 8] = 4*nn[2*node+1]+channel; //diffusion_B, right
		}


	}

}


void DEPmodel1d::UpdatePropensityTau(int* dependency_idx, double * propensity,IdxPriorityQueue& tau_q, int reaction_number , int diffusion_direction, double time, unsigned count){


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

		UpdateTau(dependency_idx,old_propensity, propensity, tau_q, reaction_number,diffusion_direction, time,count);
	}
	else if(channel == 1){
		for(int i=0;i<4;i++) old_propensity[i] = propensity[dependency_idx[9*reaction_number+i]];
		
		propensity[dependency_idx[9*reaction_number+0]] = particle_B[node]*decay_k2;
		propensity[dependency_idx[9*reaction_number+1]] = particle_A[node]*particle_B[node]*infection_k1;
		propensity[dependency_idx[9*reaction_number+2]] = particle_A[node]*diffusion_A;
		propensity[dependency_idx[9*reaction_number+3]] = particle_B[node]*diffusion_B;

		UpdateTau(dependency_idx,old_propensity, propensity, tau_q, reaction_number,diffusion_direction, time,count);
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
		UpdateTau(dependency_idx,old_propensity, propensity, tau_q, reaction_number,diffusion_direction, time,count);
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
		UpdateTau(dependency_idx,old_propensity, propensity, tau_q, reaction_number,diffusion_direction, time, count);
	}

	delete []old_propensity;

}
/*
void DEPmodel1d::UpdateTau(int* dependency_idx, double* old_propensity, double * propensity, IdxPriorityQueue& tau_q, int reaction_number,int diffusion_direction, double time){

	std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> rng(0.0,1.0);

//std::cout<<"\nUPDATE TAU======================================================\n";
//std::cout<<"old reaction number : "<<reaction_number<<", site : "<<reaction_number/4<<", channel : "<<reaction_number%4<<"\n";
	int channel = reaction_number%4;
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
//if(i==0) std::cout<<"old_tau[0] : "<<tau_q.heap_element[0].tau<<", ";
		tau_q.update(heap_idx, new_tau);
	}

	
//tau_q.print_tau();

//std::cout<<"\nnew_tau[0] : "<<tau_q.heap_element[0].tau<<"\n";
//std::cout<<"new reaction number : "<<tau_q.heap_element[0].reaction<<", site : "<<tau_q.heap_element[0].reaction/4<<", channel : "<<tau_q.heap_element[0].reaction%4<<"\n";
//std::cout<<"N_A : "<<L_A<<", N_B : "<<L_B<<"\n";
//std::cout<<"\nsite";
//for(int i =0; i< L ; i++)  std::cout<<i<<" ";
//std::cout<<"\nA : ";
//for(int i =0; i< L ; i++)  std::cout<<particle_A[i]<<" ";
//std::cout<<"\nB : ";
//for(int i =0; i< L ; i++)  std::cout<<particle_B[i]<<" ";
//std::cout<<"\n======================================================\n";

}
*/

/*
void DEPmodel1d::UpdateTau(int* dependency_idx, double* old_propensity, double * propensity, IdxPriorityQueue& tau_q, int reaction_number,int diffusion_direction, double time, unsigned count){

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
*/

void DEPmodel1d::UpdateTau(int* dependency_idx, double* old_propensity, double * propensity, IdxPriorityQueue& tau_q, int reaction_number,int diffusion_direction, double time, unsigned count){

	RNG rng;
        RNG::ctr_type c={{1,count,1,0}};
        RNG::ukey_type uk={{0,count}};
        uk[0] = seed; // some user_supplied_seed
        RNG::key_type k=uk;
	RNG::ctr_type r;

	int channel = reaction_number%4;

	if(channel < 2){
		for(int i = 0; i < 9; i++){

			if(i%4==0){
				++c[0]; 
				r = rng(c,k);
			}
			int new_reaction_number = dependency_idx[9*reaction_number+i];
			//since No Reaction is ... dependency_idx[] = -1;
			if(new_reaction_number == -1) break;

			int heap_idx = tau_q.heap_index[new_reaction_number];
			double new_tau = 2147483647;

			if(propensity[new_reaction_number]!=0){
				if(i==0)  new_tau = time - std::log(r123::u01fixedpt<double>(r.v[0]))/propensity[new_reaction_number];
				else{
					if(old_propensity[i]!=0) new_tau = time + old_propensity[i]*(tau_q.heap_element[heap_idx].tau - time)/propensity[new_reaction_number];
					else  new_tau = time - std::log(r123::u01fixedpt<double>(r.v[i%4]))/propensity[new_reaction_number];
				}
			}
			tau_q.update(heap_idx, new_tau);
		}	
	}
	else if(channel == 2){
		int * diff_index;
		if(diffusion_direction == 0) diff_index = new int[4]{0,1,2,3};
		else if(diffusion_direction == 1) diff_index = new int[4]{0,1,4,5};
				
		r = rng(c,k);

		for(int i = 0; i < 4; i++){
			int new_reaction_number = dependency_idx[9*reaction_number+diff_index[i]];
			//since No Reaction is ... dependency_idx[] = -1;

			int heap_idx = tau_q.heap_index[new_reaction_number];
			double new_tau = 2147483647;

			if(propensity[new_reaction_number]!=0){
				if(i==0)  new_tau = time - std::log(r123::u01fixedpt<double>(r.v[0]))/propensity[new_reaction_number];
				else{
					if(old_propensity[diff_index[i]]!=0) new_tau = time + old_propensity[diff_index[i]]*(tau_q.heap_element[heap_idx].tau - time)/propensity[new_reaction_number];
					else new_tau = time - std::log(r123::u01fixedpt<double>(r.v[i]))/propensity[new_reaction_number];

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
		
			if(i%4==0){
				++c[0]; 
				r = rng(c,k);
			}

			int new_reaction_number = dependency_idx[9*reaction_number+diff_index[i]];
			//since No Reaction is ... dependency_idx[] = -1;

			int heap_idx = tau_q.heap_index[new_reaction_number];
			double new_tau = 2147483647;

			if(propensity[new_reaction_number]!=0){
				if(i==0)  new_tau = time - std::log(r123::u01fixedpt<double>(r.v[0]))/propensity[new_reaction_number];
				else{
					if(old_propensity[diff_index[i]]!=0) new_tau = time + old_propensity[diff_index[i]]*(tau_q.heap_element[heap_idx].tau - time)/propensity[new_reaction_number];
					else new_tau = time - std::log(r123::u01fixedpt<double>(r.v[i%4]))/propensity[new_reaction_number];

				}
			}
			tau_q.update(heap_idx, new_tau);
		}
		delete []diff_index;

	}


}


//void DEPmodel1d::GillespieReaction(){
//	
//	std::random_device rd;
//        std::mt19937 mt(rd());
//        std::uniform_real_distribution<double> rng(0.0,1.0);
//
//
//	for(int i=0;i<L;i++){
//		if(particle_B[i] > 0){
//			double t_max = -std::log(1-rng(mt))/(particle_A[i]+particle_B[i]);
//			int tmp_B = 0, tmp_A = 0;
//
//			while(t_max>0){
//				double R_gill = infection_k1*particle_B[i]*particle_A[i] + decay_k2*particle_B[i];
//				double infection_p = infection_k1*particle_B[i]*particle_A[i];
//				infection_p /= R_gill;
//				double t_update = -std::log(1-rng(mt))/R_gill;
//
//				if(infection_p > rng(mt)){
//					tmp_B ++; tmp_A--;
//				}
//				else{
//					tmp_B--; tmp_A++;
//				}
//
//				t_max -= t_update;
//				if(particle_B[i]==-tmp_B) break;
//			}
//			particle_B[i] += tmp_B; particle_A[i] += tmp_A;
//			L_B += tmp_B; L_A += tmp_A;
//		}
//	}
//
//}

//void DEPmodel1d::SpatioTemporalPattern(int T_max,int *alltime_A, int * alltime_B){
//
//	for(int t=0;t<T_max;t++){
//		SimpleDiffusion(particle_A,diffusion_A);
//		SimpleDiffusion(particle_B,diffusion_B);
//		GillespieReaction();
//		for(int i =0; i<L;i++){
//		       	alltime_A[L*t + i] += particle_A[i];
//		       	alltime_B[L*t + i] += particle_B[i];
//		}
//	}
//}




void DEPmodel1d::DanielMaiaAlgorithm (int T_max){

	density_B = new double[T_max]{};
        survival = new int[T_max]{};
	density_Bsite = new double[T_max]{};
	radius_B = new double[T_max]{};


	double W_T = 0, W_A = diffusion_A*L_A, W_B = diffusion_B*L_B, W_R = decay_k2*L_B, W_I = 0;
	double current_time = 0;

	int a_max = particle_A[0],b_max = particle_B[0],ab_max = particle_A[0]*particle_B[0];
	for(int i=0; i<L;i++){

		if(a_max < particle_A[i])  a_max = particle_A[i];
		if(b_max < particle_B[i])  b_max = particle_B[i];
		if(ab_max < particle_A[i]*particle_B[i])  ab_max = particle_A[i]*particle_B[i];

		W_I += particle_AB[i]*infection_k1;
	}

	W_T = W_A + W_B + W_R + W_I;

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> rng;

	for(int t=0; t < T_max;){

		double rand_event = rng(mt);
		if(rand_event < W_A/W_T){
			//Diffusion of A
			while(true){
				int idx = (int)(rng(mt)*L);
				//double rand_accept = rng(mt);
				if((double)particle_A[idx]/a_max > rng(mt)){
					int tmp_A = 0;
					for(int j=0;j<particle_A[idx];j++){
						double rand_diff = rng(mt);
						if(diffusion_A/2 > rand_diff){
							tmp_A--; particle_A[nn[2*idx]]++;
						}
						else if(diffusion_A > rand_diff){
							tmp_A--; particle_A[nn[2*idx+1]]++;
						}
					}

					if(a_max == particle_A[idx]){
						int tmp_a_max = 0, tmp_ab_max = 0;
						for(int i=0;i<L;i++){
							if(i==idx) continue;
						       	if(tmp_a_max < particle_A[i]) tmp_a_max = particle_A[i];
						       	if(tmp_ab_max < particle_A[i]*particle_B[i]) tmp_ab_max = particle_B[i]*particle_A[i];
						}
						a_max = tmp_a_max; ab_max = tmp_ab_max;
					}

					particle_A[idx] += tmp_A;
					particle_AB[idx] = particle_A[idx]*particle_B[idx];
					particle_AB[nn[2*idx]] = particle_A[nn[2*idx]]*particle_B[nn[2*idx]];
					particle_AB[nn[2*idx+1]] = particle_A[nn[2*idx+1]]*particle_B[nn[2*idx+1]];
					L_A += tmp_A;

					break;
				}
			}

		}
		else if(rand_event < (W_A+W_B)/W_T){
			//Diffusion of B
			while(true){
				int idx = (int)(rng(mt)*L);
				//double rand_accept = rng(mt);
				if((double)particle_B[idx]/b_max > rng(mt)){
					int tmp_B = 0;
					for(int j=0;j<particle_B[idx];j++){
						double rand_diff = rng(mt);
						if(diffusion_B/2 > rand_diff){
							tmp_B--; particle_B[nn[2*idx]]++;
						}
						else if(diffusion_B > rand_diff){
							tmp_B--; particle_B[nn[2*idx+1]]++;
						}
					}
			
					if(b_max == particle_B[idx]){
						int tmp_b_max = 0, tmp_ab_max = 0;
						for(int i=0;i<L;i++){
							if(i==idx) continue;
						       	if(tmp_b_max < particle_B[i]) tmp_b_max = particle_B[i];
						       	if(tmp_ab_max < particle_A[i]*particle_B[i]) tmp_ab_max = particle_B[i]*particle_A[i];
						}
						b_max = tmp_b_max; ab_max = tmp_ab_max;
					}

					particle_B[idx] += tmp_B;
					particle_AB[idx] = particle_A[idx]*particle_B[idx];
					particle_AB[nn[2*idx]] = particle_A[nn[2*idx]]*particle_B[nn[2*idx]];
					particle_AB[nn[2*idx+1]] = particle_A[nn[2*idx+1]]*particle_B[nn[2*idx+1]];
					L_B += tmp_B;

					break;
				}
			}


		}
		else if(rand_event < (W_A+W_B+W_R)/W_T){
			//Decaying of B
			while(true){
				int idx = (int)(rng(mt)*L);

				if((double)particle_B[idx]/b_max > rng(mt)){
					
					int tmp_B = 0;
					for(int j=0;j<particle_B[idx];j++){
						double rand_decay = rng(mt);
						if( decay_k2 > rand_decay) tmp_B --;
					}
			
					if(b_max == particle_B[idx]){
						int tmp_b_max = 0, tmp_ab_max = 0;
						for(int i=0;i<L;i++){
							if(i==idx) continue;
						       	if(tmp_b_max < particle_B[i]) tmp_b_max = particle_B[i];
						       	if(tmp_ab_max < particle_A[i]*particle_B[i]) tmp_ab_max = particle_B[i]*particle_A[i];
						}
						b_max = tmp_b_max; ab_max = tmp_ab_max;
					}

					particle_B[idx] += tmp_B;
					particle_AB[idx] = particle_A[idx]*particle_B[idx];
					particle_AB[nn[2*idx]] = particle_A[nn[2*idx]]*particle_B[nn[2*idx]];
					particle_AB[nn[2*idx+1]] = particle_A[nn[2*idx+1]]*particle_B[nn[2*idx+1]];
				
					L_B += tmp_B;

					break;

				}

			}

		}
		else{
			//Infection
			while(true){
				int idx = (int)(rng(mt)*L);

				if( (double)particle_AB[idx]/ab_max > rng(mt)){

					int tmp_A = 0, tmp_B = 0;
					for(int j=0;j<particle_A[idx];j++){
						if(infection_k1 > rng(mt)){
						       	tmp_A --;
							tmp_B ++;
						}
					}
					
					if(a_max == particle_A[idx]){
						int tmp_a_max = 0, tmp_ab_max = 0;
						for(int i=0;i<L;i++){
							if(i==idx) continue;
						       	if(tmp_a_max < particle_A[i]) tmp_a_max = particle_A[i];
						       	if(tmp_ab_max < particle_A[i]*particle_B[i]) tmp_ab_max = particle_B[i]*particle_A[i];
						}
						a_max = tmp_a_max; ab_max = tmp_ab_max;
					}
				
					if(b_max == particle_B[idx]){
						int tmp_b_max = 0, tmp_ab_max = 0;
						for(int i=0;i<L;i++){
							if(i==idx) continue;
						       	if(tmp_b_max < particle_B[i]) tmp_b_max = particle_B[i];
						       	if(tmp_ab_max < particle_A[i]*particle_B[i]) tmp_ab_max = particle_B[i]*particle_A[i];
						}
						b_max = tmp_b_max; ab_max = tmp_ab_max;
					}

					particle_A[idx] += tmp_A; particle_B[idx] += tmp_B;
				
					particle_AB[idx] = particle_A[idx]*particle_B[idx];
					particle_AB[nn[2*idx]] = particle_A[nn[2*idx]]*particle_B[nn[2*idx]];
					particle_AB[nn[2*idx+1]] = particle_A[nn[2*idx+1]]*particle_B[nn[2*idx+1]];
					
					L_A += tmp_A; L_B += tmp_B;

					break;
				}
			}

		}

		current_time += 1/W_T;
		if(current_time >= 1){
		      	density_B[t] = (double)L_B/L;
			t++;
			current_time = 0;
		}

	}

}



template <typename T>
double avg_sample (T * sample, int T_max){

	double avg = 0;
	for (int i = 0; i < T_max ; i++){
		avg+=sample[i];
	}
	return avg/T_max;
}

template <typename T>
void binning (T * sample, int T_max){

	 int layer_max = std::log(T_max)/std::log(2) - 1;
	 if(layer_max <= 0) std::cout<<"ERROR, setting layer size"<<std::endl;

	 std::cout<<"layer_max = "<<layer_max<<std::endl;

	 double* err = new double[layer_max]{};

	 for(int layer = 0; layer < layer_max; layer++){

		 double avg = avg_sample(sample, T_max);
		 for (int i=0;i<T_max;i++){
			err[layer] += (sample[i] -avg)*(sample[i]-avg);
			if(i%2==1)  sample[i/2] = (sample[i] + sample[i-1])/2.0; 
		 }

		 err[layer] /= T_max*(T_max-1);
		 err[layer] = std::sqrt(err[layer]);
		 T_max /= 2;

		 std::cout << err[layer] << std::endl;
	 }

	 delete []err;

}


