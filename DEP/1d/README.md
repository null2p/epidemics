# Diffusive Epidemic Process (DEP)

DEP는 DEP_algorithm.cpp내에 구현되어 있습니다. 여기서는 1d DEP model을 이용하여 2022년에 출판된 논문[^1]의 그림들을 구현하도록 하겠습니다. Random Number Generator로는 D. E. Shaw Research의 Random123 라이브러리 내의 Philox를 사용합니다. Random123의 RNG들이 어떻게 구현되는지는 다음 논문에 자세히 나와 있습니다 [^2]. 또한 Github 링크도 존재하여 간단한 rng 원리와 설치법도 설명되어 있습니다 [^3].


DEP_algorithm.h
-----

```c++
class DEPmodel
{
    public:
        const double discrete_time = 1.0;
        const double diffusion_A = 2, diffusion_B = 1, infection_k1 = 0.2, decay_k2 = 1;
        //const double diffusion_A = 1*discrete_time, diffusion_B = 0.5*discrete_time, infection_k1 = 0.2*discrete_time, decay_k2 = 1*discrete_time;
        DEPmodel(int L, int total_number, int z, std::string initial_dist, int seed);
};
```
diffusion_A, diffusion_B, infection_k1, decay_k2는 각각 건강한 입자(A)와 감염된 입자(B)의 확산계수, 감염률, 치료율을 나타냅니다. 논문[^1]에서는 diffusion_A = 1, diffusion_B = 0.5, infection_k1 = 0.2, decay_k2 = 1 입니다. 하지만 Gillespie Algorithm내에서 Propensity Function으로 사용될때 diffusion은 방향까지 고려해줘야하므로 1d 상에서는 2를 곱해줘야합니다. 따라서 위 코드 상에서 const double diffusion_A = 2, diffusion_B = 1, infection_k1 = 0.2, decay_k2 = 1;로 설정되었습니다. 

discrete_time은 Gillespie Algorithm이 아닌 Discrete time method를 사용하는 경우 필요합니다. dicrete time method 에서는 propensity function을 이용하지 않으므로 기존에 설정된 값들에 discrete time을 곱해주는 꼴로 설정됩니다.

DEPmodel 클래스 내에서는 Initialization에 집중하고 실제 함수들은 DEPmodel1d에 구현되어 있습니다. DEPmodel1d라는 이름에서 볼 수 있듯이 2d, 3d, network...의 확장을 고려하여 만들었습니다. 현재는 1d만 존재합니다. DEPmodel1d는 차원에 관계없는 DEPmodel을 상속받습니다.

DEPmodel1d -> L, total_number, z, T_max, init_dist
NextReactionMethod -> T_max

[^1]: Polovnikov, B., Wilke, P., & Frey, E. (2022). [Subdiffusive Activity Spreading in the Diffusive Epidemic Process. Physical Review Letters, 128(7), 078302.](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.128.078302)

[^2]: Salmon, J. K., Moraes, M. A., Dror, R. O., & Shaw, D. E. (2011, November). [Parallel random numbers: as easy as 1, 2, 3. In Proceedings of 2011 international conference for high performance computing, networking, storage and analysis (pp. 1-12).](https://dl.acm.org/doi/pdf/10.1145/2063384.2063405)

[^3]: https://github.com/DEShawResearch/random123
