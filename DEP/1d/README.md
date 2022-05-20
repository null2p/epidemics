# Diffusive Epidemic Process (DEP)

DEP는 DEP_algorithm.cpp내에 구현되어 있습니다. 여기서는 1d DEP model을 이용하여 2022년에 출판된 논문[^1]의 그림들을 구현하도록 하겠습니다. Random Number Generator로는 D. E. Shaw Research의 Random123 라이브러리 내의 Philox를 사용합니다. Random123의 RNG들이 어떻게 구현되는지는 다음 논문에 자세히 나와 있습니다 [^2]. 또한 Github 링크도 존재하여 간단한 rng 원리와 설치법도 설명되어 있습니다 [^3].


DEP_algorithm.h
-----
DEP_algorithm.h 파일부터 보겠습니다. 

'''
                const double discrete_time = 1.0;
                //const double diffusion_A = 2/3.0, diffusion_B = 2/3.0, infection_k1 = 1, decay_k2 = 0.03278;
                //const double diffusion_A = 0.5, diffusion_B = 0.25, infection_k1 = 0.5, decay_k2 = 0.5;
                //const double diffusion_A = 1, diffusion_B = 0.5, infection_k1 = 1.0/6, decay_k2 = 5.0/6;
                //const double diffusion_A = 1, diffusion_B = 0.5, infection_k1 = 0.2, decay_k2 = 1-std::exp(-1);
                //const double diffusion_A = 2, diffusion_B = 1, infection_k1 = 0.2, decay_k2 = 1;
                const double diffusion_A = 0.5*discrete_time, diffusion_B = 0.25*discrete_time, infection_k1 = 0.5*discrete_time, decay_k2 = 0.5*discrete_time;

'''


[^1]: Polovnikov, B., Wilke, P., & Frey, E. (2022). [Subdiffusive Activity Spreading in the Diffusive Epidemic Process. Physical Review Letters, 128(7), 078302.](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.128.078302)

[^2]: Salmon, J. K., Moraes, M. A., Dror, R. O., & Shaw, D. E. (2011, November). [Parallel random numbers: as easy as 1, 2, 3. In Proceedings of 2011 international conference for high performance computing, networking, storage and analysis (pp. 1-12).](https://dl.acm.org/doi/pdf/10.1145/2063384.2063405)

[^3]: https://github.com/DEShawResearch/random123
