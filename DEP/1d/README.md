# Diffusive Epidemic Process (DEP)

DEP는 DEP_algorithm.cpp내에 구현되어 있습니다. 여기서는 1d DEP model을 이용하여 2022년에 출판된 논문[^1]의 그림들을 구현하도록 하겠습니다. Random Number Generator로는 D. E. Shaw Research의 Random123 라이브러리 내의 Philox를 사용합니다. Random123의 RNG들이 어떻게 구현되는지는 다음 논문에 자세히 나와 있습니다 [^2]. 또한 Github 링크도 존재하여 간단한 rng 원리와 설치법도 설명되어 있습니다 [^3].

컴파일 방법
-----

현재 디렉토리가 erwin 폴더일 때의 컴파일 방법입니다

```
mpic++ fig4.cpp ../DEP_algorithm.cpp
```

실행 방법
-----
```
mpirun -np $cpu스레드개수 a.out $L
```

결과물 예시 (L = 128)
-----
```
#t      density_B       P_survival
10      1.26278 1
12.5893 1.17769 1
15.8489 1.08269 1
19.9526 0.993012        1
25.1189 0.902051        1
31.6228 0.834759        1
39.8107 0.769541        0.9999
50.1187 0.7067  0.9998
63.0957 0.655439        0.998503
79.4328 0.60032 0.992315
100     0.54985 0.97994
125.893 0.504933        0.96008
158.489 0.463058        0.926447
199.526 0.419279        0.881737
251.189 0.383869        0.823353
316.228 0.353765        0.756088
398.107 0.318706        0.682735
501.187 0.28477 0.603892
630.957 0.257258        0.52515
794.328 0.228167        0.450798
1000    0.19727 0.381637
1258.93 0.171153        0.317864
1584.89 0.137741        0.25509
1995.26 0.106121        0.194711
2511.89 0.0747107       0.137026
3162.28 0.0513442       0.0928144
3981.07 0.0303565       0.0563872
5011.87 0.0159626       0.0289421
6309.57 0.00701722      0.0129741
7943.28 0.00263146      0.00469062
10000   0.000585548     0.0010978
12589.3 0.000137226     0.000199601
15848.9 0       0
```
![image](https://user-images.githubusercontent.com/68416208/169491209-671137bd-fa2b-48a8-b07a-6c478c8a95c5.png)

위 그림에서 diamond 들은 제 코드로 얻은 결과이며 회색의 원들은 Erwin 논문[^1]의 fig 4 입니다. 

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

DEPmodel은 L, total_number, z, initial_dist, seed가 필요하며 각각 의미하는 바는 다음과 같습니다.
  - L : Lattice size = 노드 수
  - total_number : 총 입자 개수 ( A + B )
  - z : 각 노드마다 연결된 선 (Edge) 수
  - initial_dist : 초기 입자 분포, single 또는 homogeneous
  - seed : RNG에 들어갈 시드

이제 erwin 폴대 내의 fig4.cpp 를 이용해 1d DEP를 실행해보겠습니다.

erwin/fig4.cpp
----

DEPmodel 클래스에서 필요한 값들과 NextReactionMethod 함수에서 필요한 T_max를 정해주면 됩니다. 논문 [^1]에 따라...
  - 10000개의 샘플 (parallelization = 10000)
  - 일차원 (z = 2)
  - &rho; = 6.765 (total_number = 6.765*L)
  - Homogeneous 시작 조건 (init = "homogeneous")
  - 최대 시간 (T_max = 100000)

```c++
    int parallelization = 10020, run_by_onecore = parallelization/NUM_WORKERS;
    int L = 128, z = 2,T_max = 100000;
    int total_number= 6.765*L;
    std::string init = "homogeneous";
```

[^1]: Polovnikov, B., Wilke, P., & Frey, E. (2022). [Subdiffusive Activity Spreading in the Diffusive Epidemic Process. Physical Review Letters, 128(7), 078302.](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.128.078302)

[^2]: Salmon, J. K., Moraes, M. A., Dror, R. O., & Shaw, D. E. (2011, November). [Parallel random numbers: as easy as 1, 2, 3. In Proceedings of 2011 international conference for high performance computing, networking, storage and analysis (pp. 1-12).](https://dl.acm.org/doi/pdf/10.1145/2063384.2063405)

[^3]: https://github.com/DEShawResearch/random123
