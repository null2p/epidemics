# Next Reaction method on 1d DEP

| Reaction | a<sub>&mu;</sub> | DepOn(a<sub>&mu;</sub>)| Affects(&mu;)|
| ------------ | ------------- | -----------------|----------|
| A+B :arrow_right: B+B |  k<sub>1</sub> #A<sub>i</sub> #B<sub>i</sub>| A<sub>i</sub>, B<sub>i</sub>|A<sub>i</sub>,B<sub>i</sub>|
| B :arrow_right: A |  k<sub>2</sub> #B<sub>i</sub>  |B<sub>i</sub>|A<sub>i</sub>,B<sub>i</sub> |
| A :arrow_right: 0 |  D<sub>A</sub> #A<sub>i</sub>  |A<sub>i</sub>| A<sub>i-1</sub>,A<sub>i</sub>,A<sub>i+1</sub> |
| B :arrow_right: 0 |  D<sub>B</sub> #B<sub>i</sub>  |B<sub>i</sub>| B<sub>i-1</sub>,B<sub>i</sub>,B<sub>i+1</sub> |

위 표는 Dependency Graph를 그리기 위해 감염, 치료, A의 확산, B의 확산의 
- Propensity function은 무엇인지 : a<sub>&mu;</sub>
- Reaction을 계산할때 필요한 것들은 무엇인지 : DepOn(a<sub>&mu;</sub>)
- Reaction이 일어난 후 영향을 받은 것들은 무엇인지 : Affects(&mu;)

정리한 표이다.


Dependency Graph
-----

![Dark image](https://user-images.githubusercontent.com/68416208/169630782-7e81ea8f-7ffb-4343-8c54-166a28aa7426.png#gh-dark-mode-only)
![Light image](https://user-images.githubusercontent.com/68416208/169630794-747f81f2-dbbc-4a90-bfed-97456ba82621.png#gh-light-mode-only)

위는 i번째 위치 와 (i+1)번째 위치의 Dependency Graph이다. A와 B의 확산이 i에서 i+1로 일어났다 가정하고 Dependency Graph를 그렸다. a<sub>i,j</sub>의 i는 노드 위치이며 j는 Reaction 타입을 의미한다. j=0 : 감염, j=1 : 치료, j=2 : A의 확산, j=3 : B의 확산이다.


Dependency Graph in code
-----

실제 1d DEP에서 코드가 어떻게 작성되었는지 보일 것이다. 모든 노드(L)마다 4개의 Reaction이 존재하기 때문에 L * 4 개의 설정이 필요하다. 또한 1차원 어레이를 이용해 순차적으로 Reaction에 영향을 받는 propensity function들을 가리키게 하였다.

```c++
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
```
