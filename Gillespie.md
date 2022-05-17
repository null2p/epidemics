Gillespie Algorithm
===================
이 글은 Gillespie Algorithm의 Direct Method, Next Reaction Method를 설명하는 글입니다.
Michael A. Gibson, Jehosua Bruck의 논문을 기반으로 요약 작성되었습니다[^1]. 이해가 어렵거나 설명이 부족한 부분은 원문을 참고하시기 바랍니다.

Direct Method
-------------
한 시스템이 어떤 상태에 있을 때, 다음 스텝으로 진행하기 위해선 두가지 정보가 필요하다.
1. 어떤 Reaction이 일어날지
2. 언제 Reaction이 일어날지

다음 Reaction이 &mu;이고 &tau;시간에 일어날 확률 밀도 <img src="https://render.githubusercontent.com/render/math?math={\large{P(\mu,\tau)}}##gh-light-mode-only"><img src="https://render.githubusercontent.com/render/math?math={\color{white}\large{P(\mu,\tau)}}#gh-dark-mode-only">가 확률적으로 정해진다.
Gillespie는 확률밀도가 다음과 같다고 밝혔다.[^2]
<img src="https://render.githubusercontent.com/render/math?math={\Large{P(\mu,\tau)d\tau=a_{\mu}\exp(-\tau \sum_j a_j)d\tau}}##gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\Large{P(\mu,\tau)d\tau=a_{\mu}\exp(-\tau \sum_j a_j)d\tau}}#gh-dark-mode-only">

위 식을 &tau; 에 대해 0부터 무한대까지 적분해 어떤 Reaction이 일어날지 확률적으로 알 수 있다:

<img src="https://render.githubusercontent.com/render/math?math={\Large{P(\mu)=a_{\mu}/\sum_j a_j}}##gh-light-mode-only"><img src="https://render.githubusercontent.com/render/math?math={\color{white}\Large{P(\mu)=a_{\mu}/\sum_j a_j}}#gh-dark-mode-only">

그리고 &mu;에 대해 적분해 시간에 대한 확률분포도 알 수 있다:

<img src="https://render.githubusercontent.com/render/math?math={\Large{P(\tau)d\tau = \sum_j a_j\exp(-\tau \sum_j a_j)d\tau}}##gh-light-mode-only"><img src="https://render.githubusercontent.com/render/math?math={\color{white}\Large{P(\tau)d\tau = \sum_j a_j\exp(-\tau \sum_j a_j)d\tau}}#gh-dark-mode-only">


위 두 확률분포를 이용한 것이 Direct Method 이다.

**Algorithm**
1. 초기 설정을 한다. t=0, 입자 수 설정
2. 모든 i에 대해 propensity function a<sub>i</sub> 을 계산한다.
3. <img src="https://render.githubusercontent.com/render/math?math={{P(\mu)=a_{\mu}/\sum_j a_j}}##gh-light-mode-only"><img src="https://render.githubusercontent.com/render/math?math={\color{white}{P(\mu)=a_{\mu}/\sum_j a_j}}#gh-dark-mode-only">를 이용해 어떤 Reaction &mu;가 일어날지 정한다
4. <img src="https://render.githubusercontent.com/render/math?math={{P(\tau)d\tau = \sum_j a_j\exp(-\tau \sum_j a_j)d\tau}}##gh-light-mode-only"><img src="https://render.githubusercontent.com/render/math?math={\color{white}{P(\tau)d\tau = \sum_j a_j\exp(-\tau \sum_j a_j)d\tau}}#gh-dark-mode-only">를 이용해 &tau;를 정한다. (propensity function a<sub>i</sub>들의 합을 계수로 하는 지수분포 이용한다)
5. Reaction &mu;를 반영해 반응한 입자수들을 업데이트하고 t는 t+&tau;로 바꿔준다.
6. 2번으로 돌아가서 알고리듬을 반복한다.

First Reaction Method
----
Direct Method는 &mu;와 &tau;를 직접 구한다는 점에서 'Direct'라는 이름과 통하는 점이 있다. 이번의 방법론은 각각의 Reaction에 대해 예상되는 시간 $tau;<sub>i</sub>를 정한다. $tau;<sub>i</sub>는 어떤 Reaction이 일어나고 난 후에 다음 첫번째 Reaction이 일어날 때까지 걸리는 시간이다. 

**Algorithm**
1. 초기 설정을 한다. t=0, 입자 수 설정
2. 모든 i에 대해 propensity function a<sub>i</sub>을 계산한다.
3. 모든 i에 대해 예상되는 시간(putative time) &tau;<sub>i</sub>를 계수가 a<sub>i</sub>인 지수함수 분포를 이용해 만든다.
4. &tau;<sub>&mu;</sub>가 최소인 Reaction &mu;를 찾는다.
5. Reaction &mu;를 시행해 입자 수들을 업데이트하고 t=t+&tau;<sub>&mu;</sub>로 시간을 바꾼다.
6. 2번으로 돌아가 알고리듬을 반복한다.

Direct Method와 First Reaction Method는 확률적으로 동등하다. &mu;와 &tau;를 구할 때의 확률분포는 같기 때문이다. 여기서 랜덤넘버를 매번 Reaction 수 만큼 불러들이고 최소의 &tau;<sub>&mu;</sub>를 찾는데 많은 시간을 쓰게 된다. 이를 개선한 방법론이 Next Reaction Method이다.

Next Reaction Method
----
First Reaction Method는 한 Reaction이 일어나면 3가지의 반복작업이 항상 필요하다. (1) 모든 Reaction에 대한 Propensity function, a<sub>i</sub>; (2) 모든 Reaction에 대한 다음 예상 시간 &tau;<sub>i</sub>; (3) 가장 작은 예상시간 &tau;<sub>&mu;</sub> 찾기.
Next Reaction Method는 몇몇 반복작업을 없애서 계산 시간을 개선한 버전이다. 주요한 변화는 다음과 같다.

* a<sub>i</sub>뿐만 아니라 &tau;<sub>i</sub>를 저장한다
* a<sub>i</sub>와 &tau;<sub>i</sub>를 계산할 때 Reaction에 영향을 받는 것들만 업데이트 한다. 한 Reaction이 어떤 입자들에 영향을 미치는지 미리 Dependency Graph를 그려서 파악해둬야 한다.

위 과정에서 한 Reaction에 영향을 받는 입자들의 a<sub>i</sub>와 &tau;<sub>i</sub>들만 업데이트 하기 때문에 재사용 되는 &tau;가 존재하게 된다. 이는 랜덤넘버를 재사용 하는 것과 마찬가지이기 때문에 부적절한 몬테카를로 시뮬레이션이 아닌가 싶지만, 두번의 간단한 변환을 이용하면 수학적으로 적절하다는 것을 밝힐 수 있다.

* 한 Reaction이 일어나고 난 다음에 상대적인 시간을 이용하는 것이 아니라 절대적인 시간을 이용한다.
* Indexed Priroity Queue를 이용해 a<sub>i</sub>와 &tau;<sub>i</sub>를 저장한다.

**Algorithm**
1. 초기 설정을 한다. 
    - t=0, 입자 수 설정, Dependency Graph를 만든다
    - 모든 i에 대해 propensity function a<sub>i</sub>을 계산한다.
    - 모든 i에 대해 예상되는 시간(putative time) &tau;<sub>i</sub>를 계수가 a<sub>i</sub>인 지수함수 분포를 이용해 만든다.
    - 모든 &tau;<sub>i</sub>를 Indexed Priority Queue에 저장한다.
2. Indexed Priority Queue에 저장된 &tau;중 가장 작은 &tau;를 고른다. 이것이 다음 Reaction &mu;가 된다.
3. Reaction &mu;에 영향을 받는 입자들의 수를 바꿔주고 시간은 t = &tau;로 설정한다.
4. Dependency Graph 내에서 영향을 받는 입자들을 고려해,
   - a<sub>i</sub>들을 업데이트 한다.
   - i==&mu;라면, 새로운 랜덤넘버가 필요하다. a<sub>&mu;</sub>를 계수로 하는 지수함수분포를 이용해 &rho;를 정하고 &tau;<sub>&mu;</sub> = &rho; + t;로 업데이트 한다.
   - i ≠ &mu;라면, &tau;<sub>i</sub> = (a<sub>i,old</sub>/a<sub>i,new</sub>)(&tau;<sub>i</sub>-t) + t
   - 위에서 새로 정한 &tau;<sub>i</sub>들을 업데이트 한다.
5. 2번 스텝으로 돌아간다.

[^1]:Gibson, M. A., & Bruck, J. (2000). Efficient exact stochastic simulation of chemical systems with many species and many channels. [The journal of physical chemistry A, 104(9), 1876-1889.](https://pubs.acs.org/doi/pdf/10.1021/jp993732q)
[^2]:Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. [The journal of physical chemistry, 81(25), 2340-2361.](https://pubs.acs.org/doi/pdf/10.1021/j100540a008)
