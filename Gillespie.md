Gillespie Algorithm
===================
이 글은 Gillespie Algorithm의 Direct Method, Next Reaction Method를 설명하는 글입니다.
Michael A. Gibson, Jehosua Bruck의 논문을 기반으로 요약 작성되었습니다[^1]. 이해가 어렵거나 설명이 부족한 부분은 원문을 참고하시기 바랍니다.

Direct Method
-------------
한 시스템이 어떤 상태에 있을 때, 다음 스텝으로 진행하기 위해선 두가지 정보가 필요하다.
1. 어떤 Reaction이 일어날지
2. 언제 Reaction이 일어날지

다음 Reaction이 &mu;이고 &tau;시간에 일어날 확률 밀도 ${\large{P(\mu,\tau)}}$가 확률적으로 정해진다.
Gillespie는 확률밀도가 다음과 같다고 밝혔다.[^2]

$${\Large{P(\mu,\tau)d\tau=a_{\mu}\exp(-\tau \sum_j a_j)d\tau}}$$

위 식을 &tau; 에 대해 0부터 무한대까지 적분해 어떤 Reaction이 일어날지 확률적으로 알 수 있다:

$${\Large{P(\mu)=a_{\mu}/\sum_j a_j}}$$

그리고 &mu;에 대해 적분해 시간에 대한 확률분포도 알 수 있다:

$${\Large{P(\tau)d\tau = \sum_j a_j\exp(-\tau \sum_j a_j)d\tau}}$$


위 두 확률분포를 이용한 것이 Direct Method 이다.

**Algorithm**
1. 초기 설정을 한다. t=0, 입자 수 설정
2. 모든 i에 대해 propensity function a<sub>i</sub> 을 계산한다.
3. ${{P(\mu)=a_{\mu}/\sum_j a_j}}$를 이용해 어떤 Reaction &mu;가 일어날지 정한다
4. ${{P(\tau)d\tau = \sum_j a_j\exp(-\tau \sum_j a_j)d\tau}}$를 이용해 &tau;를 정한다. (propensity function a<sub>i</sub>들의 합을 계수로 하는 지수분포 이용한다)
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

Re-using &tau;<sub>i</sub>
-----
 Next Reaction Method에서는 한 Reaction이 일어나고 나면 영향을 받는 입자들의 a<sub>i</sub>, &tau;<sub>i</sub>들만 업데이트 하게 된다. 그 외의 값들은 모두 재사용하게 되는데 이게 수학적으로 타당한 것인지 의문이 들것이다. 이번 챕터에서는 &tau;<sub>i</sub>를 재사용 하는 것이 왜 타당한지 수학적으로 밝힌다.
 
 첫번째로, First Reaction Method와 Next Reaction Method는 time이 상대적(First)이냐, 절대적(Next)이냐의 차이가 존재한다. random variable을 R<sub>i</sub>라 하자. R<sub>i</sub>=exp(a<sub>i</sub>) 이며 R<sub>i</sub>의 density는 P<sub>R<sub>i</sub></sub>(&tau;) = &theta;(&tau;)a<sub>i</sub>exp(-a<sub>i</sub>&tau;)이다. 절대적 시간은 n번째 iteration에서의 시간에 상대적 시간을 더한 꼴이다. T<sub>i</sub> = R<sub>i</sub> + t<sub>n</sub>. 따라서 T<sub>i</sub>의 density는 RVT(random variable transformation)을 이용해 구할 수 있다. 
 
${\Large{P_{T_i}(\tau)=\int_{-\infty}^{\infty}P_{R_i}(\tau')\delta(\tau - [ \tau' + t_n ])d\tau' = P_{R_i}(\tau-t_n) = \theta(\tau - t_n)a_i \exp(-a_i (\tau-t_n)) }}$...(7)

따라서, First Reaction Method의 절대적인 시간은 Direct method의 상대적인 시간과 동등하다. 이제 Next Reaction Method에서도 위 7번 식을 만족하는지 보일 것이다. Next Reaction Method의 알고리듬 1번은 First Reaction Method의 알고리듬 1~3번까지와 같기 때문에 위 7번 식을 만족할 것이다. 우리는 그 이후에도 이를 만족한다고 가정하고 이를 증명한다.

**Theorem 1.** Next Reaction Method의 알고리듬 Step 2가 시작할 때 7번식을 만족한다고 가정하자, 모든 i ≠ &mu;에 대해, &tau;<sub>i</sub>는 다음과 같이 분포되어 있다.

${\Large{P_{T_i}(\tau) =  \theta(\tau - t_{n+1})a_i \exp(-a_i (\tau-t_{n+1})) }}$
... (8)

**Proof.** Next Reaction Method 알고리듬 내에서 &tau;<sub>i</sub>는 식 7번에 따라 분포되어있고, 그 중 최소의 &tau;인 &tau;<sub>&mu;</sub>를 찾는다. T<sub>&mu;</sub>는 &tau;<sub>&mu;</sub>가 되고 다른 &tau;<sub>i</sub>들은 모두 &tau;<sub>&mu;</sub>보다 클 것이다. 그러므로 다른 T<sub>i</sub>들은 Pr(T<sub>i</sub> > u |T<sub>i</sub> >  &tau;<sub>&mu;</sub>)에 따라 분포 될 것이다. 여기서 u가 &tau;<sub>&mu;</sub>보다 큰지 작은지에 따라 분류할 수 있다. 따라서 첫번째로 u가 &tau;<sub>&mu;</sub>보다 큰 경우를 보면, 분모는 Pr$(T_i > u)$가 되고 $\exp(-a_{i,n}(u-t_n))/\exp(-a_{i,n}(\tau_{\mu}-t_n)) = \exp(-a_{i,n}(i-\tau_{\mu}))$를 얻게 된다. 다음으로 $u$가 $\tau_{\mu}$보다 작거나 같은 경우이다. 분모는 Pr$(T_i>\tau_{\mu})$가 되어서 1이 된다. 따라서 위 Theroem을 만족하게 된다.


[^1]:Gibson, M. A., & Bruck, J. (2000). Efficient exact stochastic simulation of chemical systems with many species and many channels. [The journal of physical chemistry A, 104(9), 1876-1889.](https://pubs.acs.org/doi/pdf/10.1021/jp993732q)
[^2]:Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. [The journal of physical chemistry, 81(25), 2340-2361.](https://pubs.acs.org/doi/pdf/10.1021/j100540a008)
