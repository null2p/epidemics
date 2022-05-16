Gillespie Algorithm
===================
이 글은 Gillespie Algorithm의 Direct Method, Next Reaction Method를 설명하는 글입니다.
Michael A. Gibson, Jehosua Bruck의 논문을 기반으로 작성되었습니다[^1].

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
<img src="https://render.githubusercontent.com/render/math?math={\Large{P(\mu)=a_{\mu}/\sum_j a_j}}##gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\Large{P(\mu)=a_{\mu}/\sum_j a_j}}#gh-dark-mode-only">

그리고 &mu;에 대해 적분해 시간에 대한 확률분포도 알 수 있다:
<img src="https://render.githubusercontent.com/render/math?math={\Large{P(\tau)d\tau = \sum_j a_j\exp(-\tau \sum_j a_j)d\tau}}##gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\Large{P(\tau)d\tau = \sum_j a_j\exp(-\tau \sum_j a_j)d\tau}}#gh-dark-mode-only">


위 두 확률분포를 이용한 것이 Direct Method 이다.

**Algorithm**
1. 초기 설정을 한다. t=0, 입자 수 설정
2. 모든 i에 대해 propensity function a<sub>i</sub> 을 계산한다.
3. <img src="https://render.githubusercontent.com/render/math?math={{P(\mu)=a_{\mu}/\sum_j a_j}}##gh-light-mode-only"><img src="https://render.githubusercontent.com/render/math?math={\color{white}{P(\mu)=a_{\mu}/\sum_j a_j}}#gh-dark-mode-only">를 이용해 어떤 Reaction &mu;가 일어날지 정한다
4. <img src="https://render.githubusercontent.com/render/math?math={\Large{P(\tau)d\tau = \sum_j a_j\exp(-\tau \sum_j a_j)d\tau}}##gh-light-mode-only"><img src="https://render.githubusercontent.com/render/math?math={\color{white}\Large{P(\tau)d\tau = \sum_j a_j\exp(-\tau \sum_j a_j)d\tau}}#gh-dark-mode-only">를 이용해 &tau;를 정한다.
5. Reaction &mu;를 반영해 반응한 입자수들을 업데이트하고 t는 t+&tau;로 바꿔준다.
6. 2번으로 돌아가서 알고리듬을 반복한다.

Next Reaction Method
----




[^1]:Gibson, M. A., & Bruck, J. (2000). Efficient exact stochastic simulation of chemical systems with many species and many channels. [The journal of physical chemistry A, 104(9), 1876-1889.](https://pubs.acs.org/doi/pdf/10.1021/jp993732q)
[^2]:Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. [The journal of physical chemistry, 81(25), 2340-2361.](https://pubs.acs.org/doi/pdf/10.1021/j100540a008)
