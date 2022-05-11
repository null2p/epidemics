Gillespie Algorithm
===================
이 글은 Gillespie Algorithm의 Direct Method, Next Reaction Method를 설명하는 글입니다.
Michael A. Gibson, Jehosua Bruck의 논문을 기반으로 작성되었습니다[^1].

Direct Method
-------------
한 시스템이 어떤 상태에 있을 때, 다음 스텝으로 진행하기 위해선 두가지 정보가 필요하다.
1. 어떤 Reaction이 일어날지
2. 언제 Reaction이 일어날지

<img src="https://render.githubusercontent.com/render/math?math={\Large{P(\mu,\tau)d\tau=a_{\mu}\exp(-\tau \sum_j a_j)d\tau}}##gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\Large{P(\mu,\tau)d\tau=a_{\mu}\exp(-\tau \sum_j a_j)d\tau}}#gh-dark-mode-only">



[^1]:Gibson, M. A., & Bruck, J. (2000). Efficient exact stochastic simulation of chemical systems with many species and many channels. [The journal of physical chemistry A, 104(9), 1876-1889.](https://pubs.acs.org/doi/pdf/10.1021/jp993732q)
