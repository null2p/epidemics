# Next Reaction method on 1d DEP

| Reaction | a<sub>&mu;</sub> | DepOn(a<sub>&mu;</sub>)| Affects(&mu;)|
| ------------ | ------------- | -----------------|----------|
| A+B :arrow_right: B+B |  k<sub>1</sub> #A #B| A, B|A,B|
| B :arrow_right: A |  k<sub>2</sub> #B  |B|A,B |
| A :arrow_right: 0 |  D<sub>A</sub> #A  |A<sub>i</sub>| A<sub>i-1</sub>,A<sub>i</sub>,A<sub>i+1</sub> |
| B :arrow_right: 0 |  D<sub>B</sub> #B  |B<sub>i</sub>| B<sub>i-1</sub>,B<sub>i</sub>,B<sub>i+1</sub> |

위 표는 Dependency Graph를 그리기 전, 감염, 치료, A의 확산, B의 확산의 Propensity function은 무엇인지(a<sub>&mu;</sub>), 
