# Contact Process with Temporal Feedback

How to Compile
-----

```sh
make density
```

How to run
------
Modify Bash file (CPTF_density.bash)

```sh
#!/bin/bash

n="30"                    #Number of CPU to run
L="1000"                  #size of lattice
p="0.46436"               #1-p is infection rate
theta="1.5"               #theta is related with the life time of B particles
R_btoa="0.5"              #the probability of B -> A
maxT="100000"             #maximal time for this dynamics
Db="0.5"                  #Probability of direction bias
seed="single"             #seed type : 'single' or 'half'
parallelization="30000"   #Number of Samples(Ensembles)
```

Run bash file
```
bash CPTF_density.bash
```

Result Files path
```
./rho_scaling/seed_*/theta*/L*T*ens*p*theta*R*Db*.txt
```

How to Draw figures
----

modify a file ./rho_scaling/rho*.py

You need to change the name of text file to read.

```python3
np.loadtxt("./seed_single/theta1/L1000T100000ens3000p0.49598theta1R0.5Db0.5.txt", usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)
```

![rho_t_example](https://user-images.githubusercontent.com/68416208/170169590-fd893f6b-5f7e-410d-9f23-ccd6fd02b998.png)


