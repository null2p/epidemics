# Spatio-Temporal Pattern for Contact Process with Temporal Feedback

How to Compile
-----

```sh
make density
```

How to run
------
Modify Bash file (CPTF_pattern.bash)

```sh
#!/bin/bash

n="4"                       #mpi size
L="128"                     #size of lattice
p="0.49598"                 #1-p is infection rate
theta="1"                   #theta is related with the life of B particles
R_btoa="0.5"                #the probability of B -> A at the first time
newR="0.5"                  #the probability of B -> A at the middle of time
initT="0"                   #initial time for saving state
midT="0"                    #middle of time : you can change R_btoa at the middle of time
finalT="256"                #maximal time for this dynamics
D_b='0.5'                   #direction bias -> branching bias
seed_type='single'          #seed type = 'single' or 'half'
```

Run bash file
```
bash CPTF_pattern.bash
```

Result Files path
```
./seed_*/theta*/Apartilce_L*T*ens*p*theta*R*Db*.txt
./seed_*/theta*/Bpartilce_L*T*ens*p*theta*R*Db*.txt
./seed_*/theta*/Scaling_L*T*ens*p*theta*R*Db*.txt
```


How to Draw figures
----

modify a file ./rho_scaling/rho*.py

You need to change the name of text file to read.

```python3
np.loadtxt("./seed_single/theta1/L1000T100000ens3000p0.49598theta1R0.5Db0.5.txt", usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)
```

![pattern_example](https://user-images.githubusercontent.com/68416208/170174754-78d2c372-6c30-4e0e-9e2d-6a6a6ed4ae66.png)

