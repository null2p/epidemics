#!/bin/bash

n="30"                  #mpi size
L="1000"                  #size of lattice
#p="0.4495"              #1-p is infection rate
p="0.46436"              #1-p is infection rate
#p="0.461"              #1-p is infection rate
theta="1.5"                #theta is related with the life time of B particles
R_btoa="0.5"             #the probability of B -> A
maxT="100000"            #maximal time for this dynamics
Db="0.5"		 #Probability of direction bias
seed="single"	 	 #seed type : 'single' or 'half'
parallelization="3000000"


#Db_list=(0.75 1 0.5 0.5 1)
#R_list=(0.5 0.75 0.75 0.25 0.25)

Db_list=(0.5 1)
R_list=(0.25 0.25)

for idx in {0..1}
do
        echo R:${R_list[idx]} Db:${Db_list[idx]}
        time mpirun -np $n --oversubscribe ./a.out $L 0.51065 1 ${R_list[idx]}  $maxT ${Db_list[idx]} $seed $parallelization &
        wait
        time mpirun -np $n --oversubscribe ./a.out $L 0.473 1.5 ${R_list[idx]}  $maxT ${Db_list[idx]} $seed $parallelization &
        wait
        time mpirun -np $n --oversubscribe ./a.out $L 0.4516 2.5 ${R_list[idx]}  $maxT ${Db_list[idx]} $seed $parallelization &
        wait
        time mpirun -np $n --oversubscribe ./a.out $L 0.45 2.5 ${R_list[idx]}  $maxT ${Db_list[idx]} $seed $parallelization &
        wait
        time mpirun -np $n --oversubscribe ./a.out $L 0.46436 1.5 ${R_list[idx]}  $maxT ${Db_list[idx]} $seed $parallelization &
        wait
        time mpirun -np $n --oversubscribe ./a.out $L 0.49598 1 ${R_list[idx]}  $maxT ${Db_list[idx]} $seed $parallelization &
        wait
done

