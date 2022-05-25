#!/bin/bash

n="4"                      #mpi size
L="128"                     #size of lattice
p="0.49598"                  #1-p is infection rate
theta="1"                 #theta is related with the life of B particles
R_btoa="0.5"                #the probability of B -> A
newR="0.5"                #the probability of B -> A
initT="0"                   #initial time for saving state
midT="0"
finalT="256"                #maximal time for this dynamics
D_b='0.5'                   #direction bias -> branching bias
seed_type='single'                 #seed type = 'single' or 'half'

time mpirun -np $n --oversubscribe ./pattern $L $p $theta $R_btoa $newR $initT $midT $finalT $D_b $seed_type & 
#wait
