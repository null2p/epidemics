#!/bin/bash

time mpirun -np 30 --oversubscribe a.out 64 > fig3/L64.txt
wait
time mpirun -np 30 --oversubscribe a.out 128 > fig3/L128.txt
wait
time mpirun -np 30 --oversubscribe a.out 256 > fig3/L256.txt
wait
time mpirun -np 30 --oversubscribe a.out 512 > fig3/L512.txt
wait
time mpirun -np 30 --oversubscribe a.out 1024 > fig3/L1024.txt
wait
time mpirun -np 30 --oversubscribe a.out 2048 > fig3/L2048.txt
wait
time mpirun -np 30 --oversubscribe a.out 4096 > fig3/L4096.txt
wait
