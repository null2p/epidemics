#!/bin/bash

time mpirun -np 30 --oversubscribe a.out 32 > fig4/L32.txt
wait
time mpirun -np 30 --oversubscribe a.out 64 > fig4/L64.txt
wait
time mpirun -np 30 --oversubscribe a.out 128 > fig4/L128.txt
wait
time mpirun -np 30 --oversubscribe a.out 256 > fig4/L256.txt
wait
time mpirun -np 30 --oversubscribe a.out 512 > fig4/L512.txt
wait
time mpirun -np 30 --oversubscribe a.out 1024 > fig4/L1024.txt
wait
time mpirun -np 30 --oversubscribe a.out 2048 > fig4/L2048.txt
wait
