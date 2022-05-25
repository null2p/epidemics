import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stat
from math import exp
import os
import argparse
import itertools
import sys
import copy
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = 'ensemble number')
parser.add_argument('--ens',type=int,help='ENS NUMBER')
args = parser.parse_args()

'''
tau_L10000T_i1000T_f10000ens4p0.9theta1.5R0.93.txt
'''

df = pd.read_table('./seed_half/theta1.5/tau_L10000T_i1000T_f10000ens4p0.9theta1.5R0.93.txt',delimiter=' ', nrows=1)
#tau = np.loadtxt("./seed_half/theta1.5/tau_L10000T_i1000T_f10000ens"+str(args.ens)+"p0.9theta1.5R0.93.txt")

print(df)
#print(np.size(df))

'''
L = np.size(tau[0])
T_max =  np.size(tau,0)

T = np.linspace(0,T_max-1,num=T_max)
A = np.zeros(T_max)
B = np.zeros(T_max)
R = np.zeros(T_max)

for t in range(T_max) :
  for x in range(L):
    if(tau[t,x]) > 0 : 
      A[t] += 1
      #R_tmp = (L/2 - x)**2
      #if R_tmp > R[t] : R[t] = copy.deepcopy(R_tmp)
    elif (tau[t,x]<0) :
      B[t] += 1
    #if R[t] < R[t-1] : R[t] =copy.deepcopy(R[t-1])

A = A/L
B = B/L



fig = plt.figure(figsize=(10,7))
ax1 = fig.add_subplot(3, 1, 1)
ax2 = fig.add_subplot(3, 1, 2)
ax3 = fig.add_subplot(3, 1, 3)


tau = tau.transpose()
#colors = ["yellow","white","blue","red"]
colors = ["yellow","lightgray","blue"]
colormap = mpl.colors.ListedColormap(colors)
im = ax1.imshow(tau, cmap=colormap,aspect='auto')
#plt.colorbar(ticks=range(len(colors)))
#plt.colorbar(im, ax = ax1)
#ax1.set_xlabel('Time',fontsize=20)
ax1.set_ylabel('x',fontsize=20)
ax1.tick_params(labelsize = 20)
#ax1.set_xlim(0,512)
#ax1.set_ylim(0,128)

ax2.plot(T,A,label='I')
ax2.plot(T,B,label='Q')
ax2.plot(T,A+B,label='I+Q')
#ax2.set_xlabel('t',fontsize=20)
ax2.set_ylabel(r'$\rho$',fontsize=20)
ax2.tick_params(labelsize = 20)
ax2.legend(fontsize=20)
#ax2.set_xlim(0,512)
#ax2.set_ylim(0,0.3)
plt.tight_layout()
plt.show()
'''

