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
import glob

parser = argparse.ArgumentParser(description = 'ens for drawing pattern')
parser.add_argument('--p',type=float,help='healing probability')
parser.add_argument('--Db',help='direction bias')
parser.add_argument('--seed',type=str,help='seed')
#parser.add_argument('--theta',type=float,help='exponent theta')
parser.add_argument('--theta',help='exponent theta')
#parser.add_argument('--ens',type=int,help='ENS NUMBER')
args = parser.parse_args()

t_min = 0; t_max = 1000;

'''
../static_pattern/seed_half/theta1/Aparticle_L10000T1000000ens64p0.49598theta1R0.5Db0.5.txt
seed_half/theta1.5/pattern_L60000T_i10000T_f100000ens0p0.6theta1.5R0.65.txt
seed_half/theta1.5/pattern_L10000T_i10000T_f100000ens0p0.6theta1.5R0.65.txt
seed_half/theta1/pattern_L100T_i0T_f10000ens0p0.56theta1R0.5.txt
seed1/theta1/pattern_L100T_i0T_f10000ens4p0.48theta1R0.5Db0.5.txt
patternA = np.loadtxt("./seed_single/theta"+str(theta)+"/Aparticle_L1024T10000ens0p"+str(p)+"theta"+str(theta)+"R0.5Db"+str(Db)+".txt")
'''
txtdirname= ' '; Aname = ' '; pngdirname = ' ' ; Bname= ' ' ;


if(args.seed == 'half') :  
	txtdirname = "./seed_half/theta"+str(args.theta)
	Aname =  "/pattern_L100T_i0T_f10000ens*p"+str(args.p)+"theta*R0.5Db"+str(args.Db)+".txt" #if Db=0.5, Db is null in file name, file name is differnt in seed_half folder...
	pngdirname = "./seed_half/theta"+str(args.theta)+"/Db"+str(args.Db)+"/r0.5/p"+str(args.p)
elif(args.seed == 'single') : 
	txtdirname = "./seed_single/theta"+str(args.theta)
	Aname =  "/Aparticle_L1024T10000ens*p"+str(args.p)+"theta*R0.5Db"+str(args.Db)+".txt" 
	Bname =  "/Bparticle_L1024T10000ens*p"+str(args.p)+"theta*R0.5Db"+str(args.Db)+".txt" 
	pngdirname = "./seed_single/theta"+str(args.theta)+"/Db"+str(args.Db)+"/r0.5/p"+str(args.p)

#Aname =  "/pattern_L100T_i0T_f10000ens"+str(args.ens)+"p"+str(args.p)+"theta"+str(args.theta)+"R0.5Db1"
#Aname =  "/pattern_L100T_i0T_f10000ens*p"+str(args.p)+"*Db0.5*"

#Aparticle_L10000T1000ens10, 13, 38, 6p0.5theta1R0.6Db0.5.txt

Aname_array = glob.glob(txtdirname+Aname)
Bname_array = glob.glob(txtdirname+Bname)
print( Aname_array )
print( Bname_array )

os.makedirs(pngdirname,exist_ok=True)


for each_Aname, each_Bname in zip(Aname_array, Bname_array) :
	patternA = np.loadtxt( each_Aname )
	patternB = np.loadtxt( each_Bname )

	L = np.size(patternA[0])
	T_max =  np.size(patternA,0)
	print(L,T_max)

	T = np.linspace(0,T_max-1,num=T_max)
	A = np.zeros(T_max)
	B = np.zeros(T_max)
	R = np.zeros(T_max)

	
	for t in range(T_max):
	  A[t] = np.sum(patternA[t])
	  B[t] = np.sum(patternB[t])

#	for t in range(T_max) :
#	  for x in range(L):
#	    if(patternA[t,x]) > 0 : 
#	      A[t] += 1
#	      R_tmp = (L/2 - x)**2
#	      if R_tmp > R[t] : R[t] = copy.deepcopy(R_tmp)
#	    elif (patternA[t,x]<0) :
#	      B[t] += 1
	#if R[t] < R[t-1] : R[t] =copy.deepcopy(R[t-1])

	#A = A/L
	#B = B/L
	R = np.sqrt(R)

	fig = plt.figure(figsize=(5,10))
	ax1 = fig.add_subplot(141)
	ax2 = fig.add_subplot((142),sharey=ax1)
	ax3 = fig.add_subplot((143),sharey=ax1,sharex=ax1)
	#ax3 = fig.add_subplot(3, 1, 3)
	ax1.set_title(r"$p = $"+str(args.p)+r"$, \theta = $"+str(args.theta)+r"$, r = 0.5, D_b = $"+str(args.Db),fontsize=20,loc='left',pad=12)
	#pattern = pattern.transpose()

	colors = ["yellow","white","blue","red"]
	colormap = mpl.colors.ListedColormap(colors)
	im = ax1.imshow(pattern, cmap=colormap,aspect='auto')
	#im = ax1.imshow(patternA,aspect='auto')
	
	#plt.colorbar(ticks=range(len(colors)))
	#plt.colorbar(im, ax = ax1)
	#ax1.set_xlabel('Time',fontsize=20)
	ax1.set_xlabel('x',fontsize=20)
	ax1.set_ylabel('t',fontsize=20)
	ax1.tick_params(labelsize = 20)
	ax1.set_ylim(t_max,t_min)
	#ax1.set_xlim(4980,5040)

	ax2.plot(A+B,T,label='I+Q',c='k')
	ax2.plot(B,T,label='Q',c='yellow')
	ax2.plot(A,T,label='I',c='b')
	#ax2.plot(T,A+B,label='I+Q')
	#ax2.plot(T,B,label='Q')
	#ax2.plot(T,A,label='I')
	#ax2.set_xlabel('t',fontsize=20)
	ax2.set_xlabel(r'$\rho$',fontsize=20)
	ax2.tick_params(labelsize = 20)
	ax2.legend(fontsize=20,loc='best')
	ax2.set_ylim(t_max,t_min)
	#ax2.set_ylim(0,1)

	im = ax3.imshow(patternB,aspect='auto')
	#plt.colorbar(ticks=range(len(colors)))
	#plt.colorbar(im, ax = ax1)
	#ax1.set_xlabel('Time',fontsize=20)
	ax3.tick_params(labelsize = 20)
	ax3.set_ylim(t_max,t_min)
	#ax3.set_xlim(4980,5040)

	ax4 = fig.add_subplot((144),projection='3d')
	X_3d = np.arange(np.size(patternA[0,4980:5040]))
	Y_3d = np.arange(np.size(patternA[0:100,0]))
	
	patternA = patternA[0:100,4980:5040]
	
	X_3d, Y_3d = np.meshgrid(X_3d,Y_3d)
	ax4.plot_surface(X_3d,Y_3d, patternA, ccount = 1, rcount =59)
	#ax4.plot_wireframe(X_3d,Y_3d, patternA, cstride=0, rstride = 1)
	#ax4.set_ylim(t_max,t_min)
	#ax4.set_xlim(4980,5040)



	'''
	ax3.plot(T,R,lw=0,marker='o',c='blue',ms=5,mfc='none')
	ax3.set_xlabel('t',fontsize=20)
	ax3.set_ylabel(r'$R_A$',fontsize=20)
	ax3.tick_params(labelsize = 20)
	ax3.set_xlim(t_min,t_max)
	#ax3.set_xscale('log')
	#ax3.set_yscale('log')
	'''

	plt.tight_layout()
	plt.show()

	#fig.savefig(pngdirname + os.path.basename(each_Aname)[:-4]+".png",dpi=200)
