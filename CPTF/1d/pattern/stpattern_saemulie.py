import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as mpe
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
mark_inset)
from mpl_toolkits.axes_grid1 import make_axes_locatable
from trendline import InsertTrendLine
import glob
from samulie import DrawPattern

L = 128;  Db = 0.75 ; R = 0.5;

p = [0.45, 0.46436, 0.49598] #Db=0.5
#p = [0.4516, 0.473, 0.51065] #Db=1.0


patternAname25 = "./seed_single/theta2.5/Aparticle_L"+str(L)+"T1024ens*p"+str(p[0])+"theta2.5R"+str(R)+"Db"+str(Db)+".txt"
patternBname25 = "./seed_single/theta2.5/Bparticle_L"+str(L)+"T1024ens*p"+str(p[0])+"theta2.5R"+str(R)+"Db"+str(Db)+".txt"

patternAname15 = "./seed_single/theta1.5/Aparticle_L"+str(L)+"T1024ens*p"+str(p[1])+"theta1.5R"+str(R)+"Db"+str(Db)+".txt"
patternBname15 = "./seed_single/theta1.5/Bparticle_L"+str(L)+"T1024ens*p"+str(p[1])+"theta1.5R"+str(R)+"Db"+str(Db)+".txt"

patternAname10 = "./seed_single/theta1/Aparticle_L"+str(L)+"T1024ens*p"+str(p[2])+"theta1R"+str(R)+"Db"+str(Db)+".txt"
patternBname10 = "./seed_single/theta1/Bparticle_L"+str(L)+"T1024ens*p"+str(p[2])+"theta1R"+str(R)+"Db"+str(Db)+".txt"


Aname25 = glob.glob(patternAname25); Bname25 = glob.glob(patternBname25);
Aname15 = glob.glob(patternAname15); Bname15 = glob.glob(patternBname15);
Aname10 = glob.glob(patternAname10); Bname10 = glob.glob(patternBname10);

Aname25.sort() ;Bname25.sort();
Aname15.sort() ;Bname15.sort();
Aname10.sort() ;Bname10.sort();

print(np.size(Aname25))
print(np.size(Aname15))
print(np.size(Aname10))

for i in range(4):
    patternA25 = np.loadtxt(Aname25[i]); patternB25 = np.loadtxt(Bname25[i]);
    patternA15 = np.loadtxt(Aname15[i]); patternB15 = np.loadtxt(Bname15[i]);
    patternA10 = np.loadtxt(Aname10[i]); patternB10 = np.loadtxt(Bname10[i]);

    fig = plt.figure(dpi=150)
    ax0 = plt.subplot(1,3,1)
    ax1 = plt.subplot(1,3,2)
    ax2 = plt.subplot(1,3,3)

    fs=15

    #time_max = patternA25.shape[0]
    time_max = 1024
#ax.set_xlabel(r'$x$', fontsize=fs)
#ax.set_ylabel(r'$t$', fontsize=fs)
    #ax0.tick_params(labelsize=fs)

    DrawPattern(ax0, patternA25, patternB25)
    DrawPattern(ax1, patternA15, patternB15)
    DrawPattern(ax2, patternA10, patternB10)

    ax0.set_ylim(time_max,0)
    ax1.set_ylim(time_max,0)
    ax2.set_ylim(time_max,0)

    plt.tight_layout()
    plt.show()







