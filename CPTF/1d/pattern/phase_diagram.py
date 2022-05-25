import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
theta = [1.0,1.1,1.2,1.5,2.5]
theta2 = [1.0,1.1,1.5,1.7,2.5]

pc= [0.49598, 0.48654,0.47887,0.46436,0.4516]
pc2= [0.51065, 0.50105, 0.47300, 0.46483, 0.45160]

x = [1,1,1]
p = [0.57, 0.53, 0.49]
size=17


ax=plt.gca()
line0, = ax.plot(theta,pc,c='k',marker='o',ms=8,lw=0,ls=':',mfc='none',label = r'$D_b=0.5$')
#line1, = ax.plot(theta2,pc2,c='g',marker='s',ms=8,lw=0,mfc='none',label=r'$D_b=1$')
#line2, = ax.plot(x,p,c='r',marker='x',ms=8,lw=2,ls=':',mfc='none')
ax.tick_params(labelsize=int(size))
ax.set_xlabel(r'$\theta$',fontsize=int(size))
ax.set_ylabel(r"$P_c$",fontsize=int(size))
ax.set_xlim(0.8,2.6)
ax.set_ylim(0.44,0.52)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=int(size))
#ax.annotate("ferromagnetic",xy=(0.5,0.65),fontsize=int(size))
#ax.annotate("paramagnetic",xy=(1.2,1.6),fontsize=int(size))


#plt.scatter([1.96],[0.609],c='b',marker='^',s=200)
#ax.annotate('tricritical',xy=(2,0.61),fontsize=int(size))
plt.tight_layout()
#plt.savefig('phase.pdf')
plt.show()
