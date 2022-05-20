import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
mark_inset)



#t, rho_B,surv = np.loadtxt("./fig4/L64_rho4.75.txt", usecols=(0,1,2), unpack=True)
t2, rho_B2, surv2 = np.loadtxt("./fig4/L128_rho4.75.txt", usecols=(0,1,2), unpack=True)
t3, rho_B3, surv3 = np.loadtxt("./fig4/L256_rho4.75.txt", usecols=(0,1,2), unpack=True)
t4, rho_B4, surv4 = np.loadtxt("./fig4/L512_rho4.75.txt", usecols=(0,1,2), unpack=True)
t5, rho_B5, surv5 = np.loadtxt("./fig4/L1024_rho4.75.txt", usecols=(0,1,2), unpack=True)
t6, rho_B6, surv6 = np.loadtxt("./fig4/L2048_rho4.75.txt", usecols=(0,1,2), unpack=True)
#t, rho_B = np.loadtxt("./temporal_fluctuation/L400_40p_1e+8Tmax.txt", usecols=(0,1), unpack=True)

#rho_B = 400*rho_B
#rho_B = rho_B/4.75

#temp_avg = np.average(rho_B[19000:20000])
#print(temp_avg)

fig, ax1 = plt.subplots()

#ax2 = plt.axes([0.7,0.7,2,2])
#ip = InsetPosition(ax1,[0.7,0.7,0.25,0.25])
#ax2.set_axes_locator(ip)
#mark_inset(ax1, ax2, loc1=2, loc2=4, fc='none', ec='0.5')

#beta_nu = 0.336
#rho_B2 = rho_B2*t2**(beta_nu);
#rho_B3 = rho_B3*t3**(beta_nu);
#rho_B4 = rho_B4*t4**(beta_nu);
#rho_B5 = rho_B5*t5**(beta_nu);
#rho_B6 = rho_B6*t6**(beta_nu);
#
#z = 2
#t2 = t2/(128**z); 
#t3 = t3/(256**z); 
#t4 = t4/(512**z); 
#t5 = t5/(1024**z); 
#t6 = t6/(2048**z); 
#

fs = 15
ax1.set_xlabel(r'$t$', fontsize=fs)
#ax1.set_ylabel(r'$\rho_B(t) t^{\beta/\nu}$', fontsize=fs)
ax1.set_ylabel(r'$\rho_B(t), P_{surv}(t)$', fontsize=fs)
ax1.tick_params(labelsize=fs)
ax1.set_xscale('log')
ax1.set_yscale('log')
#ax1.set_xlim(1,100000)
#ax1.set_ylim(1,1000)
#ax1.plot(t,rho_B,lw=0,ms=8,marker='o',mfc='none')
ax1.plot(t2,rho_B2,lw=0,ms=8,marker='o',mfc='none',label='L=128')
ax1.plot(t3,rho_B3,lw=0,ms=8,marker='o',mfc='none',label='L=256')
ax1.plot(t4,rho_B4,lw=0,ms=8,marker='o',mfc='none',label='L=512')
ax1.plot(t5,rho_B5,lw=0,ms=8,marker='o',mfc='none',label='L=1024')
ax1.plot(t6,rho_B6,lw=0,ms=8,marker='o',mfc='none',label='L=2048')
#ax1.plot(t,surv,lw=3,ms=1,marker='o',mfc='none')
ax1.plot(t2,surv2,lw=3,ms=1,marker='o',mfc='none')
ax1.plot(t3,surv3,lw=3,ms=1,marker='o',mfc='none')
ax1.plot(t4,surv4,lw=3,ms=1,marker='o',mfc='none')
ax1.plot(t5,surv5,lw=3,ms=1,marker='o',mfc='none')
ax1.plot(t6,surv6,lw=3,ms=1,marker='o',mfc='none')

#ax2.plot(t[100000:],rho_B[100000:],lw=0,ms=2,marker='o',mfc='none')
#ax2.tick_params(which='minor',width = 3)
#ax2.set_xscale('log')
#ax2.set_yscale('log')

plt.legend()
plt.tight_layout()
plt.show()

