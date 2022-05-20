import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as mpe
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
mark_inset)

t, rho_B, surv = np.loadtxt("./fig4/L32.txt", usecols=(0,1,2), unpack=True)
t2, rho_B2, surv2 = np.loadtxt("./fig4/L64.txt", usecols=(0,1,2), unpack=True)
t3, rho_B3, surv3 = np.loadtxt("./fig4/L128.txt", usecols=(0,1,2), unpack=True)
t4, rho_B4, surv4 = np.loadtxt("./fig4/L256.txt", usecols=(0,1,2), unpack=True)
t5, rho_B5, surv5 = np.loadtxt("./fig4/L512.txt", usecols=(0,1,2), unpack=True)
#t6, rho_B6, surv6 = np.loadtxt("./fig4/L1024.txt", usecols=(0,1,2), unpack=True)
#t7, rho_B7, surv7= np.loadtxt("./fig4/L2048.txt", usecols=(0,1,2), unpack=True)
#t8, rho_B8, surv8= np.loadtxt("./fig4/L4096.txt", usecols=(0,1,2), unpack=True)

#t = t[1:]
#t2 = t2[1:]
#t3 = t3[1:]
#t4 = t4[1:]
#t5 = t5[1:]
#t6 = t6[1:]
#t7 = t7[1:]
#t8 = t8[1:]

#rho_B = rho_B[:-1]
#rho_B2 = rho_B2[:-1]
#rho_B3 = rho_B3[:-1]
#rho_B4 = rho_B4[:-1]
#rho_B5 = rho_B5[:-1]
#rho_B6 = rho_B6[:-1]
#rho_B7 = rho_B7[:-1]
#rho_B8 = rho_B8[:-1]

#surv = surv[:-1]
#surv2 = surv2[:-1]
#surv3 = surv3[:-1]
#surv4 = surv4[:-1]
#surv5 = surv5[:-1]
#surv6 = surv6[:-1]
#surv7 = surv7[:-1]
#surv8 = surv8[:-1]

#fig, (ax0,ax1) = plt.subplots(1,2,sharey=True)
fig = plt.figure()
ax0 = plt.subplot(1,2,1)
ax1 = plt.subplot(1,2,2)

#ax1 = plt.axes([0.7,0.7,2,2])
#ip = InsetPosition(ax0,[0.7,0.7,0.25,0.25])
#ax1.set_axes_locator(ip)
#mark_inset(ax0, ax1, loc1=2, loc2=4, fc='none', ec='0.5')

beta_nu_pall = 0.46
beta_nu_perp = 0.5

rho_B_sc = rho_B*t**(beta_nu_pall);
rho_B2_sc = rho_B2*t2**(beta_nu_pall);
rho_B3_sc = rho_B3*t3**(beta_nu_pall);
rho_B4_sc = rho_B4*t4**(beta_nu_pall);
rho_B5_sc = rho_B5*t5**(beta_nu_pall);
#rho_B6_sc = rho_B6*t6**(beta_nu_pall);
#rho_B7_sc = rho_B7*t7**(beta_nu_pall);
#rho_B8_sc = rho_B8*t8**(beta_nu_pall);

z_s = 3
z = 2

t_sc = t/(32**z)
t2_sc = t2/(64**z)
t3_sc = t3/(128**z)
t4_sc = t4/(256**z)
t5_sc = t5/(512**z)
#t6_sc = t6/(1024**z)
#t7_sc = t7/(1024**z)
#t8_sc = t8/(2048**z)

#t2 = t2/(128**z/np.log(128)); 
#t3 = t3/(256**z/np.log(256)); 
#t4 = t4/(512**z/np.log(512)); 
#t5 = t5/(1024**z/np.log(1024)); 
#t6 = t6/(2048**z/np.log(2048)); 


fs = 15

ax0.plot(t,rho_B,lw=0,ms=8,marker='d',mfc='none',label='L=32',c='C0')
ax0.plot(t2,rho_B2,lw=0,ms=8,marker='d',mfc='none',label='L=64',c='C1')
ax0.plot(t3,rho_B3,lw=0,ms=8,marker='d',mfc='none',label='L=128',c='C2')
ax0.plot(t4,rho_B4,lw=0,ms=8,marker='d',mfc='none',label='L=256',c='C3')
ax0.plot(t5,rho_B5,lw=0,ms=8,marker='d',mfc='none',label='L=512',c='C4')
#ax0.plot(t6,rho_B6,lw=0,ms=8,marker='d',mfc='none',label='L=1024',c='C5')
#ax0.plot(t7,rho_B7,lw=0,ms=8,marker='d',mfc='none',label='L=2048',c='C6')
#ax0.plot(t8,rho_B8,lw=0,ms=8,marker='d',mfc='none',label='L=4096',c='C7')

ax0.plot(t,surv,lw=3,ms=1,marker='d',mfc='none',c='C0')
ax0.plot(t2,surv2,lw=3,ms=0,marker='d',mfc='none',c='C1')
ax0.plot(t3,surv3,lw=3,ms=0,marker='d',mfc='none',c='C2')
ax0.plot(t4,surv4,lw=3,ms=0,marker='d',mfc='none',c='C3')
ax0.plot(t5,surv5,lw=3,ms=0,marker='d',mfc='none',c='C4')
#ax0.plot(t6,surv6,lw=3,ms=0,marker='d',mfc='none',c='C5')
#ax0.plot(t7,surv7,lw=3,ms=0,marker='d',mfc='none',c='C6')
#ax0.plot(t8,surv8,lw=3,ms=0,marker='d',mfc='none',c='C7')


ax1.plot(t_sc,rho_B_sc,lw=0,ms=8,marker='d',mfc='none',c='C0')
ax1.plot(t2_sc,rho_B2_sc,lw=0,ms=8,marker='d',mfc='none',c='C1')
ax1.plot(t3_sc,rho_B3_sc,lw=0,ms=8,marker='d',mfc='none',c='C2',label='L=128')
ax1.plot(t4_sc,rho_B4_sc,lw=0,ms=8,marker='d',mfc='none',c='C3',label='L=256')
ax1.plot(t5_sc,rho_B5_sc,lw=0,ms=8,marker='d',mfc='none',c='C4',label='L=512')
#ax1.plot(t6_sc,rho_B6_sc,lw=0,ms=8,marker='d',mfc='none',c='C5',label='L=1024')
#ax1.plot(t7_sc,rho_B7_sc,lw=0,ms=8,marker='d',mfc='none',c='C6',label='L=2048')
#ax1.plot(t8_sc,rho_B8_sc,lw=0,ms=8,marker='d',mfc='none',c='C7',label='L=4096')



ax0.tick_params(labelsize=fs)
ax1.tick_params(labelsize=fs)

ax0.set_xlim(3.8,14000000)
ax0.set_ylim(0.000007,2.7)

ax1.set_xlim(1e-5,20)
ax1.set_ylim(0.0008,8)

ax0.set_xscale('log'); ax0.set_yscale('log');
ax1.set_xscale('log'); ax1.set_yscale('log');

ax0.set_xlabel(r'$t$', fontsize=fs)
ax0.set_ylabel(r'$\langle \rho_B(t) \rangle, \langle P_{surv} \rangle$', fontsize=fs)

ax1.set_xlabel(r'$t/L^{(z=2.0)}$', fontsize=fs)
ax1.set_ylabel(r'$\langle \rho_B \rangle t^{\beta/\nu_{\parallel} = 0.46}$', fontsize=fs)


plt.legend(fontsize=fs)
plt.tight_layout()
plt.show()

