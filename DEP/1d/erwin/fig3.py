import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as mpe
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
mark_inset)


#t, rho_B, surv, rho_Bsite, radius = np.loadtxt("./fig3/L32.txt", usecols=(0,1,2,3,4), unpack=True)
t2, rho_B2, surv2, rho_Bsite2, radius2 = np.loadtxt("./fig3/L64.txt", usecols=(0,1,2,3,4), unpack=True)
t3, rho_B3, surv3, rho_Bsite3, radius3 = np.loadtxt("./fig3/L128.txt", usecols=(0,1,2,3,4), unpack=True)
t4, rho_B4, surv4, rho_Bsite4, radius4 = np.loadtxt("./fig3/L256.txt", usecols=(0,1,2,3,4), unpack=True)
t5, rho_B5, surv5, rho_Bsite5, radius5 = np.loadtxt("./fig3/L512.txt", usecols=(0,1,2,3,4), unpack=True)
t6, rho_B6, surv6, rho_Bsite6, radius6 = np.loadtxt("./fig3/L1024.txt", usecols=(0,1,2,3,4), unpack=True)
#t7, rho_B7, surv7, rho_Bsite7, radius7 = np.loadtxt("./fig3/L2048.txt", usecols=(0,1,2,3,4), unpack=True)
#t6, rho_B6, surv6 = np.loadtxt("./fig3/L2048.txt", usecols=(0,1,2), unpack=True)
#t, rho_B = np.loadtxt("./temporal_fluctuation/L400_40p_1e+8Tmax.txt", usecols=(0,1), unpack=True)

#rho_B = 400*rho_B
#rho_B = rho_B/5

#temp_avg = np.average(rho_B[19000:20000])
#print(temp_avg)

#fig, (ax0,ax1) = plt.subplots(1,2,sharey=True)
fig = plt.figure()
ax0 = plt.subplot(2,2,1)
ax1 = plt.subplot(2,2,2)
ax2 = plt.subplot(2,2,3)
ax3 = plt.subplot(2,2,4)

#ax1 = plt.axes([0.7,0.7,2,2])
#ip = InsetPosition(ax0,[0.7,0.7,0.25,0.25])
#ax1.set_axes_locator(ip)
#mark_inset(ax0, ax1, loc1=2, loc2=4, fc='none', ec='0.5')

beta_nu_pall = 0.087
beta_nu_perp = 0.5
#rho_B2 = rho_B2*t2**(beta_nu_perp);
#rho_B3 = rho_B3*t3**(beta_nu_perp);
#rho_B4 = rho_B4*t4**(beta_nu_perp);
#rho_B5 = rho_B5*t5**(beta_nu_perp);
#rho_B6 = rho_B6*t6**(beta_nu_perp);


#radius = radius/surv
radius2 = radius2/surv2
radius3 = radius3/surv3
radius4 = radius4/surv4
radius5 = radius5/surv5
radius6 = radius6/surv6
#radius7 = radius7/surv7

z_s = 3
#radius_sc = radius[1:]/t[1:]**(2/z_s)
radius2_sc = radius2/t2**(2/z_s)
radius3_sc = radius3/t3**(2/z_s)
radius4_sc = radius4/t4**(2/z_s)
radius5_sc = radius5/t5**(2/z_s)
radius6_sc = radius6/t6**(2/z_s)
#radius7_sc = radius7[1:]/t7[1:]**(2/z_s)


z = 2

#t = t/(32**z)
#t2 = t2/(64**z); 
#t3 = t3/(128**z); 
#t4 = t4/(256**z); 
#t5 = t5/(512**z); 
#t6 = t6/(1024**z); 
#t7 = t7/(2048**z); 


#t2 = t2/(128**z/np.log(128)); 
#t3 = t3/(256**z/np.log(256)); 
#t4 = t4/(512**z/np.log(512)); 
#t5 = t5/(1024**z/np.log(1024)); 
#t6 = t6/(2048**z/np.log(2048)); 


fs = 15
#ax0.set_xlabel(r'$t$', fontsize=fs)
#ax0.set_ylabel(r'$\rho_B(t), P_{surv}(t)$', fontsize=fs)
#ax0.set_xlabel(r'$t$', fontsize=fs)
#ax0.set_ylabel(r'$\rho_B(t)t^{\beta/\nu = 0.087}, P_{surv}(t)$', fontsize=fs)
#ax0.set_ylabel(r'$\rho_B(t)L^{\beta/\nu_{\perp} = 0.5}/P_{surv}(t)$', fontsize=fs)
ax0.tick_params(labelsize=fs)
ax0.set_xscale('log')
ax0.set_yscale('log')
#ax0.set_xlim(1,100000)
#ax0.set_ylim(1,1000)
#ax0.plot(t,rho_B,lw=0,ms=8,marker='o',mfc='none')
#ax0.plot(t2,rho_B2,lw=0,ms=8,marker='o',mfc='none',label='L=128',c='r')
#ax0.plot(t3,rho_B3,lw=0,ms=8,marker='o',mfc='none',label='L=256',c='k')
#ax0.plot(t4,rho_B4,lw=0,ms=8,marker='o',mfc='none',label='L=512',c='y')
#ax0.plot(t5,rho_B5,lw=0,ms=8,marker='o',mfc='none',label='L=1024',c='g')
#ax0.plot(t6,rho_B6,lw=0,ms=8,marker='o',mfc='none',label='L=2048',c='b')

#rho_B2 = (rho_B2/surv2)*128**beta_nu_perp;
#rho_B3 = (rho_B3/surv3)*256**beta_nu_perp;
#rho_B4 = (rho_B4/surv4)*512**beta_nu_perp;
#rho_B5 = (rho_B5/surv5)*1024**beta_nu_perp;
#rho_B6 = rho_B6/surv6;
#ax0.plot(t,rho_B,lw=0,ms=8,marker='o',mfc='none',label='L=32',c='magenta')
ax0.plot(t2,rho_B2,lw=0,ms=8,marker='o',mfc='none',label='L=64',c='C1')
ax0.plot(t3,rho_B3,lw=0,ms=8,marker='o',mfc='none',label='L=128',c='C2')
ax0.plot(t4,rho_B4,lw=0,ms=8,marker='o',mfc='none',label='L=256',c='C3')
ax0.plot(t5,rho_B5,lw=0,ms=8,marker='o',mfc='none',label='L=512',c='C4')
ax0.plot(t6,rho_B6,lw=0,ms=8,marker='o',mfc='none',label='L=1024',c='C5')
#ax0.plot(t7,rho_B7,lw=0,ms=8,marker='o',mfc='none',label='L=2048',c='C6')


outline = [mpe.Stroke(linewidth=6,foreground='k'), mpe.Stroke(foreground='w',alpha=0.1),mpe.Normal()]
#ax1.set_xlabel(r'$t$', fontsize=fs)
ax1.tick_params(labelsize=fs)
#ax1.plot(t,surv,lw=3,ms=1,marker='o',mfc='none',c='magenta')
ax1.plot(t2,surv2,lw=3,ms=0,marker='o',mfc='none',c='C1')
ax1.plot(t3,surv3,lw=3,ms=0,marker='o',mfc='none',c='C2')
ax1.plot(t4,surv4,lw=3,ms=0,marker='o',mfc='none',c='C3')
ax1.plot(t5,surv5,lw=3,ms=0,marker='o',mfc='none',c='C4')
ax1.plot(t6,surv6,lw=3,ms=0,marker='o',mfc='none',c='C5')
#ax1.plot(t7,surv7,lw=3,ms=0,marker='o',mfc='none',c='C6')

#ax1.plot(t[100000:],rho_B[100000:],lw=0,ms=2,marker='o',mfc='none')
#ax1.tick_params(which='minor',width = 3)
#ax1.tick_params(labelsize=fs)
ax1.set_xscale('log')
ax1.set_yscale('log')

#ax2.set_xlabel(r'$t$', fontsize=fs)
ax2.tick_params(labelsize=fs)
#ax2.plot(t,rho_Bsite,lw=0,ms=8,marker='o',mfc='none',c='magenta')
#ax2.plot(t2,rho_Bsite2,lw=0,ms=8,marker='o',mfc='none',c='red')
#ax2.plot(t3,rho_Bsite3,lw=0,ms=8,marker='o',mfc='none',c='k')
#ax2.plot(t4,rho_Bsite4,lw=0,ms=8,marker='o',mfc='none',c='y')
#ax2.plot(t5,rho_Bsite5,lw=0,ms=8,marker='o',mfc='none',c='g')
#ax2.plot(t6,rho_Bsite6,lw=0,ms=8,marker='o',mfc='none',c='b')

#ax2.plot(t,rho_Bsite,lw=0,ms=8,marker='o',mfc='none',c='C1')
ax2.plot(t2,rho_Bsite2,lw=0,ms=8,marker='o',mfc='none',c='C1')
ax2.plot(t3,rho_Bsite3,lw=0,ms=8,marker='o',mfc='none',c='C2')
ax2.plot(t4,rho_Bsite4,lw=0,ms=8,marker='o',mfc='none',c='C3')
ax2.plot(t5,rho_Bsite5,lw=0,ms=8,marker='o',mfc='none',c='C4')
ax2.plot(t6,rho_Bsite6,lw=0,ms=8,marker='o',mfc='none',c='C5')
#ax2.plot(t7,rho_Bsite7,lw=0,ms=8,marker='o',mfc='none',c='C6')


ax2.set_xscale('log')
ax2.set_yscale('log')

#t = t[1:];
#t2 = t2[1:]; 
#t3 = t3[1:]; 
#t4 = t4[1:]; 
#t5 = t5[1:]; 
#t6 = t6[1:]; 
#
#t_sc = t/(32**z)
#t2_sc = t2/(64**z); 
#t3_sc = t3/(128**z); 
#t4_sc = t4/(256**z); 
#t5_sc = t5/(512**z); 
#t6_sc = t6/(1024**z); 




#ax3.set_xlabel(r'$t$', fontsize=fs)
ax3.tick_params(labelsize=fs)
#ax3.plot(t,radius,lw=0,ms=8,marker='o',mfc='none',c='magenta')
ax3.plot(t2,radius2,lw=0,ms=8,marker='o',mfc='none',c='C1')
ax3.plot(t3,radius3,lw=0,ms=8,marker='o',mfc='none',c='C2')
ax3.plot(t4,radius4,lw=0,ms=8,marker='o',mfc='none',c='C3')
ax3.plot(t5,radius5,lw=0,ms=8,marker='o',mfc='none',c='C4')
ax3.plot(t6,radius6,lw=0,ms=8,marker='o',mfc='none',c='C5')
#ax3.plot(t7[10:],radius7[10:],lw=0,ms=8,marker='o',mfc='none',c='C6')

#ax3.plot(t2[1:],radius2_sc,lw=0,ms=8,marker='o',mfc='none',c='C1')
#ax3.plot(t3[1:],radius3_sc,lw=0,ms=8,marker='o',mfc='none',c='C2')
#ax3.plot(t4[1:],radius4_sc,lw=0,ms=8,marker='o',mfc='none',c='C3')
#ax3.plot(t5[1:],radius5_sc,lw=0,ms=8,marker='o',mfc='none',c='C4')
#ax3.plot(t6[1:],radius6_sc,lw=0,ms=8,marker='o',mfc='none',c='C5')
#ax3.plot(t7[1:],radius7_sc,lw=0,ms=8,marker='o',mfc='none',c='C6')



ax3.set_xscale('log')
ax3.set_yscale('log')


ax0.set_ylabel(r'$\rho_B(t)$', fontsize=fs)
ax1.set_ylabel(r'$P_{surv}(t)$', fontsize=fs)
ax2.set_ylabel(r'$\rho_{B site}(t)$', fontsize=fs)
#ax3.set_ylabel(r'$\langle R^2 \rangle/t^{(2/z = 2/2)}$', fontsize=fs)
ax3.set_ylabel(r'$\langle R^2 \rangle$', fontsize=fs)

ax0.set_xlabel(r'$t/L^{(z=2.0)}$', fontsize=fs)
ax1.set_xlabel(r'$t/L^{(z=2.0)}$', fontsize=fs)
ax2.set_xlabel(r'$t/L^{(z=1.75)}$', fontsize=fs)
#ax3.set_xlabel(r'$t/L^{(z=2.0)}$', fontsize=fs)
ax3.set_xlabel(r'$t$', fontsize=fs)


log_tF = np.log(t6[20:35])
log_tI = np.log(t2[1:10])
log_RF = np.log(radius6[20:35])
log_RI = np.log(radius2[1:10])

slopeF, interceptF = np.polyfit(log_tF,log_RF,1)
slopeI, interceptI = np.polyfit(log_tI,log_RI,1)
#y_trend = slope*log_t + intercept

#slopeF = 2/3
#slopeI = 2/2
#
y_trendF = (t4**slopeF)*np.exp(interceptF+0.6)
y_trendI = (t2**slopeI)*np.exp(interceptI-0.8)
y_trendI = (t2**slopeI)*np.exp(interceptF-0.9)
#ax3.plot(t5[1:10000], np.exp(y_trend), lw=4, ls=':',ms=0)

ax3.plot(t4[20:35], y_trendF[20:35], lw=2, ls='dashed',ms=0,c='k')
ax3.plot(t2[1:10], y_trendI[1:10], lw=2, ls='dashed',ms=0,c='k')

text5 = ("z="+str(round(2/slopeF,2)))
text4 = ("z="+str(round(2/slopeI,2)))
ax3.text(1000,10000,text5,fontsize = fs,wrap=True)
ax3.text(10,300,text4,fontsize = fs,wrap=True)

#fig.legend(bbox_to_anchor=(0.5,1),loc='upper center',ncol=6,fontsize=fs)
fig.legend(bbox_to_anchor=(0.5,1),loc='upper center',ncol=1,fontsize=fs)

plt.tight_layout()
plt.show()

