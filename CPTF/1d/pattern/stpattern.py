import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as mpe
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
mark_inset)
from mpl_toolkits.axes_grid1 import make_axes_locatable
from trendline import InsertTrendLine

rho_c = 10
folder = "rho7.5_avoidDA1DB1"
rho_c = 7.5
folder = "rho7.5_antibiasDA1DB0.75"
#rho_c = 6.76
#fold = "rho676/t8170"
#rho_c = 6.765
#fold = "rho6765/t4686"
#rho_c = 6.77
#fold = "rho677/t3317"

L = 128
theta = 1.5
p = 0.473
Db = 1
ens = 7

patternA = np.loadtxt("./seed_single/theta"+str(theta)+"/Aparticle_L"+str(L)+"T256ens"+str(ens)+"p"+str(p)+"theta"+str(theta)+"R0.5Db"+str(Db)+".txt")
patternB = np.loadtxt("./seed_single/theta"+str(theta)+"/Bparticle_L"+str(L)+"T256ens"+str(ens)+"p"+str(p)+"theta"+str(theta)+"R0.5Db"+str(Db)+".txt")
A_max = np.max(patternA)

t,rho_A, rho_B, surv,rho_Asite, rho_Bsite,radiusA ,radiusB,NAcluster, NBcluster = np.loadtxt("./seed_single/theta"+str(theta)+"/Scailing_L"+str(L)+"T256ens"+str(ens)+"p"+str(p)+"theta"+str(theta)+"R0.5Db"+str(Db)+".txt", usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)


fs = 15 #fontsize

'''
#############################################################################################
###################################### A particle ###########################################
fig = plt.figure(figsize=(15,15))
ax0 = plt.subplot(2,2,1)
ax1 = plt.subplot(2,2,2)
ax2 = plt.subplot(2,2,3)
ax3 = plt.subplot(2,2,4)

ax0.tick_params(labelsize=fs)
ax0.set_xscale('log')
ax0.set_yscale('log')
ax0.plot(t,rho_A,lw=0,ms=8,marker='o',mfc='none',label='L=32',c='C0')


#outline = [mpe.Stroke(linewidth=6,foreground='k'), mpe.Stroke(foreground='w',alpha=0.1),mpe.Normal()]
ax1.tick_params(labelsize=fs)
ax1.plot(t,NAcluster,lw=0,ms=8,marker='o',mfc='none',c='C0')

ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.tick_params(labelsize=fs)
ax2.plot(t,rho_Asite,lw=0,ms=8,marker='o',mfc='none',c='C0')

ax2.set_xscale('log')
ax2.set_yscale('log')


ax3.tick_params(labelsize=fs)
ax3.plot(t,radiusA,lw=0,ms=8,marker='o',mfc='none',c='C0')


ax3.set_xscale('log')
ax3.set_yscale('log')


ax0.set_ylabel(r'$\rho_A(t)$', fontsize=fs)
ax1.set_ylabel(r'$\langle N_A(t)\rangle$', fontsize=fs)
ax2.set_ylabel(r'$\rho_{A site}(t)$', fontsize=fs)
#ax3.set_ylabel(r'$\langle R^2 \rangle/t^{(2/z = 2/2)}$', fontsize=fs)
ax3.set_ylabel(r'$\langle R_A^2 \rangle$', fontsize=fs)

ax0.set_xlabel(r'$t$', fontsize=fs)
ax1.set_xlabel(r'$t$', fontsize=fs)
ax2.set_xlabel(r'$t$', fontsize=fs)
ax3.set_xlabel(r'$t$', fontsize=fs)

ti = 1000; tf = 7000;
#ti = 1; tf = np.size(t)-1;

#InsertTrendLine (ax0, ti, tf,t, rho_A, 'theta')
#InsertTrendLine (ax1, ti, tf,t, NAcluster, 'slope')
#InsertTrendLine (ax2, ti, tf,t, rho_Asite, 'slope')
#InsertTrendLine (ax3, ti, tf,t, radiusA, 'z')

#fig.legend(bbox_to_anchor=(0.5,1),loc='upper center',ncol=6,fontsize=fs)

#plt.show()



###################################### A particle ###########################################
#############################################################################################




#############################################################################################
###################################### B particle ###########################################
fig = plt.figure(figsize=(15,15))
ax0 = plt.subplot(2,2,1)
ax1 = plt.subplot(2,2,2)
ax2 = plt.subplot(2,2,3)
ax3 = plt.subplot(2,2,4)

ax0.tick_params(labelsize=fs)
ax0.set_xscale('log')
ax0.set_yscale('log')
ax0.plot(t,rho_B,lw=0,ms=8,marker='o',mfc='none',label='L=32',c='C0')


#outline = [mpe.Stroke(linewidth=6,foreground='k'), mpe.Stroke(foreground='w',alpha=0.1),mpe.Normal()]
ax1.tick_params(labelsize=fs)
ax1.plot(t,NBcluster,lw=0,ms=8,marker='o',mfc='none',c='C0')

ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.tick_params(labelsize=fs)
ax2.plot(t,rho_Bsite,lw=0,ms=8,marker='o',mfc='none',c='C0')

ax2.set_xscale('log')
ax2.set_yscale('log')


ax3.tick_params(labelsize=fs)
ax3.plot(t,radiusB,lw=0,ms=8,marker='o',mfc='none',c='C0')


ax3.set_xscale('log')
ax3.set_yscale('log')


ax0.set_ylabel(r'$\rho_B(t)$', fontsize=fs)
ax1.set_ylabel(r'$\langle N_B(t)\rangle$', fontsize=fs)
ax2.set_ylabel(r'$\rho_{B site}(t)$', fontsize=fs)
#ax3.set_ylabel(r'$\langle R^2 \rangle/t^{(2/z = 2/2)}$', fontsize=fs)
ax3.set_ylabel(r'$\langle R_B^2 \rangle$', fontsize=fs)

ax0.set_xlabel(r'$t$', fontsize=fs)
ax1.set_xlabel(r'$t$', fontsize=fs)
ax2.set_xlabel(r'$t$', fontsize=fs)
ax3.set_xlabel(r'$t$', fontsize=fs)

ti = 1000; tf = 7000
#ti = 1; tf = np.size(t)-1;

#InsertTrendLine (ax0, ti, tf,t, rho_B, 'theta')
#InsertTrendLine (ax1, ti, tf,t, NBcluster, 'slope')
#InsertTrendLine (ax2, ti, tf,t, rho_Bsite, 'slope')
#InsertTrendLine (ax3, ti, tf,t, radiusB, 'z')

#fig.legend(bbox_to_anchor=(0.5,1),loc='upper center',ncol=6,fontsize=fs)



#plt.show()



###################################### B particle ###########################################
#############################################################################################
'''



#fig = plt.figure(dpi=150)
fig = plt.figure()
ax0 = plt.subplot(1,2,1)
ax1 = plt.subplot(1,2,2)


#Set 0 value to the black regadless of colormap
#patternA = np.ma.masked_where(patternA<6.765,patternA)
patternA = np.ma.masked_where(patternA<0.1,patternA)
patternB = np.ma.masked_where(patternB<0.1,patternB)

A_max = np.max(patternA)
B_max = np.max(patternB)

cmap_cool = mpl.cm.get_cmap("brg",A_max).copy()
cmap_wis = mpl.cm.get_cmap("Wistia",B_max).copy()
cmap_cool.set_bad(color='none')
cmap_wis.set_bad(color='none')


time_max = patternA.shape[0]
x_max = patternA.shape[1]

centers = [-x_max/2,x_max/2-1,patternA.shape[0],1]

dx, = np.diff(centers[:2])/(patternA.shape[1]-1)
dy, = -np.diff(centers[2:])/(patternA.shape[0]-1)

extent = [centers[0]-dx/2, centers[1]+dx/2, centers[2]+dy, centers[3]-dy]
x_range = np.linspace(-x_max/2,+x_max/2-1,num=x_max)
y_range = np.linspace(0,patternA.shape[0],num=patternA.shape[0])


#im = ax0.pcolormesh(x_range,y_range,patternA,cmap=cmap_cool,shading='auto')
#dividerA = make_axes_locatable(ax0)
#caxA = dividerA.append_axes("right",size ="5%", pad=0.05)
#cbarA = plt.colorbar(im,cax=caxA)
#caxA.tick_params(labelsize = fs)

im = ax0.pcolormesh(x_range,y_range,patternB,cmap=cmap_wis,shading='auto',alpha=1)
im = ax0.pcolormesh(x_range,y_range,patternA,cmap=cmap_cool,shading='auto',alpha=1)
dividerB = make_axes_locatable(ax0)
caxB = dividerB.append_axes("right",size ="5%", pad=0.05)
cbarB = plt.colorbar(im,cax=caxB)
caxB.tick_params(labelsize = fs)




#im = ax1.pcolormesh(x_range,y_range,patternB,vmin=1,vmax=20,cmap=cmap_wis,shading='auto')
im = ax1.pcolormesh(x_range,y_range,patternB,cmap=cmap_wis,shading='auto',alpha=1)
im = ax1.pcolormesh(x_range,y_range,patternA,cmap=cmap_cool,shading='auto',alpha=1)
dividerB = make_axes_locatable(ax1)
caxB = dividerB.append_axes("right",size ="5%", pad=0.05)
cbarB = plt.colorbar(im,cax=caxB)
caxB.tick_params(labelsize = fs)


ax0.set_xlabel(r'$x$', fontsize=fs)
ax0.set_ylabel(r'$t$', fontsize=fs)
ax0.tick_params(labelsize=fs)

ax1.set_xlabel(r'$x$', fontsize=fs)
ax1.set_ylabel(r'$t$', fontsize=fs)
ax1.tick_params(labelsize=fs)


ax0.set_xlim(-50,50)
ax0.set_ylim(time_max,0)
#ax0.set_ylim(500,0)

#ax1.set_xlim(-128,128)
ax1.set_ylim(time_max,0)
#ax1.set_ylim(500,0)

plt.tight_layout()
plt.show()



