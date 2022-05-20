import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as mpe
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
mark_inset)
from mpl_toolkits.axes_grid1 import make_axes_locatable


#patternA = np.loadtxt("./fig5/A_subdiff_single.txt")
#patternB = np.loadtxt("./fig5/B_subdiff_single.txt")

patternA = np.loadtxt("./pattern/rho676/A.txt")
patternB = np.loadtxt("./pattern/rho676/B.txt")

fs = 15 #fontsize

#rho_B = 400*rho_B
#rho_B = rho_B/5

#temp_avg = np.average(rho_B[19000:20000])
#print(temp_avg)

#fig, (ax0,ax1) = plt.subplots(1,2,sharey=True)
fig = plt.figure(dpi=150)
ax0 = plt.subplot(1,2,1)
#ax0_A1 = plt.subplot(2,3,2)
#ax0_A2 = plt.subplot(2,3,3)

ax1 = plt.subplot(1,2,2)
#ax1_B1 = plt.subplot(2,3,5)
#ax1_B2 = plt.subplot(2,3,6)


#Set 0 value to the black regadless of colormap
#patternA = np.ma.masked_where(patternA<6.765,patternA)
patternA = np.ma.masked_where(patternA<0.1,patternA)
patternB = np.ma.masked_where(patternB<0.1,patternB)

cmap_cool = mpl.cm.get_cmap("inferno").copy()
cmap_wis = mpl.cm.get_cmap("hsv").copy()
cmap_cool.set_bad(color='white')
cmap_wis.set_bad(color='black')

im = ax0.imshow(patternA,aspect='auto',cmap=cmap_cool)
#im = ax0_A1.imshow(patternA,aspect='auto',cmap=cmap_cool)
#im = ax0_A2.imshow(patternA,aspect='auto',cmap=cmap_cool)

#dividerA = make_axes_locatable(ax0_A2)
dividerA = make_axes_locatable(ax0)
caxA = dividerA.append_axes("right",size ="5%", pad=0.05)
cbarA = plt.colorbar(im,cax=caxA)
caxA.tick_params(labelsize = fs)

im = ax1.imshow(patternB,aspect='auto',cmap=cmap_wis)
#im = ax1_B1.imshow(patternB,aspect='auto',cmap=cmap_wis)
#im = ax1_B2.imshow(patternB,aspect='auto',cmap=cmap_wis)

#dividerB = make_axes_locatable(ax1_B2)
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

#ax0_A1.tick_params(labelsize=fs)
#ax0_A2.tick_params(labelsize=fs)
#ax1_B1.tick_params(labelsize=fs)
#ax1_B2.tick_params(labelsize=fs)

ax0.set_xlim(216,296)
ax0.set_ylim(250,0)
ax1.set_xlim(216,296)
ax1.set_ylim(250,0)

#ax0_A1.set_xlim(216,296)
#ax0_A1.set_ylim(300,150)
#ax1_B1.set_xlim(216,296)
#ax1_B1.set_ylim(300,150)
#
#ax0_A2.set_xlim(216,296)
#ax0_A2.set_ylim(150,0)
#ax1_B2.set_xlim(216,296)
#ax1_B2.set_ylim(150,0)



plt.tight_layout()
plt.show()

