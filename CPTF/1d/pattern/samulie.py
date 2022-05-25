import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
mark_inset)
from mpl_toolkits.axes_grid1 import make_axes_locatable

def DrawPattern(ax, patternA, patternB):
    fs = 15 #fontsize
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

    im = ax.pcolormesh(x_range,y_range,patternB,cmap=cmap_wis,shading='auto',alpha=1)
    im = ax.pcolormesh(x_range,y_range,patternA,cmap=cmap_cool,shading='auto',alpha=1)

    #dividerB = make_axes_locatable(ax)
    #caxB = dividerB.append_axes("right",size ="5%", pad=0.05)
    #cbarB = plt.colorbar(im,cax=caxB)
    #caxB.tick_params(labelsize = fs)


