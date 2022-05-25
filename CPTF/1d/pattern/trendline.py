import numpy as np

def InsertTrendLine (ax, ti, tf,time, observable, exponent) :
    fs = 20
    #time = np.linspace(1,10000,10000)
    log_tI = np.log(time[ti:tf])
    log_RI = np.log(observable[ti:tf])

    slopeI, interceptI = np.polyfit(log_tI,log_RI,1)
    #slopeI = 2/2
    y_trendI = (time**slopeI)*np.exp(interceptI)

    ax.plot(time[ti:tf], y_trendI[ti:tf], lw=2, ls='dashed',ms=0,c='k')
    y_text = observable[tf]
    if(exponent == 'z') :
        text_z = ("z="+str(round(2/slopeI,2)))
        ax.text(time[ti],y_text,text_z,fontsize = fs,wrap=True)
    if(exponent == 'theta') :
        text_theta = (r"$\theta=$"+str(round(slopeI,2)))
        ax.text(time[ti],y_text,text_theta,fontsize = fs,wrap=True)
    if(exponent == 'delta') :
        text_delta = (r"$\delta=$"+str(round(-slopeI,2)))
        ax.text(time[ti],y_text,text_delta,fontsize = fs,wrap=True)
    if(exponent == 'beta') :
        text_beta = (r"$\beta/\nu_{\parallel}=$"+str(round(-slopeI,2)))
        ax.text(time[ti],y_text,text_beta,fontsize = fs,wrap=True)
    if(exponent == 'slope') :
        text_slope = ("slope = "+str(round(slopeI,2)))
        ax.text(time[ti],y_text,text_slope,fontsize = fs,wrap=True)
