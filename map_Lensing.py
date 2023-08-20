import numpy as np
import matplotlib.pyplot as plt
import pylab as py
from matplotlib import rcParams
from numpy import ma
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
import math
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.ticker import StrMethodFormatter


#nn=int(455);
dec1=float(-72.50); 
dec2=float(-64.88); 
rig1=float(71.31);
rig2=float(87.27);
#ddec=abs(dec2-dec1)/nn/1.0;
#drig=abs(rig2-rig1)/nn/1.0;

###=============================================================================
nam=[r"$Right(\rm{degree})$", r"$Declination(\rm{degree})$", r"$\epsilon_{\rm{LSST}}$[\%]", r"$t_{\rm{E}}(\rm{days})$", r"$R_{\rm{E}}(\rm{AU})$" , r"$D_{\rm{s}}(\rm{kpc})$", r"$D_{\rm{l}}(\rm{kps})$", r"$v_{\rm{rel}}(\rm{km/s})$" , r"$v_{\rm{l}}(\rm{km/s})$" ,r"$v_{\rm{s}}(\rm{km/s})$", r"$m_{\rm{l}}(M_{\odot})$", r"$u_{0}$", r"$\log_{10}[N_{\star}(\rm{deg}^{-2})]$",  r"$\log_{10}[\rho_{\star}(M_{\odot}/\rm{deg}^{2})]$", r"$\log_{10}[\epsilon_{tE}/t_{\rm{E}}]$", r"$\log_{10}[\Gamma_{\rm{obs}}(\rm{star}^{-1}.\rm{year}^{-1})]$", r"$\log_{10}[\tau\times 10^{6}]$", r"$\log_{10}[N_{\rm{event}}]$", r"$\log_{10}[N_{\star, \rm{LSST}} (\rm{deg}^{-2})]$", r"$\rm{color}(\rm{mag})$", r"$f_{\rm{giant}}[\%]$", r"$\log_{10}[\tau_{\rm{obs}}\times 10^{6}]$",  r"$m_{I}(mag)$", r"$Av$", r"$Avv$", r"$longtitude$", r"$Latitude$"]
print (len(nam))

N=int(31)

nrig, ndec= 24, 31
num=int(nrig * ndec)

dat=np.zeros((num,N))
dat=np.loadtxt('res270102A.txt') 
i=int(0)
mapp=np.zeros((27, ndec, nrig))
for k in range(nrig):
    for j in range(ndec):    
        right, decli,nstar,nbb, nevent=  dat[i,0], dat[i,1], dat[i,2], dat[i,3],  dat[i,4]
        te,re, Ds, Dl, vt, vl=           dat[i,5], dat[i,6], dat[i,7], dat[i,8],  dat[i,9], dat[i,10] 
        vs,ml, u0, Nstart, Rostart=      dat[i,11],dat[i,12],dat[i,13],dat[i,14], dat[i,15]
        tei, ratio,opt1,nfin,Nmicro=     dat[i,16],dat[i,17],dat[i,18],dat[i,19], dat[i,20]
        nstart,col,fred, opt2, Imag,nstarti=dat[i,21],dat[i,22],dat[i,23], dat[i,24], dat[i,25], dat[i,26]
        Av, Avv, longi, lati=            dat[i,27], dat[i,28], dat[i,29], dat[i,30]
        if(Rostart<0.0 or Nstart<0.0 or te<0.0 or re<0.0 or Dl>Ds or vt<0.0 or u0<0.0 or nevent>nbb or ml<0.0 or ratio<0.0 or fred>100.0 
        or fred<0.0 or Av<0.0  or Avv<0.0 ): 
            print("Error:  te, re, Ds, Dl , vt, vl, vs:  ",  te, re, Ds, Dl,  vt, vl,  vs)
            print("Error:  Rostart, NStart,   nevent, nbb, ratio:  ",  Rostart, Nstart,  nevent,  nbb, ratio )
            print("u0,  ml, fred, Imag, Av, Avv :  ",  u0,   ml, fred, Imag, Av, Avv)
            input("Enter a number ")
            
        effi=float(nevent*100.0/nbb)
        Rostart=np.log10(Rostart)
        Nstart= np.log10( Nstart)
        ratio= np.log10(ratio*opt2/opt1)
        opt1 = np.log10(opt1)
        opt2 = np.log10(opt2)
        nfin = np.log10(nfin)
        nstart=np.log10(nstart/(0.25*0.25))
        nstarti=np.log10(nstarti)
        tei=   np.log10(tei)
        mapp[:,j,k]=right,decli,effi,te,re,Ds,Dl,vt,vl,vs, ml, u0,Nstart,Rostart, tei, ratio,opt1, nfin, nstart, col, fred, opt2, Imag, Av, Avv, longi, lati 
        i+=1
print ("number of rows: ", i)
################################################################################  
print (mapp[0,0,:],  mapp[0,2,:])      
tt=int(8)
plt.clf()
fig= plt.figure(figsize=(8,7))
for k in range(27):
    v=np.zeros((tt))
    plt.cla()
    plt.clf()
    plt.imshow(mapp[k,:,:],cmap='viridis',extent=(rig1, rig2, dec2, dec1),interpolation='nearest',aspect='auto')
    plt.clim()
    plt.title(nam[k],fontsize=19)
    
    minn=np.min(mapp[k,:,:])
    maxx=np.max(mapp[k,:,:])
    step=float((maxx-minn)/(tt-1.0));
    for m in range(tt):
        v[m]=round(float(minn+m*step),1)
    cbar= plt.colorbar(orientation='horizontal',shrink=1.0,pad=0.1,ticks=v)
    cbar.ax.tick_params(labelsize=18)
    plt.clim(v[0]-0.05*step,v[tt-1]+0.05*step)
    plt.xlim( rig1, rig2)
    plt.ylim( dec1, dec2)
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    plt.xlabel(r"$Right Ascention(\rm{deg})$",fontsize=18,labelpad=0.1)
    plt.ylabel(r"$Declination(\rm{deg})$",fontsize=18,labelpad=0.1)
    plt.subplots_adjust(hspace=.0)
    fig.tight_layout(pad=0.2)
    mpl.rcParams['axes.formatter.use_mathtext'] = 'True'
    mpl.rcParams['axes.formatter.useoffset'] = 'False'
    mpl.rcParams['figure.subplot.left'] = 1.0
    mpl.rcParams['figure.subplot.right'] = .95
    mpl.rcParams['figure.subplot.bottom'] = .20
    mpl.rcParams['figure.subplot.top'] = .90
    plt.gca().invert_xaxis()
    fig=plt.gcf()
    fig.savefig("Map{0:d}.jpg".format(k),dpi=200)
    plt.clf()
    plt.cla()
    py.clf()
    print ("map is plotted  i ", k )
    print ("***************************************************" )


################################################################################
