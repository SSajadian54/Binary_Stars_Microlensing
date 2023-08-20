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
###=============================================================================
num1=[r"$b(\rm{degree})$",r"$l(\rm{degree})$",r"$\rm{Filter}$",r"$\epsilon_{Lensing}[\%]$", r"$m_{\star,detected}(\rm{mag})$",r"$m_{bg,~detected}(\rm{mag})$",r"$m_{bg,~lensed}(\rm{mag})$",r"$A(\rm{mag})$",r"$f_{\rm{b}}$",r"$\delta_{m}(\rm{mag})$",r"$N_{\rm{PSF}}$"] 
num2=[r"$~,u-\rm{band}$",r"$~,g-\rm{band}$",r"$~,r-\rm{band}$",r"$~,i-\rm{band}$",r"$~,z-\rm{band}$",r"$~,y-\rm{band}$"]


N=int(14)
num=int(25*5*33)
dat=np.zeros((num,N))
dat=np.loadtxt('res120102Bb.txt') 
nrows, ncols=25,33

i=int(0)
mapp=np.zeros((11 , 6 , nrows , ncols))
for l in range(nrows):
    for j in range(ncols):
        for k in range(6): 
            lat, lon, fil, nsim, nbb= dat[i,0], dat[i,1], dat[i,2], dat[i,3], dat[i,4]
            nsdet, neve, Map,mbagD, mbagL=dat[i,5], dat[i,6], dat[i,7], dat[i,8], dat[i,9]
            extin ,blend , dmag, nblend = dat[i,10],dat[i,11],dat[i,12], dat[i,13]
            if(lon>180.0):  lon=float(lon-360.0)
            effiL=float(neve*100.0/nsdet)
            
            if(abs(lat)>3.0  or lon<-9.0  or lon>-1.0 or fil>5.0 or nsim<800.0 or effiL>100.0 or Map<0.0  or extin<0.0 or 
            blend>1.0 or blend<0.0 or nblend<1.0 or abs(lat-(-3.0+l*0.25))>0.2 or abs(lon-(-9.0+j*0.25))>0.2): 
                print ("lat, lon, fil, nsdet, neve: " ,  lat, lon,  fil,  nsdet, neve)
                print ("Map,  mbagD, mbagL,  extin:  ",  Map, mbagD, mbagL, extin)
                print ("blend, dmag,  nblend:  ",  blend, dmag,  nblend)
                input("Enter a number ")
            mapp[:,k,l,j]= lat, lon,  fil, effiL, Map, mbagD, mbagL,  extin,  blend, dmag, nblend  
            i+=1   
           
               
print ("number of total rows: ", int(i) )
################################################################################        
plt.clf()
fig= plt.figure(figsize=(8,6))


tt=int(8)
for i in range(11):
    print ("**********************************************")
    print ("map is plotted  i ", i, "the : ", str(num1[i]) )
    for j in range(6):
        v=np.zeros((tt))
        plt.cla()
        plt.clf()
        plt.imshow(mapp[i,j,:,:],cmap='viridis',extent=(-9.0,-1.0,3.0,-3.0),interpolation='nearest',aspect='auto')
        plt.clim()
        plt.title(str(num1[i])+ str(num2[j]),fontsize=19)
        minn=np.min(mapp[i,j,:,:])
        maxx=np.max(mapp[i,j,:,:])
        step=float((maxx-minn)/(tt-1.0));
        for m in range(tt):
            v[m]=round(float(minn+m*step),1)
        cbar= plt.colorbar(orientation='horizontal',shrink=1.0,pad=0.1,ticks=v)        
        cbar.ax.tick_params(labelsize=18)
        plt.clim(v[0]-0.05*step,v[tt-1]+0.05*step)
        plt.xlim(-9.0,-1.0)
        plt.ylim(-3.0,3.0)
        plt.xticks(fontsize=18, rotation=0)
        plt.yticks(fontsize=18, rotation=0)
        plt.xlabel(r"$l(\rm{deg})$",fontsize=18,labelpad=0.1)
        plt.ylabel(r"$b(\rm{deg})$",fontsize=18,labelpad=0.1)
        plt.subplots_adjust(hspace=.0)
        fig.tight_layout(pad=0.2)
        mpl.rcParams['axes.formatter.use_mathtext'] = 'True'
        mpl.rcParams['axes.formatter.useoffset'] = 'False'
        mpl.rcParams['figure.subplot.left'] = 1.0
        mpl.rcParams['figure.subplot.right'] = .95
        mpl.rcParams['figure.subplot.bottom'] = .20
        mpl.rcParams['figure.subplot.top'] = .90
        fig=plt.gcf()
        fig.savefig("SmapB{0:d}_{1:d}.jpg".format(i,j),dpi=200)
        plt.clf()
        plt.cla()
        py.clf()
################################################################################



