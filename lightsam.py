from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
import numpy as np
import matplotlib as mpl
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import pylab
import Image
rcParams["font.size"] = 18
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
#from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('light.pdf')
import os.path

###=============================================================================


for i in range(1):
    ni=int(np.random.rand(1)*9999)
    ni=3500
    if(os.path.exists('./L-10.00_328.00_{0:d}.dat'.format(ni))==True): 
        print "file exit, number: ", ni, "step: ", i
        with open('./L-10.00_328.00_{0:d}.dat'.format(ni)) as f:
            num=int(len(f.readlines()))
        data=np.zeros((num,5))
        data=np.loadtxt('./L-10.00_328.00_{0:d}.dat'.format(ni)) 
        u0=float(data[num-7,0]); t0=float(data[num-7,1])
        tE=float(data[num-7,2]); 
        Fs= np.power(10.0,-0.4*float(data[num-7,3])); 
        mbg=float(data[num-7,4])
        Fbg= np.power(10.0,-0.4*mbg)
        blend=float(Fs/Fbg); 
        tmin=float(t0-3.5*tE)
        tmax=float(t0+3.5*tE)
        Ndat=int(10000)
        dt=float(7.0*tE/Ndat)
        model=np.zeros((Ndat,2))
        dete=int(0)
        stat=str()
        if(abs(data[num-7+3,1])>=100.0 and int(data[num-7+3,2])==1 and int(data[num-7+3,3])==1):
            stat=str('DETECTED')
        else:
            stat=str('NOTDETECTED')
        for j in range(Ndat):
            t=float(tmin+dt*j)
            u=np.sqrt(u0*u0+(t-t0)*(t-t0)/(tE*tE))
            magni=float(u*u+2.0)/(u*np.sqrt(u*u+4.0))
            model[j,0]=float(t);
            model[j,1]=-2.5*np.log10(Fs*(magni-1.0)+Fbg)
        mmin=np.min(model[:,1]); mmax=np.max(model[:,1])
        emin=np.min(data[:(num-7),2]-data[:(num-7),3]); emax=np.max(data[:(num-7),2]+data[:(num-7),3]);
        if(mmin<emin):
            emin=mmin
        if(mmax>emax):
            emax=mmax
         ###print "min of magnitude: "mmin, mmax  
       

        plt.clf()
        plt.errorbar(data[:(num-7),0]*364.5,data[:(num-7),2],yerr=data[:(num-7),3],fmt='mo',capsize=0)
        plt.plot(model[:,0],model[:,1],color="g",label="light curve", lw=2, alpha=0.95)
        plt.ylabel('r_{LSST}-magnitude',fontsize=18)
        plt.xlabel('time (day)',fontsize=18)
        nam=[r"$u_{0}$",r"$,~t_{0}$",r"$,~t_{E}$",r"$,~m_{bg,r}$",r"$,~b$"]
        plt.title(str(nam[0])+'={0:.2f}'.format(u0)+ str(nam[1]) + '={0:.1f}'.format(t0)+str(nam[2])+'={0:.1f}'.format(tE)+str(nam[3])+ '={0:.1f}'.format(mbg)+str(nam[4]) + '={0:.2f}'.format(blend),fontsize=16,color='b')
        plt.text(tmax-2.0*tE,emin+(emax-emin)*0.15,str(stat), color='g',fontsize=16)
        pylab.xlim([tmin,tmax])
        pylab.ylim([emin-0.05*(emax-emin),emax+0.05*(emax-emin)])
        plt.gca().invert_yaxis()
        fig4=plt.gcf()
        fig4.savefig("./light_example/dLight{0:d}.png".format(ni),dpi=200)
        Image.open('./light_example/dLight{0:d}.png'.format(ni)).save('./light_example/dLight{0:d}.jpg'.format(ni),'JPEG')




