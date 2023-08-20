from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
import numpy as np
import emcee
import matplotlib as mpl
import pylab
rcParams["font.size"] = 18
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
import os.path
###=============================================================================
nam=[r"$u_{0, 1}$",r"$,~u_{0, 2}$",r"$,~\Delta t_{0}$",r"$,~\Delta mag$", r"$,~\Delta \chi^{2}$"]

nf=444241
par=np.zeros(( nf , 29 ))
par=np.loadtxt("./files/distribution/P_0.dat")   

Nums=np.zeros((14))
Numd=np.zeros((14))

Nsim=0
case=0.0
Ntot=0.0
for i in range(nf):
    Nsim,strucl,Ml,RE, Dl, vl, tetE, ksi  =      int(par[i,0]), par[i,1], par[i,2], par[i,3], par[i,4], par[i,5], par[i,6], par[i,7]  
    strucs, cl, Ds, vs, Av, AI, Map1, Map2=par[i,8], par[i,9],par[i,10], par[i,11],par[i,12],par[i,13], par[i,14], par[i,15]
    semi, Vt, tE, u0, u02, t0 , t02, opt  =par[i,16],par[i,17],par[i,18],par[i,19],par[i,20],par[i,21], par[i,22], par[i,23]
    dchi, det, peak, Dmag, jd =  par[i,24], par[i,25], par[i,26], par[i,27], par[i, 28]
    Ntot+=1.0
    print("****************************************************")
    print("Counter:  ",  Nsim)
    print("t0,  t02,  u0,   u02:   ",  t0,  t02,  u0, u02 )
    print("semi, Dl, Ds:   ", semi,    Dl,   Ds)
    print("dchi,   det,   peak, dmag:    ",dchi, det, peak , Dmag, jd)
    print("****************************************************")
    Dt0=t02-t0
    #Dmag=abs(Map1-Map2)
    '''
    nn=-1
    if(Dmag>deltam[16] or Dmag==deltam[16]):  nn=16;  
    else: 
        for k in range(16): 
            if(float((Dmag-deltam[k])*(Dmag-deltam[k+1]))<0.0 or Dmag==deltam[k]): 
                nn=k
                break 
    if(nn<0): 
        print("Big error nn:  ", nn,   Dmag  )
    '''
    Nums[int(jd)-1]+=1; 
    if(dchi>=300.0 and det>0):  Numd[int(jd)-1]+=1    
    
    
    
    nd=0;  nl=0;  
    if(Nsim%1000==0):  
    #if ( open('./files/light/d_{0:d}.dat'.format(i))==True): 
        f1=open('./files/light/d_{0:d}.dat'.format(Nsim))
        nd=int(len(f1.readlines()))
    #if ( open('./files/light/l_{0:d}.dat'.format(i))==True ): 
        f2= open('./files/light/l_{0:d}.dat'.format(Nsim))
        nl=int(len(f2.readlines()))  
         
    if(nd>1 and nl>1): 
        mod=np.zeros((nl,5))
        mod=np.loadtxt('./files/light/l_{0:d}.dat'.format(Nsim)) 
    
        dat=np.zeros((nd,6))
        dat=np.loadtxt('./files/light/d_{0:d}.dat'.format(Nsim)) 
        
        stat=[]
        if(dchi>=300.0 and det>0):
            stat=str('DETECTED')
            col='g'
            case+=1        
        else:
            stat=str('NOT-DETECTED')
            col='r'
        '''   
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
        '''        
        
        emax, emin= np.max(mod[:,2]), np.min(mod[:,2])
        plt.clf()
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        
        plt.plot(mod[:,1]*tE,mod[:,2],'g-',label="Binary stars", lw=1.5, alpha=0.95)
        plt.plot(mod[:,1]*tE,mod[:,3],'r--',label="Star 1", lw=1.5, alpha=0.95)
        plt.plot(mod[:,1]*tE,mod[:,4],'b:',label="Star 2", lw=1.5, alpha=0.95)
        
        plt.errorbar(dat[:,1]*tE,dat[:,2],yerr=dat[:,3], fmt=".", markersize=10.8,  color='m', ecolor='gray', elinewidth=0.3, capsize=0)
       
        plt.ylabel(r"$r_{\rm{LSST}}-\rm{magnitude}$",fontsize=18)
        plt.xlabel(r"$\rm{time}(\rm{days})-t_{0, 1}$",fontsize=18)

        plt.title(str(nam[0])+'={0:.2f}'.format(u0)+ str(nam[1])+'={0:.1f}'.format(u02)+str(nam[2])+ '={0:.1f}'.format(Dt0) + str(nam[3])+ '={0:.1f}'.format(Dmag) +str(nam[4])+ '={0:.1f}'.format(dchi) , fontsize=16,color='b')

        
        plt.text(np.min(mod[:,1]*tE)+0.5,emax-(emax-emin)/2.0,str(stat), color=col , fontsize=13)

        
        pylab.ylim([emin-0.02*(emax-emin),emax+0.02*(emax-emin)])
        pylab.xlim([np.min(mod[:,1]*tE),  np.max(mod[:,1]*tE)  ])
        
        plt.gca().invert_yaxis()
        #ax1.grid("True")
        #ax1.grid(linestyle='dashed')
        ax1.legend(prop={"size":12.5})
        fig=plt.gcf()
        fig.savefig("./files/light/light_example/Lbinarys{0:d}.jpg".format(Nsim),dpi=200)
        print("*************************************")
        #input("Enter a number ")
  
print(case*100.0/Ntot )

for i in range(14):  
    Numd[i]=float(Numd[i]*100.0/(Nums[i]+0.00000034756234564) )
    print ("effi:  ", i,  Numd[i],   Nums[i])
    
######################################################################    
dd=np.array([0.3,0.6,0.9,1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.2])
plt.clf()
fig=plt.figure(figsize=(8,6))    
plt.scatter(dd , Numd, s=34.0,marker="*", color="g")
plt.ylabel(r"$\rm{Efficiency}[\%]$",fontsize=18)
plt.xlabel(r"$\Delta m (\rm{mag})$",fontsize=18)
#plt.xlim([0.0,10.0])
#plt.yscale('log')
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.grid(True)
plt.grid(linestyle='dashed')  
fig=plt.gcf()
#fig.tight_layout()
fig.savefig("./effi.jpg", dpi=200)
        

        


    
    

