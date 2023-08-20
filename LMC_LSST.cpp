#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
using namespace std;

const int Num=10000;
const double MaxD=65.0;///kpc
const double step=MaxD/(double)Num/1.0;///step in kpc

const double RA=180.0/M_PI;
const double pi= M_PI; 
const double binary_fraction=double(2.0/3.0);
const double velocity=299792458.0;//velosity of light
const double Msun=1.98892*pow(10.,30); //in [kg].
const double Rsun=6.957*pow(10.0,8.0); ///solar radius [meter]
const double KP=3.08568025*pow(10.,19); // in meter.
const double G= 6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double AU=1.4960*pow(10.0,11.0);
const double logT_sun=log10(5778.0);
const double vro_sun=226.0;
const double VSunR =11.1;
const double VSunT =vro_sun*(1.00762+0.00712)+ 12.24;
const double VSunZ =7.25;
const double year=365.2425;//days 
const int nrd=10000; 
const int NM=20000;///number of rows in file 'Monte.txt'
const int M=6;///No. of filter  ugrizy
const int LL=5;
const int GG=900;///number of bins of tE distribution
const double tE_min=0.0;
const double tE_max=500.0;  
const double stept=double((tE_max-tE_min)/GG); 

///============================ Besancon constant ==============///

const double Dsun=8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[7]={3.1,2.5,3.1,3.1,3.1,3.1,2.5};
const int    N1=27004, N2=36558, N3=2612, N4=492;///CMD_BESANCON, thinD, bulge, thickD, halo

///===============LSST constant===============================///

const double Tobs=10.0*year;///LSST observational time 10 years
const double cadence=3.0;//in days 
const double texp=30.0;///in second
const double Delta2=0.005;///systematic errors
const double gama[M]={0.037,0.038,0.039,0.039,0.040,0.040};///a function of sky brightness and airmass in different wavelengths.
const double seeing[M]={0.77,0.73,0.70,0.67,0.65,0.63};///seeing
const double msky[M]={22.9,22.3,21.2,20.5,19.6,18.6};
const double Cm[M]=  {22.92,24.29,24.33,24.20,24.07,23.69};
const double Dci[M]= {0.67,0.21,0.11,0.08,0.05,0.04};
const double km[M]=  {0.451,0.163,0.087,0.065,0.043,0.138};///sky extinction
const double wav[M+1]={0.3543,0.477,0.6231,0.7625,0.9134,1.004,0.7865};///in[micrometer] u g r i z y Ic
const double sigma[M+1]={0.022,0.02,0.017,0.017,0.027,0.027,0.017};
const double detect[M]={23.4,24.6,24.3,23.6,22.9,21.7};///depth of single visit in ugriz
const double satu[M]=  {15.2,16.3,16.0,15.3,14.6,13.4};///saturation limit of single visit in ugriz
const double M50[M]=   {23.68,24.89,24.43,24.00,24.45,22.60};
const double FWHM[M+1]={1.22087,1.10136,0.993103,0.967076,0.951766,0.936578,0.967076*1.65};//LSST [arcsec]
const int YZ=3578;//// No.yzma.txt rows

////========================== LMC =============================

const double bLMC= -32.9; 
const double lLMC= 280.5; 
const double RaLMC= 80.89375;
const double DecLMC=-68.2438888888889;  
const double sLMC= 5.0;///degree (half of size of LMC)
const double DLMC= 49.97;///KPC 

const int nn=int(455);
const double dec1= -72.5; 
const double dec2= -64.8; 
const double rig1=float((4.0+27.0/60.0)*15.0 );
const double rig2=float((5.0+50.0/60.0)*15.0 );
const double ddec=5.0*fabs(dec2-dec1)/nn/1.0;
const double drig=5.0*fabs(rig2-rig1)/nn/1.0;


const int nex=23104; 

///======================================================
struct source{
    int nums,struc, cl;
    double Ds,TET,FI, lat, lon;
    double right, decli;
    double od_disk,od_ThD,od_bulge,od_halo,opt;///optical  depth
    double od_dlmc, od_blmc, od_hlmc; 
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double rho_hlmc[Num], rho_dlmc[Num],rho_blmc[Num]; 
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs, Romins,nstart, nstarti;
    double nsdis[Num],nddis[Num];
    double Nblend[M],blend[M], Fluxb[M], magb[M], Ai[M], Mab[M], Map[M], nsbl[M]; 
    double logl, logt, col, vs;
    double u0, ro_star, Rs;
    double FIsource, Isource, blendI, MIsource, magni; 
    double xv, yv, zv; 
    double Avv, Av, AI, A[M]; 
};

struct lens{
    double Ml, Dl, vl, Vt, xls;
    double rhomaxl, tE, RE;
    int numl, struc;
};

struct detection{
   int det,peak;
   double deltam[M],Dcm[M],dchi;
   double t,tmin,tmax,t0;
   double ave_re,ave_vl,ave_vs;
   double ave_opt,ave_opt2,ave_te,ave_u0,ave_Ds,ave_Dl;
   double ave_Av, ave_Avv; 
   double ave_vt,ave_ml,ave_col,aveI; 
   double ave_aps[M],ave_ex[M],ave_dmag[M],neve[M],ave_apbd[M],ave_bl[M],ave_apb[M];
   double ratio,nsdet[M],nstar,nevent,nfin; 
   double ave_npsf[M]; 
   double m5c[NM]; int filter[NM];
   double ave_tei;///<effi/tE>
   double te[GG+1],effs[GG+1],effd[GG+1]; 
   double effs_com[GG+1],effd_com[GG+1];
   double ave_fred;  
   double weight[M],ave_w;
   int N_sim,N_det;
};
struct CMD{
    double logt_d[N1],logl_d[N1],Mab_d[M][N1],col_d[N1],MI_d[N1]; int cl_d[N1];  ///thin disk
    double logt_b[N2],logl_b[N2],Mab_b[M][N2],col_b[N2],MI_b[N2]; int cl_b[N2];  /// bulge
    double logt_t[N3],logl_t[N3],Mab_t[M][N3],col_t[N3],MI_t[N3]; int cl_t[N3];  ///thick disk
    double logt_h[N4],logl_h[N4],Mab_h[M][N4],col_h[N4],MI_h[N4]; int cl_h[N4];  /// halo
};
 
struct galactic{
   double l[nrd] , b[nrd], RA[nrd], Dec[nrd]; 
};

struct extin{

  double right[nex], dec[nex], ext[nex];
};

///===================== FUNCTION ==============================================
void read_cmd(CMD & cm);
void func_source(source & s, CMD & cm);
void func_lens(lens & l,   source & s);
void vrel(source & s,   lens & l);
void Disk_model(source & s, int);
void optical_depth(source & s);
void lensing(source & s, lens & l,detection & d,int);
double RandN(double , double);
double RandR(double , double);
/////////////////////////////////////////////
    time_t _timeNow;
      unsigned int _randVal;
    unsigned int _dummyVal;
    FILE * _randStream;
////////////////////////////////////////////
///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//
int main()
{
///****************************************************************************
//gettimeofday(&newTime,NULL);
//ftime(&tp);
   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
///****************************************************************************

   CMD cm;  
   source s;
   lens l;
   detection d;
   galactic ga;
   extin ex; 
   read_cmd(cm);        
   
   char filename1[40], filename2[40]; 
   FILE* result1; FILE* result2;  FILE* result3; 
   FILE* teefi;  FILE* TEE;   
   FILE* fil1;   FILE* fil2; 
   result1=fopen("./files/res270102A.txt", "w");
   result2=fopen("./files/res270102B.txt", "w");
   result3=fopen("./files/res270102C.txt", "w");
   fclose(result1); 
   fclose(result2);
   fclose(result3);
     
     
   for(int i=0; i<M; ++i){
   if(texp!=30.0) d.Dcm[i]=Dci[i]-1.25*log10(1.0+(pow(10.0,0.8*Dci[i])-1.0)/(texp/30.0));
   else d.Dcm[i]=0.0;}
   
   FILE* convert; 
   convert=fopen("./files/convert_coordinate_2.dat","r");
   int coo=0; 
   for(int i=0; i<100; ++i){
   for(int j=0; j<100; ++j){
   fscanf(convert,"%lf  %lf   %lf  %lf \n",&ga.RA[coo], &ga.Dec[coo], &ga.l[coo], &ga.b[coo]);  coo++; }}   
   fclose(convert); 
   
   

   FILE* Extlmc; 
   Extlmc=fopen("./files/extinctionf.txt","r");
   for(int i=0; i<nex;   ++i){
   fscanf(Extlmc,"%lf    %lf     %lf \n",&ex.right[i], &ex.dec[i], &ex.ext[i] );}
   fclose(Extlmc); 
   
   
   
   FILE* Mont;
   Mont=fopen("./files/Monte2.txt","r");
   for(int i=0; i<NM; ++i){
   fscanf(Mont,"%lf   %d\n",&d.m5c[i],&d.filter[i]);  
   if(d.m5c[i]<20.0 or d.m5c[i]>27.0 or d.filter[i]>=6 or d.filter[i]<0){ 
   cout<<"ERROR airmass: "<<d.m5c[i]<<"\t filter: "<<d.filter[i]<<endl;  int yye; cin>>yye; }} 
   fclose(Mont);
   
   
   
   double nbb,effe,maga[M],lte,sumd=0.0; 
   int vf,ffd[M],flag_det, l1a, l2a, hh;
   double Nmicro=0.0;
   double mindd, lonn,countt,countd, we, dist; 
   double exti[455]={0.0}; 
   
    int uue=0; 
   
//===================== Monte Carlo Simulation ================================
   for(s.right= rig1;   s.right<rig2;  s.right+=drig*3.0){
   uue+=1; 
   for(s.decli= dec1;   s.decli<dec2;  s.decli+=ddec*3.0){
   
  // s.right= RaLMC;  
  // s.decli= DecLMC; 
  
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   
   
   //l1a=int( double(rig2-s.right)/drig +1.0); 
   //l2a=int( double(dec2-s.decli)/ddec +1.0); 
   if(s.right>71.0 and s.right<90 and s.decli<-60.0 and s.decli>-76.0){
   
   cout<<">>>>>>>>>>>>>>>> NEW STEP <<<<<<<<<<<<<<<<"<<endl;
   cout<<"Right ascention: "<<s.right<<"\t declination: "<<s.decli<<endl;
   cout<<"l1a:  "<<l1a<<"\t l2a:   "<<l2a<<endl;
   cout<<"rig1:  "<<rig1<<"\t rig2:  "<<rig2<<"\t drig: "<<drig<<endl; 
   cout<<"dec1:  "<<dec1<<"\t dec2:  "<<dec2<<"\t ddec:  "<<ddec<<endl; 
 // int uue; cin>>uue;    
    
    
   hh=-1; mindd=10000.0; 
   for(int i=0; i<10000; ++i){
   dist= double(s.right- ga.RA[i])*(s.right- ga.RA[i]) + (s.decli- ga.Dec[i])*(s.decli- ga.Dec[i]); 
   dist=sqrt(dist); 
   if(dist<=mindd){ mindd=dist; hh= i;   }}
   if(hh<0){cout<<"ERROR:right:    "<<s.right<<"\t s.dec:  "<<s.decli<<endl;   int uue;  cin>>uue; }
   s.lon=double(ga.l[hh]);   
   s.lat=double(ga.b[hh]);
   s.TET=(360.0-s.lon)/RA;///radian
   //s.TET= s.lon/RA; 
   s.FI=  s.lat/RA; 
   Disk_model(s, uue); 
   


   double ddd=100000.0, dis; 
   for(int i=0;  i<nex; ++i){
   dis=double(ex.right[i]-s.right)*(ex.right[i]-s.right) + (ex.dec[i]-s.decli)*(ex.dec[i]-s.decli) ; 
   dis=sqrt(dis);  
   if(dis<ddd){ddd=dis;   s.Avv= ex.ext[i]; }}
   if(s.Avv<0.0 )  s.Avv=0.00004352475613;  
   
   
   d.nstar=0.00001745693;  nbb=0.0; 
   d.ave_Av= d.ave_Avv=0.0; 
   d.ave_Ds=d.ave_Dl=d.ave_vt=d.ave_ml=d.ave_u0=d.ave_tei=0.0;
   d.ave_re=d.ave_vl=d.ave_vs=d.ave_te=d.nevent=d.ave_col=0.0;
   d.ave_opt=d.ave_opt2=d.aveI=s.nstart=s.nstarti=d.ratio=d.nfin=0.0; d.ave_fred=0.0;
   for(int j=0; j<M; ++j){ 
   d.ave_apb[j]=d.ave_apbd[j]=d.ave_dmag[j]=d.ave_npsf[j]=0.0; 
   d.ave_bl[j]=d.ave_aps[j]=d.ave_ex[j]=0.0; d.nsdet[j]=d.neve[j]=0.00000645974;}
   for(int i=0;  i<Num; ++i){s.nsdis[i]=0.0000145776534762;   s.nddis[i]=0.0;}
   for(int i=0; i<(GG+1);++i){d.effd[i]=0.0;  d.effs[i]=0.0000000754765;  d.te[i]=double(i*stept);}
   d.N_sim=d.N_det=0; 


   
   s.Avv= exti[l2a-1] + RandR(-0.1,0.1);
   if(s.Avv<0.0) s.Avv=0.05; 
   if(s.Avv>2.0) s.Avv=2.0 ;   
   cout<<"Avv:   "<<s.Avv<<"\t l1a:  "<<l1a<<"\t l2a:  "<<l2a<<endl;
   
  

   sprintf(filename1,"./files/distribution/%c%.2lf%c%.2lf.dat",'W',s.right,'_',s.decli);
   fil1=fopen(filename1,"w");
   sprintf(filename2,"./files/distribution/%c%.2lf%c%.2lf.dat",'V',s.right,'_',s.decli);
   fil2=fopen(filename2,"w");
   if(!fil1 or !fil2){cout<<"cannot make file longtitude : "<<s.right<<"\t latitude: "<<s.decli<<endl; exit(0);}
 
   teefi=fopen("./files/tefzzz.txt","w"); 
   fclose(teefi);
   
  
   
   do{
   d.nstar+=1.0;

   do{
   func_source(s,cm);
   func_lens(l,s);
   }while(l.tE<=0.1 or l.tE>tE_max); 
   optical_depth(s);
  
  
  
   s.nsdis[s.nums]+=1.0;
   d.ave_w=0.0;  
   sumd=0.0;    
   for(int i=0; i<M; ++i){
   ffd[i]=0; 
   d.weight[i]=0.0; 
   maga[i]=-2.5*log10(fabs(pow(10.0,-0.4*s.Map[i])*(s.magni-1.0)+s.Fluxb[i]));
   if(s.magb[i]>=satu[i] and maga[i]<=detect[i]){
   d.weight[i]=fabs(1.0/s.nsbl[i]); 
   ffd[i]=1;  
   sumd+= ffd[i]; 
   d.ave_w+=fabs(d.weight[i]);  
   d.nsdet[i]   +=fabs(1.0000*   d.weight[i]); 
   d.ave_aps[i] +=fabs(s.Map[i]* d.weight[i]);
   d.ave_apb[i] +=fabs(s.magb[i]*d.weight[i]);
   d.ave_ex[i]  +=fabs(s.Ai[i] * d.weight[i]);
   d.ave_npsf[i]+=fabs(s.nsbl[i]*d.weight[i]);}}



   if(ffd[0]>0 or ffd[1]>0 or ffd[2]>0 or ffd[3]>0 or ffd[4]>0){    
   
   d.ave_w=   fabs(d.ave_w/sumd);  
   s.nddis[s.nums]+=1.0*d.ave_w; 
   nbb +=     fabs(1.00*d.ave_w);
   d.ave_opt+=fabs(s.opt*1.0e6*d.ave_w); 
   d.ave_col+=fabs(s.col*d.ave_w);
   d.aveI +=  fabs(s.Isource*d.ave_w);
   d.N_sim+=1;
   lensing(s,l,d,d.N_sim);


   flag_det=0; 
   if(d.det==1 and d.peak==1 and d.dchi>=200.0){
   flag_det= 1;     
   d.N_det +=1; 
   for(int i=0 ; i<M; ++i){
   d.neve[i] +=   fabs(1.0000000*  d.weight[i]); 
   d.ave_bl[i]+=  fabs(s.blend[i]* d.weight[i]); 
   d.ave_dmag[i]+=fabs(d.deltam[i]*d.weight[i]);
   d.ave_apbd[i]+=fabs(s.magb[i]*  d.weight[i]);}
   d.nevent+=fabs(1.00*d.ave_w); 
   d.ave_te+=fabs(l.tE*d.ave_w);
   d.ave_Ds+=fabs(s.Ds*d.ave_w);
   d.ave_Dl+=fabs(l.Dl*d.ave_w);
   d.ave_vt+=fabs(l.Vt*d.ave_w);
   d.ave_vl+=fabs(l.vl*d.ave_w);
   d.ave_vs+=fabs(s.vs*d.ave_w);
   d.ave_ml+=fabs(l.Ml*d.ave_w);
   d.ave_u0+=fabs(s.u0*d.ave_w);   
   d.ave_Av+=fabs(s.Av*d.ave_w); 
   d.ave_Avv+=fabs(s.Avv*d.ave_w); 
   d.ave_opt2+=fabs(s.opt*1.0e6*d.ave_w); 
   d.ave_re+=fabs(l.RE/AU)*d.ave_w;
   if(s.cl<=4) d.ave_fred+=1.0*d.ave_w;}
   
    vf=-1; 
    for(int j=0; j<GG; ++j){
    if((l.tE-d.te[j])*(l.tE-d.te[j+1])<0.0 or l.tE==d.te[j]){vf=j;  j=GG+3;}  }
    if(vf==-1 and (l.tE>tE_max or l.tE==tE_max)) vf=GG; 
    if(vf==-1 and (l.tE<tE_min or l.tE==tE_min)) vf=0;    
    if(vf<0 or l.tE<0.0 or vf>GG){cout<<"ERROR tE: "<<l.tE<<"\t vf: "<<vf<<endl;  int uue; cin>>uue;}
    d.effs[vf]    +=1.0;  
    d.effs_com[vf]+=1.0;
    if(flag_det>0.5){d.effd[vf] +=1.0;  d.effd_com[vf]+=1.0;}
    teefi=fopen("./files/tefzzz.txt","a+");
    fprintf(teefi,"%.7lf   %e  %d  %d\n",l.tE,d.ave_w,vf,flag_det);
    fclose(teefi);


    if(flag_det==1){
    fprintf(fil1,"%d %d %d %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.6lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %d %d   %.8lf\n",
    d.N_sim,s.struc,l.struc,l.Ml,l.RE/AU, l.Dl, s.Ds, l.Vt, l.tE, l.vl, s.vs,
    s.logl,s.logt,s.Rs,s.ro_star,s.col,s.Av,s.AI,s.u0,d.t0,s.opt*1.0e6,s.cl,flag_det,d.ave_w);
    for(int i=0; i<M; ++i)
    fprintf(fil2,"%d  %.6lf  %.8lf  %.6lf  %.4lf   %.2lf  %.4lf  %.4lf  %d   %.8lf\n",
    d.N_sim,s.magb[i],s.Ai[i],d.deltam[i],s.blend[i],s.nsbl[i],s.Map[i],s.Mab[i],flag_det,d.weight[i]);
    }
     
     /*if(int(d.N_det)%100000==0 and  flag_det==1){
     cout<<"**********************************************"<<endl;
     cout<<"detecable: "<<sumd<<"\t flag_det: "<<flag_det<<"\t d.nevent: "<<d.nevent<<endl;
     cout<<"nbb: "<<nbb<<"\t ave_w: "<<d.ave_w<<"\t N_det: "<<d.N_det<<endl;
     cout<<"weight0: "<<d.weight[0]<<"\t weight1: "<<d.weight[1]<<"\t weight2: "<<d.weight[2]<<endl;
     cout<<"weight3: "<<d.weight[3]<<"\t weight4: "<<d.weight[4]<<endl;
     cout<<"latitude: "<<s.lat<<"\t longtitude: "<<s.lon<<endl;
     cout<<"Ml (Msun): "<<l.Ml<<"\t RE (AU): "<<l.RE/AU<<"\t Dl(Kpc): "<<l.Dl<<endl;
     cout<<"relative velosity (km/s): "<<l.Vt<<"\t Vs (km/s): "<<s.vs<<"\t Vl(km/s): "<<l.vl<<endl;
     cout<<"Ds: "<<s.Ds<<"\t tE (days): "<<l.tE<<"\t ro_star"<<s.ro_star<<endl;
     cout<<"app_mag[r]: "<<s.Map[2]<<"\t app_mag(i): "<<s.Map[3]<<endl;
     cout<<"app_mag_bg[r]: "<<s.magb[2]<<"\t app_mag_bg(i): "<<s.magb[3]<<endl;
     cout<<"belnding_r(factor): "<<s.blend[2]<<"\t belnding_i: "<<s.blend[3]<<endl;
     cout<<"optical_depth(x10^6): "<<s.opd*1.0e6<<"\t u0: "<<s.u0<<endl;
     cout<<"N_blend[r_band]: "<<s.nsbl[2]<<"\t N_blend[i_band]: "<<s.nsbl[3]<<endl;
     cout<<"det: "<<d.det<<"\t peak : "<<d.peak<<"\t dci2: "<<d.dchi<<endl;
     cout<<"maga[0]:"<<maga[0]<<"\t maga[1]: "<<maga[1]<<"\t maga[2]: "<<maga[2]<<endl;
     cout<<"maga[3]: "<<maga[3]<<"\t maga[4]: "<<maga[4]<<endl;
     cout<<"*********************************************************************"<<endl; }*/
     if(sumd>M or d.ave_w>1.0001 or d.ave_w<=0.0 or s.Ds>MaxD or l.tE<0.0 or l.tE>tE_max or 
     s.blend[0]>1.0001 or s.blend[1]>1.0001 or s.blend[2]>1.0001 or s.blend[3]>1.0001 or s.blend[4]>1.0001){
     cout<<"ERROR sumd: "<<sumd<<"\t ave_w: "<<d.ave_w<<"\t Ds: "<<s.Ds<<endl;  
     cout<<"b0: "<<s.blend[0]<<"\t b1: "<<s.blend[1]<<"\t b2: "<<s.blend[2]<<"\t b3: "<<s.blend[3]<<"\t b4:"<<s.blend[4]<<endl;
     int uue;  cin>>uue; }
     }
     }while(d.N_det<100); 
     fclose(fil1); 
     fclose(fil2);



     TEE=fopen("./files/tef_last_finalb.txt","w"); 
     for(int j=0;j<(GG+1); ++j){
     d.effd[j]=double(d.effd[j]/d.effs[j]);
     fprintf(TEE,"%d   %.3lf     %.2lf    %.2lf\n",j,d.te[j],d.effs_com[j],d.effd_com[j]);}     
     fclose(TEE);
     


     d.ave_tei=0.0; countd=0.0;
     teefi=fopen("./files/tefzzz.txt","r"); 
     if(!teefi){cout<<"cannot make file teefi:  "<<endl; exit(0);}
     while(!feof(teefi)){
     fscanf(teefi,"%lf  %lf %d  %d\n",&l.tE,&we,&vf,&flag_det);
     if(flag_det>0.5 and l.tE>0.0){
     countd += we;  
     d.ave_tei+=fabs(d.effd[vf]/l.tE)*we;}}
     d.ave_tei=double(d.ave_tei/countd); 
     fclose(teefi);
 

    ///HHHHHHHHHHHHHHHHHHHHHHHHHH SAVING IN FILE HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    result1=fopen("./files/res270102A.txt", "a+");
    result2=fopen("./files/res270102B.txt", "a+"); 
    result3=fopen("./files/res270102C.txt", "a+"); 
    
    d.ave_te= d.ave_te/d.nevent;
    d.ave_re= d.ave_re/d.nevent;
    d.ave_Ds= d.ave_Ds/d.nevent;
    d.ave_Dl= d.ave_Dl/d.nevent;
    d.ave_vt= d.ave_vt/d.nevent;
    d.ave_vl= d.ave_vl/d.nevent;
    d.ave_vs= d.ave_vs/d.nevent;
    d.ave_ml= d.ave_ml/d.nevent;
    d.ave_u0= d.ave_u0/d.nevent;
    d.ave_Av= d.ave_Av/d.nevent;
    d.ave_Avv=d.ave_Avv/d.nevent;
    
    
    d.ave_opt2=d.ave_opt2/d.nevent;
    d.ave_fred=d.ave_fred*100.0/d.nevent;
    d.ave_opt=double(d.ave_opt/nbb);
    d.ave_col=double(d.ave_col/nbb); 
    d.aveI=   double(d.aveI/nbb); 
    d.ratio= 2.0*d.ave_opt*1.0e-6*d.ave_tei*year/M_PI;///[N_event/year]
    for(int i=1; i<Num; ++i){
    effe=double(s.nddis[i]/(s.nsdis[i]+0.0000064646745));
    s.nstart += s.Rostari[i]*(s.Nstart/s.Rostart)*effe*ddec*drig;}
    d.nfin=d.ratio*10.0*s.nstart;////No.event  
    Nmicro+=d.nfin;///No. event_total
    
    fprintf(result1,"%.2lf %.2lf %.1lf %.1lf %.1lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %e  %e %.8lf %e %e %e %e  %e  %.5lf  %.3lf   %.5lf   %.5lf  %e     %.7lf  %.7lf  %.6lf  %.6lf\n",
    s.right,s.decli,d.nstar,nbb,d.nevent,d.ave_te,d.ave_re,d.ave_Ds,d.ave_Dl,d.ave_vt,d.ave_vl,  
    d.ave_vs,d.ave_ml,d.ave_u0,s.Nstart,s.Rostart,d.ave_tei,d.ratio,d.ave_opt,d.nfin,Nmicro,
    s.nstart,d.ave_col,d.ave_fred,d.ave_opt2,d.aveI,s.nstarti, d.ave_Av, d.ave_Avv,  s.lon, s.lat );
    for(int j=0; j<M; ++j){
    d.ave_aps[j]= d.ave_aps[j]/d.nsdet[j];
    d.ave_apb[j]= d.ave_apb[j]/d.nsdet[j];
    d.ave_ex[j]=  d.ave_ex[j]/d.nsdet[j];
    d.ave_npsf[j]=d.ave_npsf[j]/d.nsdet[j];
    d.ave_apbd[j]=d.ave_apbd[j]/d.neve[j];
    d.ave_bl[j]=  d.ave_bl[j]/d.neve[j];
    d.ave_dmag[j]=d.ave_dmag[j]/d.neve[j];
    fprintf(result2,"%.2lf %.2lf %d %.6lf  %lf  %.1lf  %.2lf  %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.2lf\n",
    s.right,s.decli,j,d.nstar,nbb,d.nsdet[j],d.neve[j],d.ave_aps[j],d.ave_apb[j],d.ave_apbd[j],d.ave_ex[j],d.ave_bl[j],
    d.ave_dmag[j],d.ave_npsf[j]);}


///==========================================================================

    sprintf(filename1,"./files/distribution/%c%.2lf%c%.2lf.dat",'W',s.right,'_',s.decli);
    fil1=fopen(filename1,"r");
    if(!fil1){cout<<"cannot make file longtitude : "<<s.right<<"\t latitude: "<<s.decli<<endl; exit(0);}
    countt=0.0;   countd=0.0;
    double SD_Dl=0.0, SD_Ds=0.0, SD_ml=0.0,  SD_tE=0.0,  SD_vt=0.0; 
    while(!feof(fil1)){
    countt+=we;  
    if(we>0.0) countd+=1.0;
    fscanf(fil1,"%lf  %d  %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d  %lf\n",
    &nbb,&s.struc,&l.struc,&l.Ml,&l.RE,&l.Dl,&s.Ds,&l.Vt,&l.tE,&l.vl,&s.vs,
    &s.logl,&s.logt,&s.Rs,&s.ro_star,&s.col,&s.Av,&s.AI,&s.u0,&d.t0,&s.opt,&s.cl,&flag_det,&we);///24
    SD_Dl+=pow(l.Dl-d.ave_Dl,2.0)*we; 
    SD_Ds+=pow(s.Ds-d.ave_Ds,2.0)*we; 
    SD_ml+=pow(l.Ml-d.ave_ml,2.0)*we;  
    SD_tE+=pow(l.tE-d.ave_te,2.0)*we; 
    SD_vt+=pow(l.Vt-d.ave_vt,2.0)*we;}
    SD_Dl= sqrt(SD_Dl/(countt*(countd-1.0)/countd));
    SD_Ds= sqrt(SD_Ds/(countt*(countd-1.0)/countd));
    SD_ml= sqrt(SD_ml/(countt*(countd-1.0)/countd));
    SD_tE= sqrt(SD_tE/(countt*(countd-1.0)/countd));
    SD_vt= sqrt(SD_vt/(countt*(countd-1.0)/countd));
    fclose(fil1);   
    cout<<"<Dl>: "<<d.ave_Dl<<"\t SD_Dl: "<<SD_Dl<<endl;
    cout<<"<Ds>: "<<d.ave_Ds<<"\t SD_Ds: "<<SD_Ds<<endl;
    cout<<"<ml>: "<<d.ave_ml<<"\t SD_ml: "<<SD_ml<<endl;
    cout<<"<tE>: "<<d.ave_te<<"\t SD_tE: "<<SD_tE<<endl;
    cout<<"<vt>: "<<d.ave_vt<<"\t SD_vt: "<<SD_vt<<endl;
    fprintf(result3,"%.2lf %.2lf   %.5lf  %.5lf  %.5lf  %.5lf  %.5lf\n",s.right,s.decli,SD_Dl, SD_Ds, SD_ml, SD_tE, SD_vt);//7
    fclose(result1);  fclose(result2);  fclose(result3); 
///============================================================================
    
    cout<<" \n\n******  Right ascention : "<<s.right<<"\t declination: "<<s.decli<<endl;
    for (int  j=0; j<M; ++j){
    cout<<"********** LSST FILTER ***************************: "<<j<<endl; 
    cout<<"nsdet: "<<d.nsdet[j]<<endl;
    cout<<"<M_bg_detected>: "<<d.ave_apb[j]<<"\t <M_bg_lensed>: "<<d.ave_apbd[j]<<endl;
    cout<<"<Map_sources>:"<<d.ave_aps[j]<<"\t <EXtinction>: "<<d.ave_ex[j]<<endl;
    cout<<"effi star-det: "<<d.nsdet[j]/d.nstar<<"\t effi_ML: "<<d.neve[j]/d.nsdet[j]<<endl;
    cout<<"<Delta_mag>[i]: "<<d.ave_dmag[j]<<"\t <blending>: "<<d.ave_bl[j]<<endl;}
    cout<<"============================================================"<<endl;
    cout<<"<Ds>: "<<d.ave_Ds<<"\t <Dl>: "<<d.ave_Dl<<endl;
    cout<<"<vl>: "<<d.ave_vl<<"\t <vs>: "<<d.ave_vs<<endl;
    cout<<"<Ml>: "<<d.ave_ml<<"\t <u0>: "<<d.ave_u0<<endl;
    cout<<"<Vt>: "<<d.ave_vt<<"\t<tE>: "<<d.ave_te<<"\t<color>: "<<d.ave_col<<endl;
    cout<<"blending[2]: "<<d.ave_bl[2]<<"\t apb[2]: "<<d.ave_apb[2]<<endl;
    cout<<"<RE>[AU]: "<<d.ave_re<<"\t effi_Microlensing: "<<d.nevent/nbb<<endl;
    cout<<"simulated stars: "<<d.nstar<<"\t detected events: "<<d.nevent<<endl;
    cout<<"No. lensing event[l,b]: "<<d.nfin<<"\t n_start[Omega_l] "<<s.nstart<<"\t nstart(OGLE): "<<s.nstarti<<endl;
    cout<<"<opt>: "<<d.ave_opt<<"\t <opt2>: "<<d.ave_opt2<<"\t event_ra: "<<d.ratio<<endl;
    cout<<">>>>>>>>>>>>>> No. microlensing total: "<<Nmicro<<endl;
    cout<<"fraction_red_gients[%]: "<<d.ave_fred<<"\t <I-band>: "<<d.aveI<<endl;
    cout<<"cadence: "<<cadence<<"\t nbb: "<<nbb<<endl;
    cout<<"================================================================================"<<endl;
    }
     }}///end of subfield
     fclose(_randStream);
     return(0);
}
////==========================================================================
double RandR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return(p*(up-down)+down);
}
///==========================================================================//
///                                                                          //
///                   Lensing                                                //
///                                                                          //
///==========================================================================//
void lensing(source & s, lens & l,detection & d,int ei)
{
    int sign, f;
    char filenam1[40];
    d.peak=0;
    int uue;
    double As,u,Delta1,seei,prob1,prob2,ggh,dela, magg;
    double maga,maga2,x;
    double m5[M]={0.0};
    int ds[M]={0};
    double chi2,chi1;

    d.t0=RandR(0.0 , Tobs);
    d.tmin=d.t0-3.5*l.tE;
    d.tmax=d.t0+3.5*l.tE;
    if(d.tmin<0.0)   d.tmin=0.0;
    if(d.tmax>=Tobs) d.tmax=Tobs;
    if(((d.tmin/year)-int(d.tmin/year))>double(7.0/12.0)) d.tmin=double(int(d.tmin/year)+1.0)*year;
    
    int Nb=int(1.0+(d.tmax-d.tmin)/cadence);
    if(d.tmin<=d.t0 and d.tmax>=d.t0)  d.peak=1;///it means we do see just one side of the light curve
    
    
    if(Nb<=4){
    d.det=0;
    chi1=chi2=d.dchi=0.0;}
    
    
    else{
    int flag[Nb];
    int numt=0;
    for(int i=0; i<M; ++i) ds[i]=0;
    d.det=0; chi1=chi2=d.dchi=0.0;
    for(int i=0;i<Nb; ++i)  flag[i]=0;
    int randd=int(RandR(0.0,NM-1));
    if(randd>=(NM-1)) randd=0;
    
    for(d.t=d.tmin; d.t<=d.tmax; d.t=d.t+cadence){
    randd+=1;
    if(randd>=(NM-1)) randd=0;
    prob1=RandR(0.0,100.0);///probability of bad weather
    prob2=RandR(0.0,100.0);///probability of uniform observation
    
    
    if(prob1>17.0 and prob2>=10.0){
    u=sqrt(s.u0*s.u0+(d.t-d.t0)*(d.t-d.t0)/l.tE/l.tE);
    As=(u*u+2.0)/(u*sqrt(u*u+4.0));
    f=d.filter[randd];
    if(f<0 or f>5){cout<<"Big error filter: "<<f<<"\t randd: "<<randd<<endl;  int yye; cin>>yye;}
    if(f<LL){////ugrizy
    
    maga=-2.5*log10(fabs(pow(10.0,-0.4*s.Map[f])*(As-1.0)+s.Fluxb[f]));
    if(maga>=satu[f] and maga<=detect[f]){
   // dels=RandN(0.2,1);///amount of flactuation in seeing
    //dela=RandN(1.0,2);///the amount of flactuation in magnification factor
      
//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(numt==0 and int(ei)%1000==0){
    FILE *test; 
    sprintf(filenam1,"./files/light/%c%.2lf%c%.2lf%c%d.dat",'l',s.lat,'_',s.lon,'_',ei);
    test=fopen(filenam1,"w");
    fclose(test);}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
  
    x=pow(10.0 , 0.4*(maga-d.m5c[randd]));
    Delta1=sqrt(fabs((0.04-gama[f])*x+gama[f]*x*x));///random photometric error    
    if(((0.04-gama[f])*x+gama[f]*x*x)<0.0){
    cout<<"BIG ERROR negative: "<<((0.04-gama[f])*x+gama[f]*x*x)<<endl; int yye; cin>>yye; }
    d.deltam[f]= Delta2 +Delta1 ;

    dela=fabs(RandN(1.0,1.0));
    if(numt%2==0) sign=-1;
    else          sign=+1; 
    maga2 = maga+sign*dela*d.deltam[f];
    chi1 += pow((maga2-s.magb[f])/d.deltam[f],2.0);
    chi2 += pow((maga2-maga )/d.deltam[f],2.0);
    
    if(d.det==0){
    if(fabs(maga2-s.magb[f])<5.0*d.deltam[f]) flag[numt]=0;
    else{
    if(maga2<s.magb[f])  flag[numt]=+1;//upper 
    else                 flag[numt]=-1;}//downer
    if(numt>=3 and abs(flag[numt]+flag[numt-1]+flag[numt-2]+flag[numt-3])==4) d.det=1;}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(int(ei)%1000==0){
    FILE *test; 
    sprintf(filenam1,"./files/light/%c%.2lf%c%.2lf%c%d.dat",'l',s.lat,'_',s.lon,'_',ei);
    test=fopen(filenam1,"a+");
    
    m5[2]=Cm[2]+d.Dcm[2]+0.5*(msky[2]-21.0)+d.m5c[randd]-Cm[f]-d.Dcm[f]-0.5*(msky[f]-21.0); 
    magg=-2.5*log10(fabs(pow(10.0,-0.4*s.Map[2])*(As-1.0)+s.Fluxb[2]));
    
    x=pow(10.0,+0.4*(magg-m5[2]));
    Delta1=sqrt((0.04-gama[2])*x+gama[2]*x*x);///random photometric error
    d.deltam[2]=sqrt(Delta2*Delta2+Delta1*Delta1);
 
    ggh= magg+fabs(dela)*d.deltam[2]*sign;
    fprintf(test,"%.3lf    %.3lf   %.4lf  %.4lf   %d  %.5lf\n",d.t/year,(d.t-d.t0)/l.tE,ggh,d.deltam[2],f, d.deltam[f]);  
    fclose(test);}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    ++numt;
    ++ds[f];}}}///end of probability of good weather
    if(((d.t/year)-int(d.t/year))>double(7.0/12.0)){d.t= double(int(d.t/year)+1.0)*year; }
    }///end of for time
   // cout<<"********** END OF time loop ************"<<endl;
    d.dchi= fabs(chi2-chi1);


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(int(ei)%1000==0){
    FILE *test; 
    sprintf(filenam1,"./files/light/%c%.2lf%c%.2lf%c%d.dat",'l',s.lat,'_',s.lon,'_',ei);
    test=fopen(filenam1,"a+");
    fprintf(test,"%.8lf   %.5lf  %.5lf  %.5lf %.5lf\n",s.u0,d.t0,l.tE,s.Map[2],s.magb[2]);
    for(int id=0; id<M; ++id)
    fprintf(test,"%d  %lf  %d  %d   %d\n",id,d.dchi,d.det,d.peak,ds[id]);
    fclose(test); } 
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    }///end of else 
  
    for (int i=0; i<M; ++i){
    x=pow(10.0,+0.4*(s.magb[i]-M50[i]));
    Delta1=sqrt(fabs((0.04-gama[i])*x+gama[i]*x*x));///random photometric error
    d.deltam[i]=sqrt(Delta2*Delta2+Delta1*Delta1);}
}
///==============================================================//6
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s)
{
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*Msun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opt=0.0;
    s.od_dlmc=s.od_hlmc=s.od_blmc=0.0; 
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC;
    s.od_dlmc +=  s.rho_dlmc[k]*x*(1.0-x)*dx*CC;
    s.od_blmc +=  s.rho_blmc[k]*x*(1.0-x)*dx*CC; 
    s.od_hlmc +=  s.rho_hlmc[k]*x*(1.0-x)*dx*CC;  }
    s.opt= fabs(s.od_disk+s.od_ThD+s.od_bulge+s.od_halo+ s.od_dlmc + s.od_blmc + s.od_hlmc );///total
    //cout<<"total_opticalD: "<<s.opd<<"\t od_disk: "<<s.od_disk<<endl;
   // cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
}
///==============================================================//
///                                                              //
///                  func_source   Initial amounts               //
///                                                              //
///==============================================================//
void func_source(source & s, CMD  &  cm)
{
    int struc,nums, num;
    double rho,rf,N_sblend[M+1];
    double Alv;
    double Ds,Ai[M],Av;
    double Map[M];      
    double extI; 
    double MIsource,Isource;   
    double a0[M+1]={0.9429,1.0138,0.94027, 0.8139, 0.6641, 0.5703,  0.7934};
    double b0[M+1]={1.9788, 0.5575,-0.2197, -0.4982, -0.6097, -0.5236, -0.5462};
    s.FIsource=0.0,s.Isource=0.0;
    double Mab[M]={0.0};   
    s.Ds=-10.0; 


    double maxnb=0.0;
    for(int i=0; i<(M+1); ++i){
    s.Fluxb[i]=0.0;
    N_sblend[i]=s.Nstart*pow(FWHM[i]*0.4,2)*M_PI/(3600.0*3600.0);
    N_sblend[i]=N_sblend[i]+RandN(sqrt(N_sblend[i]),1);
    if(N_sblend[i]<=1.0)   N_sblend[i]=1.0;
    if(N_sblend[i]>maxnb) maxnb=N_sblend[i];}
    
   
    for(int k=1; k<=int(maxnb+0.0034756346); ++k){
    
    do{
    nums=int(RandR(5.0,Num-5.0));
    rho=RandR(s.Romins , s.Romaxs); 
    Ds=double(nums*step);
    }while(rho>s.Rostari[nums] or Ds<0.0 or Ds>MaxD);///distance larger than 50.0
    if(Ds>MaxD or Ds<0.0){cout<<"ERROR (1): Ds: "<<Ds<<"\t MaxD: "<<MaxD<<"\t step: "<<step<<"\t num: "<<nums<<endl; int yye; cin>>yye;}
    
    rf= RandR(0.0, s.Rostar0[nums]); 
     if (rf<= s.rho_disk[nums])                    struc=0;///thin disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums])) struc=3;///halo
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums]+s.rho_hlmc[nums]))struc=4;///halo_lmc
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums]+s.rho_hlmc[nums]+s.rho_dlmc[nums]))struc=5;//disk_lmc
else struc=6;///bulge_lmc


    if(struc==0 or struc==5 ){///thin disk
    num=int(RandR(0.0 , N1-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_d[i][num];}
    if(k==1){ 
    s.logt=cm.logt_d[num]; 
    s.col = cm.col_d[num]; 
    s.logl=cm.logl_d[num];
    s.cl  =  cm.cl_d[num];}
    MIsource=cm.MI_d[num];}


    if(struc==1 or struc==6){///bulge
    num=int(RandR(0.0,N2-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_b[i][num];}
    if(k==1){ 
    s.logt= cm.logt_b[num];
    s.col=   cm.col_b[num]; 
    s.logl= cm.logl_b[num];
    s.cl=     cm.cl_b[num];}
    MIsource= cm.MI_b[num];}


    if(struc==2){///thick disk
    num=int(RandR(0.0,N3-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_t[i][num];}
    if(k==1){ 
    s.logt=cm.logt_t[num];
    s.col=  cm.col_t[num]; 
    s.logl=cm.logl_t[num];
    s.cl=    cm.cl_t[num]; }
    MIsource=cm.MI_t[num];}


    if(struc==3 or struc==4){///stellar halo
    num=int(RandR(0.0,N4-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_h[i][num];}
    if(k==1){ 
    s.logt=cm.logt_h[num];
    s.col=  cm.col_h[num]; 
    s.logl=cm.logl_h[num];
    s.cl=    cm.cl_h[num]; }
    MIsource=cm.MI_h[num];}
    
    Av= double(Ds*s.Avv/MaxD); 
    
    //Av=s.Avv;
    extI=-3.0; 
    Alv=fabs(a0[M]+b0[M]/Rv[struc]);///AIc/AV
    extI=double(Av*Alv+RandN(sigma[M],1) ); //extinction in Ic-band
    if(extI<0.0) extI=0.0;
    Isource= MIsource+5.0*log10(Ds*100.0)+extI; 
    
    if(Av<0.0 or Av==0.0 or Av>s.Avv or extI<0.0 or fabs(Alv-0.6)>0.1 or Rv[struc]>3.1 or Rv[struc]<2.5 ){
    cout<<"Error(2) Av:  "<<Av<<"\t Avv:  "<<s.Avv<<"\t extI:   "<<Av*Alv<<"\t Alv:  "<<Alv<<"\t Ds:   "<<Ds<<endl;  
    cout<<"Rv[struc]:   "<<Rv[struc]<<"\t Avv:  "<<s.Avv<<endl;   int uue;  cin>>uue;}




    for(int i=0; i<M; ++i){
    Alv=fabs(a0[i]+b0[i]/Rv[struc]);///Alambda/AV
    Ai[i]= Av*Alv+RandN(sigma[i],1); //extinction in other bands
    if(Ai[i]<0.0) Ai[i]=0.0;
    Map[i]=Mab[i]+5.0*log10(Ds*100.0)+Ai[i];
    if(N_sblend[i]>=k) s.Fluxb[i]+=fabs(pow(10.0,-0.4*Map[i] ));
    if(Alv<=0.0 or Av<0.0 or Ai[i]<0.0 or Map[i]<Mab[i] or Isource<MIsource or extI<0.0 or s.Fluxb[i]<0.0 or N_sblend[i]<1.0){
    cout<<"ERROR(3) Ds: "<<Ds<<"\t Alv: "<<Alv<<"\t Av: "<<Av<<"\t i:  "<<i<<endl;
    cout<<"A[i]: "<<Ai[i]<<"\t extI: "<<extI<<"\t Nblend: "<<N_sblend[i]<<endl; 
    cout<<"Map[i]: "<<Map[i]<<"\t Mab[i]: "<<Mab[i]<<"\t Fluxb:  "<<s.Fluxb[i]<<endl;
    cout<<"MI_abs: "<<MIsource<<"\t MI_app:  "<<Isource<<endl;
    int yue; cin>>yue;}}
    if(N_sblend[M]>=k) s.FIsource+=fabs(pow(10.0,-0.4*Isource));/// I-band Johnson standard filters


    if(k==1){
    s.struc=struc;   
    s.Ds=Ds;  
    s.nums=nums;
    s.AI=extI;
    s.Av=Av; 
    s.MIsource=MIsource; 
    s.Isource= Isource; 
    s.col= s.col+s.Av-s.AI;
    for(int i=0; i<M; ++i){s.Ai[i]=Ai[i];   s.Map[i]=Map[i];   s.Mab[i]=Mab[i];}}
    
    }///loop over the stars 


    s.blendI=double(pow(10.0,-0.4*s.Isource)/s.FIsource);
    for(int i=0; i<M; ++i){
    if(s.Fluxb[i]<=0.0){cout<<"BIG ERROR Flux is negative: "<<s.Fluxb[i]<<endl; int yye; cin>>yye; }
    s.magb[i]=-2.5*log10(s.Fluxb[i]);
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    s.nsbl[i]=double(N_sblend[i]);
    
    if(int(s.nsbl[i])<0.0 or s.nsbl[i]==0.0 or s.blend[i]>1.00001 or s.blend[i]<=0.00 or 
    (s.nsbl[i]==1.0 and s.blend[i]<1.0) or s.blendI>1.0 or s.blendI<0.0 or 
    fabs(s.Map[0]-s.Mab[0]-s.Ai[0]-s.Map[1]+s.Mab[1]+s.Ai[1])>0.1 or s.Ds<0.0  or s.Ds>MaxD){ 
    cout<<"N_sblend: "<<N_sblend[M]<<"\t Ds:  "<<s.Ds<<endl; 
    cout<<"blendI: "<<s.blendI<<"\t FIsource: "<<s.FIsource<<"\t Isource: "<<s.Isource<<endl;
    cout<<"BIGG ERRROR nsbl: "<<s.nsbl[i]<<"\t N_sblend: "<<N_sblend[i]<<"\t blend: "<<s.blend[i]<<endl; int uue; cin>>uue;}}

    s.u0=RandR(1.0e-50,1.0);
    s.magni=fabs((s.u0*s.u0+2.0)/(s.u0*sqrt(s.u0*s.u0+4.0)));
    s.Rs=s.logl-4.0*(s.logt-logT_sun);
    s.Rs=pow(10.0,s.Rs/2.0);/// in Sun radius 
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts              //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s)
{
    double f,test, tt;
    double rholens[s.nums+2]={0.0};
    
    l.rhomaxl=0.0;
    for(int k=1;k<int(s.nums-1) ;++k){
    rholens[k]=0.0;
    l.Dl=double(k*step);
    l.xls=l.Dl/s.Ds;
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}



    do{
    l.numl =int(RandR(1.0,s.nums-2.0));
    test =RandR(0.0,l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){
    cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;  int ue; cin>>ue;}
    }while(test>rholens[l.numl]);
    l.Dl=double(l.numl*step);
    
   
 double randflag=RandR(0.0,s.Rostar0[l.numl]);
       if (randflag<=fabs(s.rho_disk[l.numl]) ) l.struc=0;//disk
  else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1;//bulge
  else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2;//thick
  else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl]) ) l.struc=3;//halo
  else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl]+s.rho_hlmc[l.numl]))l.struc=4;//halo_lmc
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl]+s.rho_hlmc[l.numl]+s.rho_dlmc[l.numl])) l.struc=5;//disk_lmc
else if (randflag<=fabs( s.Rostar0[l.numl])) l.struc=6;//bar_lmc
else {  cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}



  if(l.struc==0 or l.struc==5){///thin disk
  do{
  l.Ml=RandR(0.08,4.5);
  test=RandR(0.0,57.0);
  if(l.Ml<=1.0) f=pow(l.Ml,-1.6);
  if(l.Ml>1.0) f=pow(l.Ml,-3.0);
  }while(test>f);}


  if(l.struc==1 or l.struc==6){///Galactic bulge
  do{
  l.Ml=RandR(0.08,1.4);
  test=RandR(0.3,378.3);
  f=pow(l.Ml,-2.35);
  }while(test>f);}


  if(l.struc==2){///thick disk
  do{
  l.Ml=RandR(0.08,1.4);
  test=RandR(0.8,3.8);
  f=pow(l.Ml,-0.5);
  }while(test>f);}


  if(l.struc==3 or l.struc==4){//stellar halo
  do{
  l.Ml=RandR(0.08,0.8);
  test=RandR(0.8, 3.8);
  f=pow(l.Ml,-0.5);
  }while(test>f);}



    l.Dl=double(l.numl*step);///kpc
    l.xls=l.Dl/s.Ds;
    l.RE=sqrt(4.0*G*l.Ml*Msun*s.Ds*KP)/velocity;
    l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
    vrel(s , l);
    l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0);///in day
    s.ro_star=s.Rs*Rsun*l.xls/l.RE; 


    if(l.tE<0.0 or l.tE==0.0 or l.RE<0.0 or  l.Dl>s.Ds  or l.xls>=1.0 or l.Dl<0.0 or s.Ds<0.0 or s.Ds>MaxD or l.Ml<0.08 or l.Ml>5.0){ 
    cout<<"Dl:  "<<l.Dl<<"\t Ds:   "<<s.Ds<<"\t xls:  "<<l.xls<<"\t Ml:   "<<l.Ml<<endl;
    cout<<"BIG ERROR te: "<<l.tE<<"\t RE: "<<l.RE<<"\t V_rel: "<<l.Vt<<endl;l.tE=0.1; int iie; cin>>iie;} 
   // cout<<"Dl: "<<l.Dl<<"\t xls: "<<l.xls<<"\t rho_star: "<<s.ro_star<<endl;
  //  cout<<"RE: "<<l.RE/AU<<"\t tE: "<<l.tE<<"\t Vt: "<<l.Vt<<endl;
  //  cout<<"************** End of lens function  ***************"<<endl; 
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm)
{

////mass log10(T) log10(age) log10(L) log10(g) metal U G R I Z Y Bj Vj Rj Ij CL TYP
    int yye;     
    int h,g,k1,k2,uui; 
    double mass, age, gravity,metal,Bjc,Vjc,Rjc,type; 
    char filename[40];
    FILE *fp2;
    double Age[YZ]={0.0}; double B[YZ]={0.0};  double M[YZ]={0.0};   double mm[YZ]={0.0}; 
    int number[70]={0};   int count[70]={0};   double Metal[70]={0.0}; 


    FILE *meta; 
    meta=fopen("./files/metal.txt","r"); 
    for(int i=0; i<70; ++i){
    fscanf(meta,"%lf   %d  %d\n",&Metal[i],&count[i],&number[i]);    
    if((Metal[i]<Metal[i-1] and i>0) or Metal[i]<0.0 or number[i]==0 or count[i]<0 or count[i]>YZ) {
    cout<<"ERROR Metal[i]: "<<Metal[i]<<"\t count[i]: "<<count[i]<<"\t number[i]: "<<number[i]<<"\t i: "<<i<<endl; cin>>uui;} }
    fclose(meta); 
    
    FILE *yz; 
    yz=fopen("./files/yzma.txt", "r"); 
    for(int i=0; i<YZ; ++i){
    fscanf(yz,"%lf   %lf   %lf  %lf\n",&Age[i],&mm[i],&B[i],&M[i]); 
    if(Age[i]<0.0 or mm[i]<0.0 or fabs(B[i])>1.5 or M[i]<0.5 or Age[i]>18.0) {   
    cout<<"ERROR Age: "<<Age[i]<<"\t metal: "<<mm[i]<<"\t B[i]"<<B[i]<<"\t M[i]: "<<M[i]<<"\t i: "<<i<<endl; cin>>uui; }}
    fclose(yz); 




////=================================== THIN DISK ============================== 
    int j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c%c.dat",'C','M','D','T','i','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTi.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %d %lf\n",
    &mass,&cm.logt_d[j],&age,&cm.logl_d[j],&gravity,&metal,&cm.Mab_d[0][j],
    &cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.Mab_d[3][j],&cm.Mab_d[4][j],&Bjc,&Vjc,&Rjc,&cm.MI_d[j],&cm.cl_d[j],&type);

///*******************************************************
    h=-1; 
    if(metal<Metal[0] or metal==Metal[0]) h=0; 
    else if(metal>Metal[69] or  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 or metal==Metal[i-1]) { h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui; } 
    //cout<<"metal: "<<metal<<"\t MEtal[h]: "<<Metal[h]<<"\t Metal[h+1]: "<<Metal[h+1]<<"\t h: "<<h<<endl;
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>3578 or mm[k1]!=mm[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else {
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl;
    int ooe; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  or  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>3577 or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_d[5][j]= double(B[g]+M[g]*cm.Mab_d[4][j]);  

    if(fabs(cm.Mab_d[5][j]-cm.Mab_d[4][j])>1.5 or fabs(age-Age[g])>3.0 or cm.Mab_d[5][j]==0.0){
    cout<<"ERROR: Mab_d(y-band): "<<cm.Mab_d[5][j]<<"\t Mab_d(z-band): "<<cm.Mab_d[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B[g]<<"\t M[g]: "<<M[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
///********************************************************
    cm.col_d[j]=Vjc-cm.MI_d[j]; 
    if(mass<0.0||cm.logt_d[j]<0.0||cm.Mab_d[2][j]>20.0||metal>0.12||age>10.0||cm.cl_d[j]>7||type>9.0){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N1){cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
   // cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;



////=================================== BULGE ================================== 
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c%d.dat",'C','M','D','B','b',2);
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDB.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %d %lf\n",
    &mass,&cm.logt_b[j],&age,&cm.logl_b[j],&gravity,&metal,&cm.Mab_b[0][j],
    &cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.Mab_b[3][j],&cm.Mab_b[4][j],&Bjc,&Vjc,&Rjc,&cm.MI_b[j],&cm.cl_b[j],&type);
   ///cm.Mab_b[5][j]=0.0; 
///*******************************************
    h=-1; 
    if(metal<Metal[0] or metal==Metal[0]) h=0; 
    else if(metal>Metal[69] or metal==Metal[69]) h=69;
    else { 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 or metal==Metal[i-1]) { h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui; } 
    //cout<<"metal: "<<metal<<"\t MEtal[h]: "<<Metal[h]<<"\t Metal[h+1]: "<<Metal[h+1]<<"\t h: "<<h<<endl;
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>3578 or mm[k1]!=mm[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else {
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl;
    int ooe; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  or  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>3577 or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_b[5][j]= double(B[g]+M[g]*cm.Mab_b[4][j]);  

    if(fabs(cm.Mab_b[5][j]-cm.Mab_b[4][j])>1.5 or fabs(age-Age[g])>3.0  or cm.Mab_b[5][j]==0.0){
    cout<<"ERROR:   Mab_b(y-band): "<<cm.Mab_b[5][j]<<"\t Mab_b(z-band): "<<cm.Mab_b[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B[g]<<"\t M[g]: "<<M[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
///***********************************************
    cm.col_b[j]=Vjc-cm.MI_b[j]; 
    if(mass<0.0||cm.logt_b[j]<0.0||Vjc>18.0||age>10.0||metal>0.9||cm.cl_b[j]>7||type>9.0){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
   // cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;



////=================================== THICK DISK =============================
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c%c.dat",'C','M','D','T','k','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTk.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %d %lf \n",
    &mass,&cm.logt_t[j],&age,&cm.logl_t[j],&gravity,&metal,&cm.Mab_t[0][j],
    &cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.Mab_t[3][j],&cm.Mab_t[4][j],&Bjc,&Vjc,&Rjc,&cm.MI_t[j],&cm.cl_t[j],&type);
///************************************************
    h=-1; 
    if(metal<Metal[0] or metal==Metal[0]) h=0; 
    else if(metal>Metal[69] or  metal==Metal[69]) h=69;
    else { 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 or metal==Metal[i-1]) { h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui; } 
    //cout<<"metal: "<<metal<<"\t MEtal[h]: "<<Metal[h]<<"\t Metal[h+1]: "<<Metal[h+1]<<"\t h: "<<h<<endl;
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>3578 or mm[k1]!=mm[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else {
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl;
    int ooe; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  or  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>3577 or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_t[5][j]=double(B[g]+M[g]*cm.Mab_t[4][j]);  

    if(fabs(cm.Mab_t[5][j]-cm.Mab_t[4][j])>1.5 or fabs(age-Age[g])>3.0 or cm.Mab_t[5][j]==0.0){
    cout<<"ERROR:   Mab_t(y-band): "<<cm.Mab_t[5][j]<<"\t Mab_t(z-band): "<<cm.Mab_t[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B[g]<<"\t M[g]: "<<M[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
///*************************************************
    cm.col_t[j]=Vjc-cm.MI_t[j]; 
    if(mass<0.0||cm.logt_t[j]<0.0||cm.Mab_t[2][j]>20.0||metal>0.025||cm.cl_t[j]>7|| type>9.0){
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
   // cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;






////=================================== STELLAR HALO =========================== 
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c.dat",'C','M','D','H','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDH.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &mass,&cm.logt_h[j],&age,&cm.logl_h[j],&gravity,&metal,&cm.Mab_h[0][j],
    &cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.Mab_h[3][j],&cm.Mab_h[4][j],&Bjc,&Vjc,&Rjc,&cm.MI_h[j],&cm.cl_h[j],&type);
///****************************************************************************
    h=-1; 
    if(metal<Metal[0] or metal==Metal[0]) h=0; 
    else if(metal>Metal[69] or  metal==Metal[69]) h=69;
    else { 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 or metal==Metal[i-1]) { h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui; } 
    //cout<<"metal: "<<metal<<"\t MEtal[h]: "<<Metal[h]<<"\t Metal[h+1]: "<<Metal[h+1]<<"\t h: "<<h<<endl;
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>3578 or mm[k1]!=mm[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else {
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl;
    int ooe; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  or  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>3577 or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_h[5][j]= double(B[g]+M[g]*cm.Mab_h[4][j]);  

    if(fabs(cm.Mab_h[5][j]-cm.Mab_h[4][j])>1.5 or fabs(age-Age[g])>3.0 or cm.Mab_h[5][j]==0.0){
    cout<<"ERROR:   Mab_h(y-band): "<<cm.Mab_h[5][j]<<"\t Mab_h(z-band): "<<cm.Mab_h[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B[g]<<"\t M[g]: "<<M[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
///*****************************************************************************
    cm.col_h[j]=Vjc-cm.MI_h[j]; 
    if(mass<0.0 or cm.logt_h[j]<0.0 or cm.Mab_h[2][j]>20.0 or metal>0.01 or cm.cl_h[j]>7|| type>9.0){
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
  //  cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
   cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
   
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s, int yy)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nnf=0.4/0.8;
   double mBarre;///stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   s.Romins=10000000000.0;
   double fd=1.0;///see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;///0.657066;////just stars brighter than V=11.5, but we change to consider all stars  
   double fh=1.0;///No limitation 
   double Rdd=2.17;///2.53;///2.17;
   double Rhh=1.33;///1.32;//1.33;

   double alm= 2.0; ////2KPC,  Mancini 2004
   double rlm_0= 1.76*0.01;/// solar mass/ Pc^-3
   double M_dlm=2.6;//[Msun]mass of disk LMC  from Mancini 2004 paper
   double M_blm=0.15*M_dlm;///[Msun] mass of bar LMC
   double xb0, yb0, zb0;///KPC
   double frac=0.05; // fraction of halo in the form of compact objects
   double qq=0.688; 
   double Rd0=1.54;
   double xol, yol, zol, x0, y0, z0; 
   double r0,  pos1, inc1, Rlm; 

  char filename[40];
  FILE *fill;     
  FILE *fil1;
  fil1=fopen("./files/density/right.txt", "a+"); 
  int flagf=0;
  if(fabs(s.decli+68.5)<0.8){
  flagf=1; 
  sprintf(filename,"./files/density/%c%d%c%d.dat",'d',0,'_',yy);
  fill=fopen(filename,"w");
   //if(!fill){cout<<"cannot open file longtitude : "<<s.lg<<"\t latitude: "<<s.bg<<endl;  exit(0);}
  }


for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;
   s.rho_hlmc[i]=s.rho_dlmc[i]=s.rho_blmc[i]=0.0; 


   x=double(i*step);
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = Dsun- x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);
   double rsun= sqrt(zb*zb+ yb*yb+ xb*xb); 


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///Msun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nnf)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*exp(nnf)*exp(-fabs(zb)/0.8)/(1.0+0.5*nnf);///Msun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=frac*(0.932*0.00001/867.067)*pow(0.5/Dsun,-2.44);
   else            s.rho_halo[i]=frac*(0.932*0.00001/867.067)*pow(rdi/Dsun,-2.44);///Msun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///Msun/pc^3
///==================================================================



///=====================  LMC density profiles ====================== 
   x0 = -x*cos(s.decli/RA)*sin(s.right/RA-RaLMC/RA);
   y0 =  x*sin(s.decli/RA)*cos(DecLMC/RA)-x*cos(s.decli/RA)*sin(DecLMC/RA)*cos(s.right/RA-RaLMC/RA);
   z0 = DLMC- x*cos(s.decli/RA)*cos(DecLMC/RA)*cos(s.right/RA-RaLMC/RA)-x*sin(s.decli/RA)*sin(DecLMC/RA);
   if(fabs(x-DLMC)<double(step*1.5)){s.xv=x0;  s.yv=y0;   s.zv=z0;}
   r0= sqrt(x0*x0+y0*y0+z0*z0); 
 
  
///=====================  LMC density HALO ======================
   if(r0<15.0) s.rho_hlmc[i]=fabs(rlm_0*frac/(1.0+r0*r0/alm/alm) );  ///Msun/PC^3
   else        s.rho_hlmc[i]=0.0;   
  
  
///=====================  LMC density Disk ======================
   pos1=(170.0-90.0)*M_PI/180.0; 
   inc1=34.7*M_PI/180.0;  
   xol= x0*cos(pos1)+ y0* sin(pos1);
   yol=-x0*sin(pos1)*cos(inc1) + y0*cos(pos1)*cos(inc1) -z0*sin(inc1);
   zol=-x0*sin(pos1)*sin(inc1) + y0*cos(pos1)*sin(inc1) +z0*cos(inc1);    
   Rlm= sqrt(xol*xol + yol*yol ); 
   double zd0=0.3;//KP Kim's paper
   double Rd1=1.8;///KPc  Kim's paper
   s.rho_dlmc[i]=fabs(M_dlm/(4.0*M_PI*zd0*Rd1*Rd1)*exp(-Rlm/Rd1)*exp(-fabs(zol/zd0))); ////kim [Msun/Pc^3]  Kim 2000
   

///=====================  LMC density Bulge ======================
   xb0= 1.2;  yb0=zb0= 0.44;
   pos1=(110.0-90.0)*M_PI/180.0; 
   inc1=0.0*M_PI/180.0;  
   xol= x0* cos(pos1) +  y0* sin(pos1);
   yol=-x0*sin(pos1)*cos(inc1) + y0*cos(pos1)*cos(inc1) -z0*sin(inc1);
   zol=-x0*sin(pos1)*sin(inc1) + y0*cos(pos1)*sin(inc1) +z0*cos(inc1);    
   s.rho_blmc[i]=M_blm*pow(2.0*M_PI,-1.5)/(xb0*yb0*zb0)*exp(-0.5*(xol*xol/xb0/xb0+yol*yol/yb0/yb0+zol*zol/zb0/zb0));///[Msun/Pc^3]  
///===============================================================================



s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]) +fabs(s.rho_hlmc[i])+fabs(s.rho_dlmc[i])+ fabs(s.rho_blmc[i]);///[Msun/pc^3]

s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Msun/deg^2]

s.Nstari[i]= binary_fraction*((s.rho_disk[i]+s.rho_dlmc[i])*fd/0.403445+s.rho_ThD[i]*fh/0.4542+(s.rho_halo[i]+s.rho_hlmc[i])*fh/0.4542+(s.rho_bulge[i]+s.rho_blmc[i])*fb/0.308571);////[Nt/pc^3] 

s.Nstari[i]=s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]

s.Nstart  +=  s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]

if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection
if(s.Rostari[i]<s.Romins) s.Romins=s.Rostari[i];///source selection

//cout<<"Rostari[i]:  "<<s.Rostari[i]<<"\t Ds:  "<<x<<endl;

if(flagf>0)
fprintf(fill,"%e   %e   %e   %e   %e  %e  %e   %e   %e   %e\n",
  x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.rho_dlmc[i],s.rho_blmc[i],s.rho_hlmc[i], s.Rostar0[i],s.Nstari[i]); 
 // cout<<"rho_disk(LMC):  "<<s.rho_dlmc[i]<<"\t rho_bar(LMC): "<<s.rho_blmc[i]<<endl;
  //cout<<"rho_halo(LMC):  "<<s.rho_hlmc[i]<<endl;
   }

if(flagf>0) fprintf(fil1,"%.5lf  %.5lf  %.5lf  %.5lf\n",s.right, s.decli, s.TET*RA, s.FI*RA);

//cout<<"Romins:  "<<s.Romins<<"\t Romaxs:  "<<s.Romaxs<<endl;

if(flagf>0)   {fclose(fill); fclose(fil1);}
//int rre; cin>>rre;

   //cout<<"xLMC:  "<<xLMC<<"\t yLMC:  "<<yLMC<<"\t zLMC:  "<<zLMC<<endl;
   //cout<<"LMC_distance:  "<<sqrt(xLMC*xLMC + yLMC*yLMC +  zLMC*zLMC)<<endl;
   //cout<<"xol:  "<<xol<<"\t yol:  "<<yol<<"\t  zol:  "<<zol<<endl;
   //cout<<"distance_LMC_Center:  "<<rlm<<"\t projected_distance:  "<<Rlm<<endl;
  // cout<<"Nstart [Nt/deg^2]: "<<s.Nstart<<"\t Ro_star [Mass/deg^2]: "<<s.Rostart<<endl;
   //cout<<">>>>>>>>>>>>>>>>>>>>>>>>> END OF DISK MODLE <<<<<<<<<<<<<<<<<<<<"<<endl;
   //exit(0);
   
}
///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
///===========================================================================//
double RandN(double sigma, double nnd){
   double rr,f,frand;
   do{
   rr=RandR(-sigma*nnd, sigma*nnd);//[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/sigma/sigma);
   frand=RandR(0.0,1.0);
   }while(frand>f);
   return rr;
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l)
{
if(l.Dl==0.0) l.Dl=0.00034735;
 double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*l.Dl*cos(s.TET)*cos(s.FI));
 double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*s.Ds*cos(s.TET)*cos(s.FI));
 if(Rlc==0.0) Rlc=0.00034346123;
 if(Rsc==0.0) Rsc=0.0004762654134;  
 double SVT, SVR, SVZ, LVT, LVR, LVZ,SVt,LVt;
 double SVb, SVl, LVb, LVl;
 double fv, testfv, tetd;
 double betal,betas,deltal,deltas,deltao;
 double SV_n1, LV_n1, VSun_n1, SVx, LVx, VSunx, SV_n2, LV_n2, VSun_n2,  vls1, vls2; 

 double NN=2.5;
 double VSunR =-10.3;
 double VSunT =vro_sun*(1.00762+0.00712)+6.3;
 double VSunZ = 5.9;
 double sigma_R_Disk= 43.0,  sigma_T_Disk= 27.8, sigma_Z_Disk=17.5;
 double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
 double sigma_R_halo= 131.0, sigma_T_halo=106.0, sigma_Z_halo=85.0;
 double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0, sigma_Z_Bulge=100.0;
 double Rho[8]={00.0}; double maxr=0.0;
 for(int i=0; i<8; ++i){  Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}
 
  double v_R_lmc=-57.0;
  double v_T_lmc=-226.0; 
  double v_Z_lmc= 221.0;
  double sigma_LMC=20.2; 
  double err_rlmc= 13.0; ///error of global velocity
  double err_tlmc= 15.0; 
  double err_zlmc= 19.0; 
 

  double test=RandR(0.0,maxr);//total ages
     if(test<=Rho[0])                  {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0;}
else if(test<=(Rho[0]+Rho[1]))         {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]))  {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3])) {sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4])) {sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5])) {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4;}
else if(test<=maxr)        {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5;}
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}  


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
/// Generate Source velocity components in Glactocenteric cylindrical coordinates(x',y')
    SVR=SVT=SVZ=0.0;
    if(s.struc==0){///Galactic disk
    SVR= RandN(sigma_R_Disk, NN); 
    SVT= RandN(sigma_T_Disk, NN);
    SVZ= RandN(sigma_Z_Disk, NN); 
    SVT =SVT+ vro_sun *(1.00762 * pow(Rsc/Dsun,0.0394) + 0.00712);}

   else if(s.struc==1){///Galactic bulge
   SVZ=RandN(sigma_Z_Bulge, NN); 
   SVR=RandN(sigma_R_Bulge, NN);
   SVT=RandN(sigma_T_Bulge, NN);}

   else if(s.struc==2){///thick disk
   SVR= RandN(sigma_R_TDisk, NN); 
   SVT= RandN(sigma_T_TDisk, NN); 
   SVZ= RandN(sigma_Z_TDisk, NN);
   SVT =SVT+ vro_sun *(1.00762*pow(Rsc/Dsun,0.0394) + 0.00712); }
   
   else if(s.struc==3){///stellar halo
   SVR= RandN(sigma_R_halo, NN); 
   SVT= RandN(sigma_T_halo, NN); 
   SVZ= RandN(sigma_Z_halo, NN);}
   
   else if(s.struc>3){
   SVR= RandN(sigma_LMC, NN); 
   SVT= RandN(sigma_LMC, NN); 
   SVZ= RandN(sigma_LMC, NN); 
   SVZ +=  v_Z_lmc +RandR(-err_zlmc, err_zlmc);
   SVR +=  v_R_lmc +RandR(-err_rlmc, err_rlmc);
   SVT +=  v_T_lmc +RandR(-err_tlmc, err_tlmc);}
   s.vs=sqrt(SVT*SVT+SVZ*SVZ+SVR*SVR);

///======================================================================================
/// Generate Lens velocity components in Glactocenteric cylindrical coordinates(x',y'
    LVR=LVT=LVZ=0.0;
    if(l.struc==0){///Galactic disk
    LVR= RandN(sigma_R_Disk, NN); 
    LVT= RandN(sigma_T_Disk, NN);
    LVZ= RandN(sigma_Z_Disk, NN); 
    LVT =LVT+ vro_sun *(1.00762 * pow(Rlc/Dsun,0.0394) + 0.00712);}

   else if(l.struc==1){///Galactic bulge
   LVZ=RandN(sigma_Z_Bulge, NN); 
   LVR=RandN(sigma_R_Bulge, NN);
   LVT=RandN(sigma_T_Bulge, NN);}

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN); 
   LVT= RandN(sigma_T_TDisk, NN); 
   LVZ= RandN(sigma_Z_TDisk, NN);
   LVT =LVT+ vro_sun *(1.00762*pow(Rlc/Dsun,0.0394) + 0.00712); }
   
   else if(l.struc==3){///stellar halo
   LVR= RandN(sigma_R_halo, NN); 
   LVT= RandN(sigma_T_halo, NN); 
   LVZ= RandN(sigma_Z_halo, NN);}
   
   else if(l.struc>3){
   LVR= RandN(sigma_LMC, NN); 
   LVT= RandN(sigma_LMC, NN); 
   LVZ= RandN(sigma_LMC, NN); 
   LVZ +=  v_Z_lmc +RandR(-err_zlmc, err_zlmc);
   LVR +=  v_R_lmc +RandR(-err_rlmc, err_rlmc);
   LVT +=  v_T_lmc +RandR(-err_tlmc, err_tlmc); }
   l.vl=sqrt(LVT*LVT+LVZ*LVZ+LVR*LVR);
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHH BETA HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

   betal=betas=0.0;
   tetd= s.TET;
   test= double(l.Dl*cos(s.FI)*sin(tetd)/Rlc);
   if(fabs(test-1.0)<0.01)       betal= pi/2.0;
   else if(fabs(test+1.0)<0.01)  betal=-pi/2.0;
   else                          betal=asin(test);

   test= double(s.Ds*cos(s.FI)*sin(tetd)/Rsc); 
   if( fabs(test-1.0)<0.01)     betas=pi/2.0;
   else if(fabs(test+1.0)<0.01) betas=-pi/2.0;
   else                         betas=asin(test);
    
   if(Dsun < fabs(l.Dl*cos(s.FI)*cos(tetd))) betal= pi-betal; 
   if(Dsun < fabs(s.Ds*cos(s.FI)*cos(tetd))) betas= pi-betas; 

   if(fabs(l.Dl*cos(s.FI)*sin(tetd))>double(Rlc*1.0003256463) or fabs(test)>1.0063645){
   cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
   cout<<"FI: "<<s.FI<<"\t TET: "<<tetd<<"\t betal:  "<<betal<<endl;
   cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<"\t betas:   "<<betas<<endl;
   cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(tetd)/Rlc<<"\t sin(s): "<<test<<endl;
   int ew; cin>>ew;}

///HHHHHHHHHHHHHHHHHHHHHHHHHH  DELTA   HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.TET>pi)  tetd=s.TET-2.0*pi; 
    deltal= pi - fabs(tetd) -fabs(betal);
    deltas= pi - fabs(tetd) -fabs(betas);  
    if(betal<0.0)  deltal= -1.0*deltal;
    if(betas<0.0)  deltas= -1.0*deltas;
    deltao= pi-fabs(tetd);
    if(tetd<0.0)  deltao=-1.0*deltao; 
 
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH    

    
    SV_n1 =+SVR * sin(deltas)- SVT * cos(deltas);
    LV_n1 =+LVR * sin(deltal)- LVT * cos(deltal);
    VSun_n1=+VSunR*sin(deltao)-VSunT*cos(deltao);
    
    SVx= -SVR*cos(deltas)- SVT*sin(deltas);
    LVx= -LVR*cos(deltal)- LVT*sin(deltal);
    VSunx= -VSunR*cos(deltao) -VSunT*sin(deltao);
    
    SV_n2=-sin(s.FI)*(SVx) + cos(s.FI)*SVZ;
    LV_n2=-sin(s.FI)*(LVx) + cos(s.FI)*LVZ;
    VSun_n2=-sin(s.FI)*(VSunx)+cos(s.FI)*(VSunZ);
   
    vls1= l.xls*SV_n1 - LV_n1 +(1.0-l.xls)* VSun_n1;  ///Source - lens 
    vls2= l.xls*SV_n2 - LV_n2 +(1.0-l.xls)* VSun_n2;  /// Source -lens
    l.Vt=sqrt(fabs( vls1*vls1 + vls2*vls2 ) );
   
    if (l.Vt<0.0 or l.Vt>1.0e6){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<s.vs<<endl;   int yee; cin>>yee;}
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
