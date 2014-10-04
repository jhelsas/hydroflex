#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#define MYPI 3.1415926
/*#define A_d 1.*/  /*dim = 1*/
/*#define A_d (15.)/(7.*MYPI)*/ /*dim=2*/
#define A_d (3.0)/(2.0*MYPI) /*dim=3*/

#define DIM 3

using namespace std;

typedef struct PHSDparticle{
  int phsdid, tid;
  double x[4],p[4];
} PHSDparticle;

double signum(double z)
{
	if(z>=0) return 1.0;
	else return (-1.0);
}

double w_bspline(double r,double h)
{
	double R;
	if(h<=0.) exit(10);
	R=fabs(r)/h;
	if(R>=2.)
		return 0;
	else if((1.<=R)&&(R<2.))
		return (A_d)*(1./6.)*(2.-R)*(2.-R)*(2.-R)/pow(h,DIM);
	else
		return ((A_d)*((2./3.)-(R*R) + (R*R*R/2.0)))/pow(h,DIM);
}

int print_Tmn(int D,double **Tmn){
  int i,j;
  cout << "\n";
  for(i=0;i<=D;i+=1){
    for(j=0;j<=D;j+=1)
      cout << " " << Tmn[i][j] << " ";
    cout << "\n";
  }
  cout << "\n";
  return 0;
}

int show_particle_content(int D,vector <PHSDparticle> PP){
  int j;
  vector <PHSDparticle> :: iterator pp;
  for(pp=PP.begin();pp!=(PP.begin()+10);pp++){
    for(j=0;j<=D;j+=1)
      cout << pp->x[j] << " ";
    for(j=0;j<=D;j+=1)
      cout << pp->p[j] << " ";
    cout << pp->phsdid << " " << pp->tid << "\n";
  }
  
  return 0;
}

int main(){
  const int D=3;
  int i,j,k,l,m,n,N,isGood,err,direct,dTmn;
  int tid,phsdid,Nsteps,step,smoothing,show;
  double dist,t,dx[D+1],x[D+1],p[D+1],ti,h,t0=1.422649;
  double g[D+1][D+1],T[D+1][D+1],tpmn[D+1][D+1],tg[(D+1)*(D+1)],**Tprint;
  PHSDparticle Pt;
  vector <int> si,sf;
  vector <int> :: iterator sit;
  vector <PHSDparticle> PP;
  vector <PHSDparticle> :: iterator pp;
  //ifstream phsdin("fort.5001");
  ifstream phsdin("fort.5001");
  gsl_matrix_view m_gsl;  
  gsl_vector_complex *eval = gsl_vector_complex_alloc (D+1);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (D+1, D+1);
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (D+1);
  gsl_complex eval_k;
  gsl_vector_complex_view evec_k;
  
  smoothing=0;
  direct=0;
  show=0;
  h=0.5;
  
  for(i=0;i<=D;i+=1)
    for(j=0;j<=D;j+=1)
      g[i][j]=0.0;

  g[0][0]=1.0;
  for(i=1;i<=D;i+=1)
    g[i][i] = -1.0;
  
  dTmn=-1;
  if(direct==1)
    goto end;

  N=1709502;
  //N=1602865;
  //N=1;
  
  if(phsdin.is_open()){
    for(i=0;i<N;i+=1){
      phsdin >> x[0] >> x[1] >> x[2] >> x[3];
      phsdin >> p[0] >> p[1] >> p[2] >> p[3];
      phsdin >> phsdid >> tid;
            
      Pt.x[0] = x[0];Pt.x[1] = x[1];Pt.x[2] = x[2];Pt.x[3] = x[3];
      Pt.p[0] = p[0];Pt.p[1] = p[1];Pt.p[2] = p[2];Pt.p[3] = p[3];
      Pt.phsdid=phsdid;Pt.tid=tid;
      
      PP.push_back(Pt);
    }
  }
  
  if(show==0)
    err=show_particle_content(D,PP);
  
  i=0;
  ti=PP[0].x[0];
  si.push_back(i);
  for(pp=PP.begin();pp!=PP.end();pp++){
    if(pp->x[0]!=ti){
      si.push_back(i);
      sf.push_back(i);
      ti=pp->x[0];
    }
    i+=1;
  }
  sf.push_back(N);
  
  if(si.size()!=sf.size()){
    cout << "problemas! o tamanho dos limitantes nao batem\n";
    return 1;
  }
  
  Nsteps=si.size();
  for(i=0;i<Nsteps;i+=1){
    cout << "i: " << i << " ";
    cout << "si,sf:" << si[i] << " " << sf[i];
    cout << " t: "<< PP[si[i]].x[0] << " Npart: " << (sf[i]-si[i]) << "\n";
  }  
  
  cout << "\nInsira o momento desejado \n";
  cin >> x[0]; 
  cout << "\nInsira a posicao desejada \n";
  cin >> x[1] >> x[2] >> x[3];
  cout << "\n\n";

  dx[0]=0.0045803659595549107;
  dx[1]=0.5;
  dx[2]=0.5;
  //dx[3]=0.0091607319191098213/2.0;
  dx[3]=0.043428700417280197/2.0;
  step=0;

  for(i=0;i<Nsteps;i+=1){
    if(PP[si[i]].x[0]>x[0])
      break;
  }
  step=i-1;
  cout << "\nx[0]= " << x[0] << "; step = " << step << ", t_step= " << PP[si[step]].x[0] << "\n\n";
  
  for(m=0;m<=D;m+=1)
    for(n=0;n<=D;n+=1)
      T[m][n]=0.0;
  
  cout << "vcell: " << 2.0*2.0*2.0*dx[1]*dx[2]*dx[3] << "\n";

  Tprint = new double*[D+1];
  for(i = 0; i <=D; i+=1)
    Tprint[i] = new double[D+1];

  for(i=si[step];i<sf[step];i+=1){
    isGood=0;
    for(l=1;l<=D;l+=1){
      dist=fabs(PP[i].x[l]-x[l]);
      if(dist>=dx[l])
        isGood=1;
    }
    if(isGood!=0)
      continue;
    cout << i - si[step] << ": " << PP[i].p[0] << " " << PP[i].p[1] << " " << PP[i].p[2] << " " << PP[i].p[3] << " \n";
   
    if(smoothing==0){
      for(l=0;l<=D;l+=1)
        for(m=0;m<=D;m+=1)
          Tprint[l][m]=(((PP[i].p[l])*(PP[i].p[m]))/PP[i].p[0])*(1.0/(2.0*2.0*2.0*dx[1]*dx[2]*dx[3]));
    
      err=print_Tmn(D,Tprint);

      for(l=0;l<=D;l+=1)
        for(m=0;m<=D;m+=1)
          T[l][m]+=(((PP[i].p[l])*(PP[i].p[m]))/PP[i].p[0])*(1.0/(2.0*2.0*2.0*dx[1]*dx[2]*dx[3]));
    }
    else if(smoothing==1){
      dist=0.0;
      for(l=1;l<=D;l+=1)
        dist+=(x[l]-PP[i].x[l])*(x[l]-PP[i].x[l]);
      dist=sqrt(dist);
    
      for(l=0;l<=D;l+=1)
        for(m=0;m<=D;m+=1)
          T[l][m]+=(((PP[i].p[l])*(PP[i].p[m]))/PP[i].p[0])*w_gauss(dist,h);
    }
  }

  for(i=0;i<=D;i+=1)
    for(j=0;j<=D;j+=1)
      Tprint[i][j]=T[i][j];

  isGood=print_Tmn(D,Tprint); 
  
  end:
  if(direct==1){
    if(dTmn==0){
      T[0][0]=11638.87060082359857915434986;T[1][0]= 13.62797602302641664095972;T[2][0]= 11.25653370772096728558154;T[3][0]=-5819.09262053193469910183921;  
      T[0][1]=   13.62797602302641664095972;T[1][1]=  0.06382809358531429944072;T[2][1]=  0.02636059403570561662788;T[3][1]=  -13.62717353396256925179841;
      T[0][2]=   11.25653370772096728558154;T[1][2]=  0.02636059403570561662788;T[2][2]=  0.04354702632542096357726;T[3][2]=  -11.25587086202898667863792;
      T[0][3]=-5819.09262053193469910183921;T[1][3]=-13.62717353396256925179841;T[2][3]=-11.25587086202898667863792;T[3][3]=11637.49992166183801600709558;
    }else if(dTmn==1){
      T[0][0]=13759.4912109375000000000000;T[1][0]= 15.96470451354980468750000;T[2][0]=-50.63423919677734375000000;T[3][0]= -13758.85351562500000000000000; 
      T[0][1]=15.96470451354980468750000  ;T[1][1]=  0.01852334477007389068604;T[2][1]= -0.05874931439757347106934;T[3][1]=    -15.96396446228027343750000; 
      T[0][2]=-50.63423919677734375000000 ;T[1][2]= -0.05874931439757347106934;T[2][2]=  0.18633146584033966064453;T[3][2]=     50.63189697265625000000000; 
      T[0][3]=-13758.853515625000000000000;T[1][3]=-15.96396446228027343750000;T[2][3]= 50.63189697265625000000000;T[3][3]=  13758.21679687500000000000000;
    }else if(dTmn==2){
       T[0][0]=13759.49112707040512759704143  ;T[1][0]=   15.96470447148223570366099  ;T[2][0]=  -50.63423996506128332839580 ;T[3][0]= -13758.85408746870416507590562; 
       T[0][1]=  15.96470447148223570366099   ;T[1][1]=   0.01852334410538849368555   ;T[2][1]=  -0.05874931490670832340273  ;T[3][1]=  -15.96396533448335475213753  ;
       T[0][2]=  -50.63423996506128332839580  ;T[1][2]=   -0.05874931490670832340273  ;T[2][2]=    0.18633147353795234679730 ;T[3][2]=    50.63189569115159827106254 ;
       T[0][3]=-13758.85408746870416507590562 ;T[1][3]=   -15.96396533448335475213753 ;T[2][3]=    50.63189569115159827106254;T[3][3]=  13758.21707736078860762063414;
    }else if(dTmn==3){
      T[0][0]=  679.10276719510454768169438  ;T[1][0]=    1.77702177683438478084099  ;T[2][0]=    33.76839227820781985656140 ;T[3][0]=  -677.10968467056159170169849;  
      T[0][1]=  1.77702177683438478084099    ;T[1][1]=  0.00464996838164319272607    ;T[2][1]=   0.08836242663965095112122   ;T[3][1]=  -1.77180643797813286433041  ;
      T[0][2]= 33.76839227820781985656140    ;T[1][2]=  0.08836242663965095112122    ;T[2][2]=   1.67913366303117861377814   ;T[3][2]= -33.66928622860418585105435  ;
      T[0][3]= -677.10968467056159170169849  ;T[1][3]=   -1.77180643797813286433041  ;T[2][3]=   -33.66928622860418585105435 ;T[3][3]=   675.12245159641349800949683;
    }
  }

  for(l=0;l<=D;l+=1)
    for(j=0;j<=D;j+=1){
      tpmn[l][j]=0.0;
      for(k=0;k<=D;k+=1)
        tpmn[l][j] += g[l][k]*T[k][j];
    }
  
  for(l=0;l<=D;l+=1)
    for(j=0;j<=D;j+=1)
      tg[l+(D+1)*j]=tpmn[l][j];
        
  m_gsl=gsl_matrix_view_array(tg,D+1,D+1);  
        
  gsl_eigen_nonsymmv (&m_gsl.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_sort (eval, evec,GSL_EIGEN_SORT_ABS_DESC);
        
  for(k=0;k<=D;k+=1){
    eval_k = gsl_vector_complex_get (eval,k);
    evec_k = gsl_matrix_complex_column (evec,k);
	  
    cout << GSL_REAL(eval_k) << "\n";
  }
    
  gsl_eigen_nonsymmv_free(w);
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
  
  PP.clear();
  return 0;
}
