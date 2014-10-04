#include <iostream>
#include <cmath>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TF1.h"
#include "TApplication.h"

using namespace std;

const int D=1;

typedef struct SPHparticle{
  int good;
  double x[D],S;
} SPHparticle;

double f(double *x){
  return exp(-x[0]*x[0])/sqrt(acos(-1.));
}


#define DIM 1

#define A_d 1.  
// #define A_d (15.)/(7.*M_PI) 
// #define A_d (3.0)/(2.0*MYPI) 

double w_bspline(double r,double h)
{
	double R;
	if(h<=0.) return -1.;
	R=fabs(r)/h;
	if(R>=2.)
		return 0;
	else if((1.<=R)&&(R<2.))
		return (A_d)*(1./6.)*(2.-R)*(2.-R)*(2.-R)/pow(h,DIM);
	else
		return ((A_d)*((2./3.)-(R*R) + (R*R*R/2.0)))/pow(h,DIM);
}

int msplit(){
  vector <SPHparticle> lsph;
  int N=101,i,j,k,l,m;
  double x[D],xl[D],xu[D],dx[D],S,min,max;
  double h=0.1,kh=0.1*h,dS,s,dist=0.;
   
  xl[0]=-3.0; xu[0]=3.0; dx[0]=(xu[0]-xl[0])/((double)(N-1));

  TH1D *horig = new TH1D("horig","",101,-3.,3.);

  for(i=0;i<N;i+=1){
    for(l=0;l<D;l+=1)
      x[l]=xl[l]+((double)i)*dx[l];
    S=f(x);
    {
      SPHparticle sphp;
      for(l=0;l<D;l+=1)
        sphp.x[l]=x[l];
      sphp.S=S;
      sphp .good=0;
      lsph.push_back(sphp);
      horig->Fill(sphp.x[0],sphp.S);
    }
  }

  for(i=0;i<lsph.size();i+=1){
    if(i==0){min=lsph[i].S;max=lsph[i].S;}
    else{
      if(lsph[i].S<min)
        min=lsph[i].S;
      if(lsph[i].S>max)
        max=lsph[i].S;
    }
  }

  //min=0.001*max;

  for(i=0;i<lsph.size();i+=1){
    m=(int)((lsph[i].S)/min);
    dS=lsph[i].S-((double)m)*min;
    if(m<=1)
      continue;
    
    lsph[i].good=1;
    
    for(j=0;j<m;j+=1){
      for(l=0;l<D;l+=1)
        x[l]= (lsph[i].x[l]-kh)+((2.0*kh)/((double)(m-1)))*((double)j);
      
      {
        SPHparticle sphp;
        for(l=0;l<D;l+=1)
          sphp.x[l]=x[l];
        sphp.S = min+(dS/((double)m));
        sphp.good=0;
        lsph.push_back(sphp);
      }
    }
  }

  for(i=0;i<lsph.size();i+=1)
    if(lsph[i].good!=0){lsph.erase(lsph.begin()+i);i-=1;}

  for(i=0;i<lsph.size();i+=1)
    if(lsph[i].good!=0) cout<< "problemas na limpeza\n";

  s=0.;
  for(i=0;i<lsph.size();i+=1)
    s+=lsph[i].S;

  for(i=0;i<lsph.size();i+=1)
    lsph[i].S=(lsph[i].S)/s;

  TGraph *gr1 = new TGraph();
  N=50;
  xl[0]=-3.0; xu[0]=3.0; dx[0]=(xu[0]-xl[0])/((double)(N-1));
  for(i=0;i<N;i+=1){
    s=0.;
    for(l=0;l<D;l+=1)
      x[l]=xl[l]+((double)i)*dx[l];
    for(j=0;j<lsph.size();j+=1){
      dist=0.;
      for(l=0;l<D;l+=1)
        dist+=(x[l]-lsph[j].x[l])*(x[l]-lsph[j].x[l]);
      dist=sqrt(dist);

      s += (lsph[j].S)*w_bspline(dist,h);
    }
    gr1->SetPoint(i,x[0], s);
  }

  TGraph *gr = new TGraph();

  TH1D *hnew = new TH1D("hnew","",101,-3.,3.);
  for(i=0;i<lsph.size();i+=1) {
    gr->SetPoint(i,lsph[i].x[0],lsph[i].S);
    hnew->Fill(lsph[i].x[0],lsph[i].S);
  }

  horig->SetLineColor(2);
  horig->Draw();

  hnew->SetLineColor(4);
  hnew->Draw("same");
  
//  TF1 *f1 = new TF1("f1","exp(-x[0]*x[0])+0.5*exp(-pow((x[0]-1.)/0.2,2.))",-3.,3.);
  TF1 *f1 = new TF1("f1","exp(-x[0]*x[0])/sqrt(acos(-1))",-3.,3.);
  f1->Draw("same");
 
  TH1D *hdiff = (TH1D*)horig->Clone("hdiff");
  hdiff->Add(hnew,-1.);
  TCanvas *c2 = new TCanvas("c2","c2");
  hdiff->Draw();

  TCanvas *c3 = new TCanvas("c3","c3");
//  gr->SetMarkerStyle(24);
//  gr->Draw("ap");
  gr1->SetMarkerStyle(24);
  gr1->SetMarkerColor(4);
  gr1->Draw("ap");
  f1->Draw("same");

  return 0;
}

int main(int argc,char* argv[]) {

  TApplication app("app",&argc,argv);
  msplit();

  app.Run();

  return 0;
}
