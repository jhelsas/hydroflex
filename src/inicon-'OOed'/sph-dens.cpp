#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include "splitandfit.h"
#include "trial-functions.h"

using namespace std;

#define _USE_MATH_DEFINES

int main(){
  int N=0,D=0;
  int i,j,k,l,err,Npoints;
  double p[5];
  double *xp,*x,*u,*S,s,dist,h,a,b,c;
  ifstream sphfile;
  ofstream plotfile;
  FILE *sphofile;
  
  p[0]=1.; 
  p[1]=0.;
  p[2]=0.;
  p[3]=.5;
  p[4]=.5;
  
  h=0.1;    
  
  cout << "oi-1\n";
    
  err=sph_read("SPH-particles.dat",&D,&N,&x,&u,&S);if(err!=0) return err;
  double xl[D],xu[D],dx[D];
  for(l=0;l<D;l+=1){xl[l]=-1.;dx[l]=0.05;xu[l]=1.+1.01*dx[l];}
  
  cout << "oi-2\n";
  
  err=create_grid(D,&xp,xl,xu,dx,&Npoints);if(err!=0) return err;         
  
  cout << "oi-3\n";
  err=sph_dens(D,N,Npoints,xp,x,S,h,xl,xu,winicon,"ploting.dat",p);if(err!=0){ cout << err << endl; return err;}
  
  delete x;delete u; delete S;
  
  return 0;
}
