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

int main(int argc,char **argv){
  int N=0,D=0;
  int i,j,k,l,err,Npoints;
  double p[3],tau=1.0;
  double *xp,*x,*u,*S,s,dist,h,a,b,c;
  ifstream sphfile;
  ofstream plotfile;
  char infilename[100+1],outfilename[100+1];
  FILE *sphofile;
  
  p[0]=1.0; /*  s0 */ 
  p[1]=1.0; /*  q  */
  p[2]=1.0; /* tau */ 
  
  h=0.1;    
  
  double xl[D],xu[D],dx[D];
  for(l=0;l<D;l+=1){xl[l]=-3.0;dx[l]=0.15;xu[l]=3.0+1.01*dx[l];}  
  err=create_grid(D,&xp,xl,xu,dx,&Npoints);if(err!=0){cout << "erro = " << err << endl; return err;}
  
  for(tau=1.0;tau<8.0;tau+=1.0){
	p[2]=tau;
	
	cout << "tau= " << tau << endl;
	
	sprintf(infilename,"results/gubser-data/eqp-gubser-(t=%.1lf).dat",tau);
    err=sph_read(infilename,&D,&N,&x,&u,&S);if(err!=0) return err;
    
    sprintf(outfilename,"results/gubser-plots/dens-gubser-(t=%.1lf).dat",tau);
    err=sph_density_ploting(D,tau,Npoints,xp,xl,xu,N,x,u,S,h,gubser_entropy, gubser_proper_entropy,gubser_proper_energy,outfilename,p);
  }
  
  delete x;delete u; delete S;
  
  return 0;
}
