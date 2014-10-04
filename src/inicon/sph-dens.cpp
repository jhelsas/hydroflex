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
  double p[3];
  double *xp,*x,*u,*S,s,dist,h,a,b,c;
  ifstream sphfile;
  ofstream plotfile;
  FILE *sphofile;
  
  p[0]=1.0; /*  s0 */ 
  p[1]=1.0; /*  q  */
  p[2]=1.0; /* tau */ 
  
  h=0.1;    
    
  err=sph_read(argv[1],&D,&N,&x,&u,&S);if(err!=0) return err;

  double xl[D],xu[D],dx[D];
  for(l=0;l<D;l+=1){xl[l]=-7.0;dx[l]=0.05;xu[l]=7.0+1.01*dx[l];}
  
  err=create_grid(D,&xp,xl,xu,dx,&Npoints);if(err!=0) return err;         
  err=sph_dens(D,N,Npoints,xp,x,S,h,xl,xu,gubser_entropy,"ploting.dat",p);if(err!=0){ cout << err << endl; return err;}
  
  delete x;delete u; delete S;
  
  return 0;
}
