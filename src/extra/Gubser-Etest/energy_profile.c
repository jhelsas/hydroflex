#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_integration.h>

double gubser_transverse_energy(double r,void *p){
  double Q=*(double*)p,foo,lambda;
  foo = r*(3.*(1+Q*Q*(1+r*r))*(1+Q*Q*(1+r*r)) + 4.*Q*Q*Q*Q*r*r);
  lambda = 1+2.*Q*Q*(1+r*r)+Q*Q*Q*Q*(1-r*r)*(1-r*r);
  foo = foo/cbrt(gsl_sf_pow_int(lambda,7));
  return foo;
}

int main(){
  int err;
  size_t limit;
  double tau,tau0,dt,tauM;
  double *p,res,abserr;
  double eabs,erel,I,q=1.0,Q;
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  FILE *dadosout;
  
  F.function = &gubser_transverse_energy;
  
  tau0=1.0;
  tauM=8.0;
  dt=0.01;
  eabs=0.0001;
  erel=1e-7;
  limit =1000;
  dadosout=fopen("energy_profile.dat","w");
  for(tau=tau0;tau<tauM;tau+=dt){
	Q=q*tau;
	F.params = &Q;
    err=gsl_integration_qagiu(&F,0.0,eabs,erel,limit,w,&I,&abserr);

    fprintf(dadosout,"%.8f %.8f %.8f %d %lf %lf\n",tau,I,abserr,w->size,log(tau),log(I));
  }
  
  fclose(dadosout);
  gsl_integration_workspace_free (w);

  return 0;
} 
