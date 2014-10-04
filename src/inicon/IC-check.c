#include <stdio.h>

double gubser_entropy(double *x,size_t dim,void *par){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 2) {
    fprintf(stderr, "error: dim != 2");
    abort();
  }
  double *p = (double*)(lpar->p);
  double s0=p[0], q=p[1],tau=p[2];
  double r2,lambda;
  r2=x[0]*x[0]+x[1]*x[1];
  lambda = 1.+2.*q*q*(tau*tau+r2) + q*q*q*q*(tau*tau-r2)*(tau*tau-r2);
  
  return (s0*(2.*q)*(2.*q)*(1.+q*q*(tau*tau+r2)))/(lambda*sqrt(lambda));
  
}

int gubser_velocity(double *x,size_t dim,void *par,double *u){
  if(dim!=2)
    return 1;
    
  int err;
  double *p = (double*)(par);
  double s0=p[0], q=p[1],tau=p[2];
  double r2,lambda,c;
  r2=x[0]*x[0]+x[1]*x[1];
  lambda = 1.+2.*q*q*(tau*tau+r2) + q*q*q*q*(tau*tau-r2)*(tau*tau-r2);
  c=(2.0*q*q*tau)/sqrt(lambda);
  u[0]=c*x[0];
  u[1]=c*x[1];
  
  return 0;
}

int main(int argc,char **argv){
  const int D=2;
  int i,N,err;
  double nu,q,S,x[2],u[2],err,U[2],p[3],St;
  FILE *dadosin,*dadosout;
  if(argc!=3)
    return 1;
  
  N=atoi(argv[2]);
  dadosin = fopen(argv[1],"r");
  dadosout= fopen("check.dat","w");
  p[0]=1.0;p[1]=1.0;p[2]=1.0;
  St=0.0;
  for(i=0;i<N;i+=1){
    fscanf(dadosin,"%lg %lg %lg %lg %lg %lg %lg",&nu,&q,&S,&(x[0]),&(x[1]),&(u[0]),&(u[1]));
    St+=S;
    err=gubser_velocity(x,D,p,U);
    fprintf(dadosout,"%lf %lf %lf %lf",);
  }
  
  return 0;
}
