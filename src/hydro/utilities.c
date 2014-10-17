#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "SPHtypes.h"
#include "SPH-lista.h"
#include "SPH-state.h"
#include "SPH-ODEsolver.h"
#include "kernel.h"

int sph_eqMalloc(int N_sph,int D,SPHeq_list **sph_eq){

  int i,l;
    
  if(sph_eq==NULL)
    return -1;
    
  (*sph_eq)=(SPHeq_list*)malloc(N_sph*sizeof(SPHeq_list));
  if((*sph_eq)==NULL)
    return 1;
  
  for(i=0;i<N_sph;i+=1){
    (*sph_eq)[i].p.id=i;
    (*sph_eq)[i].p.ni=0.1;
    (*sph_eq)[i].p.S=0.0;
    (*sph_eq)[i].p.x=(double*)malloc((1+D)*sizeof(double));
    if((*sph_eq)[i].p.x==NULL)
      return 3;
    (*sph_eq)[i].p.u=(double*)malloc((1+D)*sizeof(double));
    if((*sph_eq)[i].p.u==NULL)
      return 4;
    (*sph_eq)[i].p.v=(double*)malloc((1+D)*sizeof(double));
    if((*sph_eq)[i].p.v==NULL)
      return 5;
    (*sph_eq)[i].p.xa=(double*)malloc((1+D)*sizeof(double));
    if((*sph_eq)[i].p.xa==NULL)
      return 6;
    (*sph_eq)[i].p.ua=(double*)malloc((1+D)*sizeof(double));
    if((*sph_eq)[i].p.ua==NULL)
      return 7;
    (*sph_eq)[i].p.dudt=(double*)malloc((1+D)*sizeof(double));
    if((*sph_eq)[i].p.dudt==NULL)
      return 8;
    (*sph_eq)[i].p.fo=0; 
    (*sph_eq)[i].p.x[0]=0.0;
    (*sph_eq)[i].p.u[0]=0.0;
    (*sph_eq)[i].p.v[0]=1.0;
    (*sph_eq)[i].p.xa[0]=0.0;
    (*sph_eq)[i].p.ua[0]=0.0;
    (*sph_eq)[i].p.dudt[0]=0.0;
    (*sph_eq)[i].p.Ta=0.0;
    (*sph_eq)[i].p.Sa=0.0;
    (*sph_eq)[i].p.rho_pa=1.0;
    for(l=1;l<=D;l+=1){
      (*sph_eq)[i].p.x[l]=0.0;
      (*sph_eq)[i].p.u[l]=0.0;
      (*sph_eq)[i].p.v[l]=0.0;
      (*sph_eq)[i].p.xa[l]=0.0;
      (*sph_eq)[i].p.ua[l]=0.0;
      (*sph_eq)[i].p.dudt[l]=0.0;
    }
  }
  return 0;
}

int sph_neqMalloc(int Nspecies,int *N,int D,SPHneq_list ***sph_neq){
  int k,i,l;
    
  if(sph_neq==NULL)
    return -1;
    
  (*sph_neq)=(SPHneq_list**)malloc(Nspecies*sizeof(SPHneq_list*));
  if((*sph_neq)==NULL)
    return 1;
  for(k=0;k<Nspecies;k+=1){
    (*sph_neq)[k]=(SPHneq_list*)malloc(N[k]*sizeof(SPHneq_list));
    if((*sph_neq)[k]==NULL)
      return 2;
  }
  
  for(k=0;k<Nspecies;k+=1){
    for(i=0;i<N[k];i+=1){
      (*sph_neq)[k][i].p.m=0.1396;
      (*sph_neq)[k][i].p.id=i;
      (*sph_neq)[k][i].p.ni=1.0/10.0;
      (*sph_neq)[k][i].p.N_p=0.0;
      (*sph_neq)[k][i].p.x=(double*)malloc((1+D)*sizeof(double));
      if((*sph_neq)[k][i].p.x==NULL)
        return 3;
      (*sph_neq)[k][i].p.u=(double*)malloc((1+D)*sizeof(double));
      if((*sph_neq)[k][i].p.u==NULL)
        return 4;
      (*sph_neq)[k][i].p.v=(double*)malloc((1+D)*sizeof(double));
      if((*sph_neq)[k][i].p.v==NULL)
        return 5;
      (*sph_neq)[k][i].p.x[0]=0.0;
      (*sph_neq)[k][i].p.u[0]=0.0;
      (*sph_neq)[k][i].p.v[0]=1.0;
      for(l=1;l<=D;l+=1){
        (*sph_neq)[k][i].p.x[l]=0.0;
        (*sph_neq)[k][i].p.u[l]=0.0;
        (*sph_neq)[k][i].p.v[l]=0.0;      
      }
    }
  }
  return 0;
}

int sph_eqFree(int N_sph,SPHeq_list **sph_eq){
  int i;
  
  for(i=0;i<N_sph;i+=1){
    free((*sph_eq)[i].p.x);
    free((*sph_eq)[i].p.u);
    free((*sph_eq)[i].p.v);
  }
  free((*sph_eq));(*sph_eq)=NULL;
  return 0;
}

int sph_neqFree(int Nspecies,int *N,SPHneq_list ***sph_neq){
  int k,i;

  for(k=0;k<Nspecies;k+=1){
    for(i=0;i<N[k];i+=1){
      free((*sph_neq)[k][i].p.x);
      free((*sph_neq)[k][i].p.u);
      free((*sph_neq)[k][i].p.v);
    }
    free((*sph_neq)[k]);
  }
  free((*sph_neq));(*sph_neq)=NULL;
  return 0;
}

int init_simulation(char *config_file,int *D_out,
                    int *N_sph_out,int *Nspecies_out,int **N_out,
                    double *h_out,double *kh_out,
                    double *t0_out,double *tf_out,double *dt_out,
                    SPHeq_list **sph_eq_out,SPHeq_list **sph_eqTemp_out,
                    SPHeq_list **f0_eq_out,SPHeq_list **f1_eq_out,
                    SPHeq_list **f2_eq_out,SPHeq_list **f3_eq_out,
                    SPHneq_list ***sph_neq_out,SPHneq_list ***sph_neqTemp_out,
                    SPHneq_list ***f0_neq_out,SPHneq_list ***f1_neq_out,
                    SPHneq_list ***f2_neq_out,SPHneq_list ***f3_neq_out,
                    Box **lbox_out)
{
  int i,k,l,D,N_sph,Nspecies,*N,err,inD,inNsph;
  double h,kh,t0,tf,dt;
  SPHeq_list *sph_eq,*sph_eqTemp, *f0_eq, *f1_eq, *f2_eq, *f3_eq;
  SPHneq_list **sph_neq,**sph_neqTemp,**f0_neq,**f1_neq,**f2_neq,**f3_neq;
  Box *lbox;
  char sphfile[100+1],checkfile[100+1];
  FILE *dadosin,*dadossph,*dadoscheck;
  
  dadosin=NULL;
  dadossph=NULL;
  dadoscheck=NULL;
  
  printf("Inicializando as Variaveis de Hidrodinamica\n");
  
  dadosin=fopen(config_file,"r");
  if(dadosin==NULL)
    return 1;
    
  err=fscanf(dadosin,"%d %d %d %lf %lf",&D,&N_sph,&Nspecies,&h,&kh);
  if(err!=5)
    return (10);
        
  N=(int*)malloc(Nspecies*sizeof(int));
  if(N==NULL)
    return (10+1);
  for(k=0;k<Nspecies;k+=1){
    err=fscanf(dadosin,"%d",&(N[k]));
    if(err!=1)
      return (10+2);
  }
  /*Teste
  N_sph=1000;
  N[0]=1000;*/
  /*Teste*/
  err=fscanf(dadosin,"%lf %lf %lf",&t0,&tf,&dt);
  if(err!=3)
    return (10+3);
  
  lbox=(Box*)malloc(sizeof(Box));
  if(lbox==NULL)
    return (20);
  lbox->Nspecies=Nspecies;
  lbox->D=D;
  lbox->Nbox=1;
  
  lbox->xmin=(double*)malloc((lbox->D)*sizeof(double));
  if(lbox->xmin==NULL)
    return (20+1);
  lbox->xmax=(double*)malloc((lbox->D)*sizeof(double));
  if(lbox->xmax==NULL)
    return (20+2);
  lbox->dx  =(double*)malloc((lbox->D)*sizeof(double));
  if(lbox->dx==NULL)
    return (20+3);
  lbox->Lbox=(int*)malloc((lbox->D)*sizeof(int));
  if(lbox->Lbox==NULL)
    return (20+4);
  lbox->k=(int*)malloc((lbox->D)*sizeof(int));
  if(lbox->k==NULL)
    return (20+5);
  lbox->ind=(int*)malloc((lbox->D)*sizeof(int));
  if(lbox->ind==NULL)
    return (20+6);
          
  for(l=0;l<D;l+=1){
    err=fscanf(dadosin,"%lf %lf %lf",&(lbox->xmin[l]),&(lbox->xmax[l]),&(lbox->dx[l]));
    if(err!=3)
      return (20+7);
    lbox->Lbox[l]=(int)((lbox->xmax[l]-lbox->xmin[l])/lbox->dx[l]);
    lbox->Nbox *= lbox->Lbox[l];  
  }
  lbox->box=(SPHeq_list**)malloc((lbox->Nbox)*sizeof(SPHeq_list*));
  if(lbox->box==NULL)
    return (20+8);
  
  lbox->neqbox=(SPHneq_list***)malloc((lbox->Nspecies)*sizeof(SPHneq_list**));
  if(lbox->neqbox==NULL)
    return (20+9);
    
  for(k=0;k<Nspecies;k+=1){
    lbox->neqbox[k]=(SPHneq_list**)malloc((lbox->Nbox)*sizeof(SPHneq_list**));
    if(lbox->neqbox[k]==NULL)
      return (20+9);
  }
    
  err=sph_eqMalloc(N_sph,D,&sph_eq);
  if(err!=0)
    return (30);
  err=sph_eqMalloc(N_sph,D,&sph_eqTemp);
  if(err!=0)
    return (30+1);
  err=sph_eqMalloc(N_sph,D,&f0_eq);
  if(err!=0)
    return (30+2); 
  err=sph_eqMalloc(N_sph,D,&f1_eq);
  if(err!=0)
    return (30+3); 
  err=sph_eqMalloc(N_sph,D,&f2_eq);
  if(err!=0)
    return (30+4);  
  err=sph_eqMalloc(N_sph,D,&f3_eq);
  if(err!=0)
    return (30+5);
  
  if(N_sph!=0){
    if(fscanf(dadosin,"%s",sphfile)!=1)
      printf("Problemas leitura sphfile\n");
    dadossph=fopen(sphfile,"r");
    if(dadossph==NULL)
      return (30+6);
  }
  /*
  err = fscanf(dadossph,"%lf %lf",&(inD),&(inNsph));
  if(err!=2)
    return (30+6);
   */
  for(i=0;i<N_sph;i+=1){
    err=fscanf(dadossph,"%lf %lf %lf",&(sph_eq[i].p.ni),&(sph_eq[i].p.q),&(sph_eq[i].p.S));
    if(err!=3)
      return (30+7);
    for(l=1;l<=D;l+=1){
      err=fscanf(dadossph,"%lf",&(sph_eq[i].p.x[l]));
      if(err!=1)
        return (30+8);
    }
    for(l=1;l<=D;l+=1){
      err=fscanf(dadossph,"%lf",&(sph_eq[i].p.u[l]));
      if(err!=1)
        return (30+9);
    }
    
    sph_eqTemp[i].p.ni=sph_eq[i].p.ni;
    sph_eqTemp[i].p.q=sph_eq[i].p.q;
  }
  if(dadossph!=NULL){fclose(dadossph);dadossph=NULL;}
  
  printf("Nspecies=%d D=%d N=%p\n",Nspecies,D,N);
  
  err=sph_neqMalloc(Nspecies,N,D,&sph_neq);
  if(err!=0)
    return (40);
  err=sph_neqMalloc(Nspecies,N,D,&sph_neqTemp);
  if(err!=0)
    return (40+1);
  err=sph_neqMalloc(Nspecies,N,D,&f0_neq);
  if(err!=0)
    return (40+2);
  err=sph_neqMalloc(Nspecies,N,D,&f1_neq);
  if(err!=0)
    return (40+3);
  err=sph_neqMalloc(Nspecies,N,D,&f2_neq);
  if(err!=0)
    return (40+4);
  err=sph_neqMalloc(Nspecies,N,D,&f3_neq);
  if(err!=0)
    return (40+5);
  
  printf("OI-2\n");
  
  for(k=0;k<Nspecies;k+=1){
    if(fscanf(dadosin,"%s",sphfile)!=1)
      printf("Problemas leitura sphfile-%d\n",k);
    dadossph=fopen(sphfile,"r");
    if(dadossph==NULL)
      return (40+6);
      
    for(i=0;i<N[k];i+=1){
      err=fscanf(dadossph,"%lf %lf %lf %lf",&(sph_neq[k][i].p.m),&(sph_neq[k][i].p.q),&(sph_neq[k][i].p.ni),&(sph_neq[k][i].p.N_p));
      sph_neqTemp[k][i].p.m=sph_neq[k][i].p.m;
      sph_neqTemp[k][i].p.q=sph_neq[k][i].p.q;
      sph_neqTemp[k][i].p.ni=sph_neq[k][i].p.ni;
      if(err!=4)
        return (40+7);
      for(l=0;l<=D;l+=1){
        err=fscanf(dadossph,"%lf",&(sph_neq[k][i].p.x[l]));
        if(err!=1)
          return (40+8);
      }
      for(l=0;l<=D;l+=1){
        err=fscanf(dadossph,"%lf",&(sph_neq[k][i].p.u[l]));
        if(err!=1)
          return (40+9);
      }
    }
    fclose(dadossph);
  }
  if(dadosin!=NULL){fclose(dadosin);dadosin=NULL;}

  (*D_out)=D;
  (*N_sph_out)=N_sph;
  (*Nspecies_out)=Nspecies;
  (*N_out)=N;
  (*t0_out)=t0;
  (*tf_out)=tf;
  (*dt_out)=dt;
  (*h_out)=h;
  (*kh_out)=kh;
  
  (*sph_eq_out)=sph_eq;
  (*sph_eqTemp_out)=sph_eqTemp;
  (*f0_eq_out)=f0_eq;
  (*f1_eq_out)=f1_eq;
  (*f2_eq_out)=f2_eq;
  (*f3_eq_out)=f3_eq;
  
  (*sph_neq_out)=sph_neq;
  (*sph_neqTemp_out)=sph_neqTemp;
  (*f0_neq_out)=f0_neq;
  (*f1_neq_out)=f1_neq;
  (*f2_neq_out)=f2_neq;
  (*f3_neq_out)=f3_neq;
  
  (*lbox_out)=lbox;
  
  return 0;
}

int init_printing_grid(int Npoints,int D,char *filename,double ***xp_out){
  int l,p;
  double **xp;
  FILE *dadosin;
  
  xp=(double**)malloc(Npoints*sizeof(double*));
  if(xp==NULL)
    return 1;
  for(p=0;p<Npoints;p+=1){
    xp[p]=(double*)malloc((D+1)*sizeof(double));
    if(xp[p]==NULL)
      return p;
  }
  
  dadosin=fopen(filename,"r");
  if(dadosin==NULL)
    return (-1);
  for(p=0;p<Npoints;p+=1){
    for(l=1;l<=D;l+=1)
      if(fscanf(dadosin,"%lf",&(xp[p][l]))!=1)
        return (-2);
  }
  fclose(dadosin);
  
  *xp_out=xp;
  return 0;
}


int create_grid(int D,double **xpo,double *xl,double *xu,
                double *dx,int *Np){
  double *xp;
  if(D<=0)
    return -1;
  else if(D==1){
    int i,Nx;
    Nx=(int)((xu[1]-xl[1])/dx[1]);
    xp=(double *)malloc(Nx*(D+1)*sizeof(double));
    for(i=0;i<Nx;i+=1){
      xp[i*(D+1)+0] = 0.0;
      xp[i*(D+1)+1] = xl[1]+((double)i)*(dx[1]);
    }
    *Np=Nx;
  }
  else if(D==2){
    int k,i,j,Nx,Ny;
    Nx=(int)((xu[1]-xl[1])/dx[1]);
    Ny=(int)((xu[2]-xl[2])/dx[2]);
    xp=(double*)malloc(Nx*Ny*(D+1)*sizeof(double));
        
    for(j=0;j<Ny;j+=1)
      for(i=0;i<Nx;i+=1){
        xp[(j*Nx+i)*(D+1)+0] = 0.0;
        xp[(j*Nx+i)*(D+1)+1] = xl[1]+((double)i)*(dx[1]);
        xp[(j*Nx+i)*(D+1)+2] = xl[2]+((double)j)*(dx[2]);
      }
    
    *Np=Nx*Ny;
  }
  else if(D==3){
    int i,j,k,Nx,Ny,Nz;
    Nx=(int)((xu[1]-xl[1])/dx[1]);
    Ny=(int)((xu[2]-xl[2])/dx[2]);
    Nz=(int)((xu[3]-xl[3])/dx[3]);
    xp=(double*)malloc(Nx*Ny*Nz*(D+1)*sizeof(double));
    for(k=0;k<Nz;k+=1)
      for(j=0;j<Ny;j+=1)
        for(i=0;i<Nx;i+=1){
          xp[((k*Ny+j)*Nx+i)*(D+1)+1] = 0.0;
          xp[((k*Ny+j)*Nx+i)*(D+1)+1] = xl[1]+((double)i)*(dx[1]);
          xp[((k*Ny+j)*Nx+i)*(D+1)+2] = xl[2]+((double)j)*(dx[2]);
          xp[((k*Ny+j)*Nx+i)*(D+1)+3] = xl[3]+((double)k)*(dx[3]);
        }
      
    *Np=Nx*Ny*Nz;
  }
  *xpo=xp;
  return 0;
}

int create_X_line(int D,double **xpo,double *xl,double *xu,
                  double *dx,int *Np)
{
  double *xp;
  if(D<=0)
    return -1;
  else if(D==1){
    int i,Nx;
    Nx=(int)((xu[1]-xl[1])/dx[1]);
    xp=(double *)malloc(Nx*(D+1)*sizeof(double));
    for(i=0;i<Nx;i+=1){
      xp[i*(D+1)+0] = 0.0;
      xp[i*(D+1)+1] = xl[1]+((double)i)*(dx[1]);
    }
    *Np=Nx;
  }
  else if(D==2){
    int k,i,j,Nx,Ny;
    Nx=(int)((xu[1]-xl[1])/dx[1]);
    xp=(double*)malloc(Nx*(D+1)*sizeof(double));
        
    for(i=0;i<Nx;i+=1){
      xp[i*(D+1)+0] = 0.0;
      xp[i*(D+1)+1] = xl[1]+((double)i)*(dx[1]);
      xp[i*(D+1)+2] = 0.0;
    }
    
    *Np=Nx;
  }
  else if(D==3){
    int i,j,k,Nx,Ny,Nz;
    Nx=(int)((xu[1]-xl[1])/dx[1]);
    xp=(double*)malloc(Nx*(D+1)*sizeof(double));
    for(i=0;i<Nx;i+=1){
      xp[i*(D+1)+0] = 0.0;
      xp[i*(D+1)+1] = xl[1]+((double)i)*(dx[1]);
      xp[i*(D+1)+2] = 0.0;
      xp[i*(D+1)+3] = 0.0;
    }
      
    *Np=Nx;
  }
  *xpo=xp;
  return 0;
}

int create_Y_line(int D,double **xpo,double *xl,double *xu,
                  double *dx,int *Np)
{
  double *xp;
  if(D<=0)
    return -1;
  else if(D==1){
    return -2;
  }
  else if(D==2){
    int k,i,j,Nx,Ny;
    Ny=(int)((xu[2]-xl[2])/dx[2]);
    xp=(double*)malloc(Ny*(D+1)*sizeof(double));
        
    for(i=0;i<Ny;i+=1){
      xp[i*(D+1)+0] = 0.0;
      xp[i*(D+1)+1] = 0.0;
      xp[i*(D+1)+2] = xl[2]+((double)i)*(dx[2]);
    }
    
    *Np=Ny;
  }
  else if(D==3){
    int i,j,k,Nx,Ny,Nz;
    Ny=(int)((xu[2]-xl[2])/dx[2]);
    xp=(double*)malloc(Ny*(D+1)*sizeof(double));
    for(i=0;i<Ny;i+=1){
      xp[i*(D+1)+0] = 0.0;
      xp[i*(D+1)+1] = 0.0;
      xp[i*(D+1)+2] = xl[2]+((double)i)*(dx[2]);
      xp[i*(D+1)+3] = 0.0;
    }
      
    *Np=Ny;
  }
  *xpo=xp;
  return 0;
}

int print_proper_entropy(int D,double h,double kh,Box *lbox,double (*w)(int,double,double),
                         char *filename,int Npoints,double **xp)
{
  int p,l;
  double x[D+1],s,dist;
  FILE *dadosout;
  SPHeq_list *inter_eq,*aux_eq;
  
  dadosout=fopen(filename,"w");
  for(p=0;p<Npoints;p+=1){
    for(l=1;l<=D;l+=1)
      x[l]=xp[p][l];
    if(fetch_inter_eq(x,lbox,kh,&inter_eq)!=0)
      return (10+1);
    s=0.0;
    while(inter_eq!=NULL){
      dist=0.0;
      for(l=1;l<=D;l+=1)
        dist+=(x[l]-(inter_eq->p).x[l])*(x[l]-(inter_eq->p).x[l]);
      dist=sqrt(dist);
          
      s+= ((inter_eq->p).ni)*(((inter_eq->p).S)/((inter_eq->p).u[0]))*w(D,dist,h);

      aux_eq=inter_eq;
      inter_eq=inter_eq->inext;
      aux_eq->inext=NULL;
    }
    for(l=1;l<=D;l+=1)
      fprintf(dadosout,"%lf ",x[l]);
    fprintf(dadosout,"%lf\n",s);
  }
  fclose(dadosout);
  return 0;
}

int print_proper_energy(int D,int Nspecies,double h,double kh,Box *lbox,
                        double (*w)(int,double,double),char *filename,
                        int Npoints,double **xp,double *displ)
{
  int k,p,l;
  double L,omega;
  double x[D+1],e_p,T,dist,r,rt,lambda;
  FILE *dadosout;
  SPHeq_list *inter_eq,*aux_eq;
  SPHneq_list **inter_neq,*aux_neq;
  
  omega=0.1;
  L=1.0;
  
  inter_neq=(SPHneq_list**)malloc(Nspecies*sizeof(SPHneq_list*));
  
  dadosout=fopen(filename,"w");
  for(p=0;p<Npoints;p+=1){
    x[0]=xp[p][0];
    for(l=1;l<=D;l+=1)
      x[l]=xp[p][l]+displ[l];
    e_p=0.0;
    T=0.0;
    if(fetch_inter_eq(x,lbox,kh,&inter_eq)!=0)
      return (10+1);
    while(inter_eq!=NULL){
      dist=0.0;
      for(l=1;l<=D;l+=1)
        dist+=(x[l]-(inter_eq->p).x[l])*(x[l]-(inter_eq->p).x[l]);
      dist=sqrt(dist);
          
      e_p+= ((inter_eq->p).ni)*(((inter_eq->p).e_p)/((inter_eq->p).rho))*w(D,dist,h);
      T += ((inter_eq->p).ni)*(((inter_eq->p).T)/((inter_eq->p).rho))*w(D,dist,h);

      aux_eq=inter_eq;
      inter_eq=inter_eq->inext;
      aux_eq->inext=NULL;
    }
    
    for(k=0;k<Nspecies;k+=1){
      if(fetch_inter_neq(x,lbox,kh,k,&inter_neq)!=0)
        return (10+1);
      while(inter_neq[k]!=NULL){
        dist=0.0;
        for(l=1;l<=D;l+=1)
          dist+=(x[l]-(inter_neq[k]->p).x[l])*(x[l]-(inter_neq[k]->p).x[l]);
        dist=sqrt(dist);
      
        e_p += ((inter_neq[k]->p).ni)*((inter_neq[k]->p).m)*((inter_neq[k]->p).N_p)*w(D,dist,h);
    
        aux_neq=inter_neq[k];
        inter_neq[k]=inter_neq[k]->inext;
        aux_neq->inext=NULL;
      }
    }
    
    r=0.0;
    for(l=1;l<=D;l+=1)
      r+=(x[l])*(x[l]);
    r=sqrt(r);
        
    rt=0.0;
    for(l=1;l<D;l+=1)
      rt+=(x[l])*(x[l]);
    rt=sqrt(rt);
    
    lambda = (L*L + (x[0]-r)*(x[0]-r))*(L*L+(x[0]+r)*(x[0]+r))-4.*omega*omega*rt*rt*L*L;
    
    for(l=1;l<=D;l+=1)
      fprintf(dadosout,"%lf ",x[l]);
    fprintf(dadosout,"%lf %lf %lf\n",e_p/(1340.713891*2.153619),T,1./(lambda*lambda));
  }
  fclose(dadosout);
  free(inter_neq);
  return 0;
}

double gubser_entropy(double *x,size_t dim,double *par){
  int err;
  
  double *p = (double*)(par);
  double s0=p[0], q=p[1],tau=p[2];
  double r2,lambda;
  r2=x[0]*x[0]+x[1]*x[1];
  lambda = 1.+2.*q*q*(tau*tau+r2) + q*q*q*q*(tau*tau-r2)*(tau*tau-r2);
  
  return (s0*(2.*q)*(2.*q)*(1.+q*q*(tau*tau+r2)))/(lambda*sqrt(lambda));
}

double gubser_proper_entropy(double *x,size_t dim,void *par){
  int err;
  
  double *p = (double*)(par);
  double s0=p[0], q=p[1],tau=p[2];
  double r2,lambda;
  r2=x[0]*x[0]+x[1]*x[1];
  lambda = 1.+2.*q*q*(tau*tau+r2) + q*q*q*q*(tau*tau-r2)*(tau*tau-r2);
  
  return (s0*(2.*q)*(2.*q))/(tau*lambda);
  
}

double gubser_proper_energy(double *x,size_t dim, void *par){
  int err;
  
  double *p = (double*)(par);
  double s0=p[0], q=p[1],tau=p[2];
  double r2,lambda,epsilon;
  r2=x[0]*x[0]+x[1]*x[1];
  lambda = 1.+2.*q*q*(tau*tau+r2) + q*q*q*q*(tau*tau-r2)*(tau*tau-r2);
  
  epsilon = 3.0*C_qg*pow(s0,4.0/3.0)*pow(2.*q,8./3.);
  epsilon = epsilon/pow(tau*lambda,4./3.);
  
  return epsilon;
}

int print_new(int D,int Nspecies,double t,double h,double kh,Box *lbox,
              double (*w)(int,double,double),char *filename,
              int Npoints,double *xp,double *displ,int mode)
{
  int k,p,l;
  double L,q[3];
  double x[D+1],r[D],s_p,s,e_p,dist;
  FILE *dadosout;
  SPHeq_list *inter_eq,*aux_eq;
  SPHneq_list **inter_neq,*aux_neq;
  
  inter_neq=(SPHneq_list**)malloc(Nspecies*sizeof(SPHneq_list*));
  
  q[0]=1.0;q[1]=1.0;q[2]=t;
  dadosout=fopen(filename,"w");
  for(p=0;p<Npoints;p+=1){
    x[0]=xp[p*(D+1)+0];
    for(l=1;l<=D;l+=1){
      x[l]=xp[p*(D+1)+l];r[l-1]=x[l];}
    
    if(fetch_inter_eq(x,lbox,kh,&inter_eq)!=0)
      return (10+1);
    s=0.0;
    s_p=0.0;
    e_p=0.0;
    while(inter_eq!=NULL){
      dist=0.0;
      for(l=1;l<=D;l+=1)
        dist+=(xp[p*(D+1)+l]-(inter_eq->p).x[l])*(xp[p*(D+1)+l]-(inter_eq->p).x[l]);
      dist=sqrt(dist);
                
      s += ((inter_eq->p).ni)*((inter_eq->p).S)*w(D,dist,h);
      s_p += ((inter_eq->p).ni)*(((inter_eq->p).S)/((inter_eq->p).u[0]))*w(D,dist,h);
      e_p += ((inter_eq->p).ni)*(((inter_eq->p).e_p)/((inter_eq->p).rho))*w(D,dist,h);
      
      aux_eq=inter_eq;
      inter_eq=inter_eq->inext;
      aux_eq->inext=NULL;
    }
    s_p = s_p/t;
    
    for(l=1;l<=D;l+=1)
      fprintf(dadosout,"%lf ",xp[p*(D+1)+l]);
    
    fprintf(dadosout,"%lf %lf ",s,gubser_entropy(r,D,q));
    fprintf(dadosout,"%lf %lf ",s_p,gubser_proper_entropy(r,D,q));
    fprintf(dadosout,"%lf %lf ",e_p,gubser_proper_energy(r,D,q));
    fprintf(dadosout,"\n");
  }
  fclose(dadosout);
  free(inter_neq);
  return 0;
}

int print_4vel_profile(int D,int Nspecies,double h,double kh,Box *lbox,
                       double (*w)(int,double,double),char *filename,int Npoints,double **xp)
{
  int k,p,l;
  double x[D+1],u[D+1],dist,rho;
  FILE *dadosout;
  SPHeq_list *inter_eq,*aux_eq;
  SPHneq_list **inter_neq,*aux_neq;
  
  inter_neq=(SPHneq_list**)malloc(Nspecies*sizeof(SPHneq_list*));
  
  dadosout=fopen(filename,"w");
  for(p=0;p<Npoints;p+=1){
    rho=0.0;
    for(l=1;l<=D;l+=1){
      x[l]=xp[p][l];
      u[l]=0.0;
    }
    if(fetch_inter_eq(x,lbox,kh,&inter_eq)!=0)
      return (10+1);
    while(inter_eq!=NULL){
      dist=0.0;
      for(l=1;l<=D;l+=1)
        dist+=(x[l]-(inter_eq->p).x[l])*(x[l]-(inter_eq->p).x[l]);
      dist=sqrt(dist);
      /*
      for(l=1;l<=D;l+=1)
        u[l] += ((inter_eq->p).ni)*(((inter_eq->p).u[l])/((inter_eq->p).rho))*w(D,dist,h);
      */
      rho+= ((inter_eq->p).ni)*w(D,dist,h);
      for(l=1;l<=D;l+=1)
        u[l] += ((inter_eq->p).ni)*((inter_eq->p).u[l])*w(D,dist,h);
      
      aux_eq=inter_eq;
      inter_eq=inter_eq->inext;
      aux_eq->inext=NULL;
    }
    for(l=1;l<=D;l+=1)
      fprintf(dadosout,"%lf ",x[l]);
    
    for(l=1;l<=D;l+=1)
      u[l]/=rho;
    for(l=1;l<=D;l+=1)
      fprintf(dadosout,"%lf ",u[l]);
      
    rho=0.0;
    for(l=1;l<=D;l+=1)
      u[l]=0.0;
      
    for(k=0;k<Nspecies;k+=1){
      if(fetch_inter_neq(x,lbox,kh,k,&inter_neq)!=0)
        return (10+1);
      
      while(inter_neq[k]!=NULL){
        dist=0.0;
        for(l=1;l<=D;l+=1)
          dist+=(x[l]-(inter_neq[k]->p).x[l])*(x[l]-(inter_neq[k]->p).x[l]);
        dist=sqrt(dist);
        /*
        for(l=1;l<=D;l+=1)
          u[l] += ((inter_neq[k]->p).ni)*(((inter_neq[k]->p).u[l])/((inter_neq[k]->p).rho))*w(D,dist,h);
        */
        rho+=((inter_neq[k]->p).ni)*w(D,dist,h);
        for(l=1;l<=D;l+=1)
          u[l]+=((inter_neq[k]->p).ni)*((inter_neq[k]->p).u[l])*w(D,dist,h);
        
        aux_neq=inter_neq[k];
        inter_neq[k]=inter_neq[k]->inext;
        aux_neq->inext=NULL;
      }
    }
    for(l=1;l<=D;l+=1)
      u[l]/=rho;
    for(l=1;l<=D;l+=1)
      fprintf(dadosout,"%lf ",u[l]);
    fprintf(dadosout,"\n");
  }
  fclose(dadosout);
  free(inter_neq);
  return 0;
}

int print_charge_density(int D,int Nspecies,double h,double kh,Box *lbox,
                         double (*w)(int,double,double),char *filename,
                         int Npoints,double **xp,double *displ)
{
  int k,p,l;
  double x[D+1],nc,ncp,T,dist;
  FILE *dadosout;
  SPHeq_list *inter_eq,*aux_eq;
  SPHneq_list **inter_neq,*aux_neq;
  
  inter_neq=(SPHneq_list**)malloc(Nspecies*sizeof(SPHneq_list*));
  
  dadosout=fopen(filename,"w");
  for(p=0;p<Npoints;p+=1){
    for(l=1;l<=D;l+=1)
      x[l]=xp[p][l]+displ[l];
    nc=0.0;
    ncp=0.0;
    if(fetch_inter_eq(x,lbox,kh,&inter_eq)!=0)
      return (10+1);
    while(inter_eq!=NULL){
      dist=0.0;
      for(l=1;l<=D;l+=1)
        dist+=(x[l]-(inter_eq->p).x[l])*(x[l]-(inter_eq->p).x[l]);
      dist=sqrt(dist);
          
      nc+= ((inter_eq->p).ni)*((inter_eq->p).Nc)*w(D,dist,h);
      ncp+=((inter_eq->p).ni)*(((inter_eq->p).Nc)/((inter_eq->p).u[0]))*w(D,dist,h);
      
      aux_eq=inter_eq;
      inter_eq=inter_eq->inext;
      aux_eq->inext=NULL;
    }
    
    for(l=1;l<=D;l+=1)
      fprintf(dadosout,"%lf ",x[l]);
    fprintf(dadosout,"%lf %lf\n",nc,ncp);
  }
  fclose(dadosout);
  free(inter_neq);
  return 0;
}

int scripting(int argc,char **argv,
              double t0,double tf,double dt){
  
  int count,secf=0,namef=0,gfcf=0,titlef=0;
  double graph_size,t;
  char filename[100+1];
  FILE *script_t,*script_r;
  
  if(argc==1){printf("problemas\n");return 1;}
  
  printf("hashname=%s\n",argv[1]);
  if(argc>=3)
    secf=atoi(argv[2]);
  printf("section format=%d -- 0=wxt | 1==qt\n",secf);
  
  if(argc>=4)
    namef=atoi(argv[3]);
  printf("name format=%d -- 0=adaptive | 1==fixed\n",namef);
    
  if(argc>=5)
    gfcf=atoi(argv[4]);  
  printf("graphics format=%d -- 0=eps/pdf tex | 1==jpg\n",gfcf);
  
  if(argc>=6)
    titlef=atoi(argv[5]);
  printf("title format=%d -- 0= no title | 1== just time | 2==full\n",titlef);

  if(argc>=7)
    printf("gfc/localization = %s\n",argv[6]);

  if(argc>1){
	  graph_size=10;  
	  
    if(argc>=7)
	    sprintf(filename,"gfc/%s/e_trans-%s.gnp",argv[6],argv[1]);
    else
      sprintf(filename,"gfc/e_trans-%s.gnp",argv[1]);
    script_t=fopen(filename,"w");
    fprintf(script_t,"set dgrid3d %d,%d\n",161,161);
    fprintf(script_t,"set pm3d\n");
    fprintf(script_t,"set hidden3d\n");
    fprintf(script_t,"set view map\n");
    fprintf(script_t,"unset key\n");
    fprintf(script_t,"set size square\n");
    fprintf(script_t,"set xl 'x (fm)'\n");
    fprintf(script_t,"set yl 'y (fm)'\n");
    fprintf(script_t,"set xr [-%.1f:%.1f]\n",graph_size,graph_size);
    fprintf(script_t,"set yr [-%.1f:%.1f]\n",graph_size,graph_size);
    fprintf(script_t,"set border 1+2+4\n");
    fprintf(script_t,"set lmargin at screen 0.0\n");
    fprintf(script_t,"set rmargin at screen 1.0\n");
    fprintf(script_t,"set bmargin at screen 0.2\n");
    fprintf(script_t,"set tmargin at screen 0.9\n");
    fflush(script_t);
    
    if(argc>=7)
	    sprintf(filename,"gfc/%s/e_react-%s.gnp",argv[6],argv[1]);
    else
	    sprintf(filename,"gfc/e_react-%s.gnp",argv[1]);
    script_r=fopen(filename,"w");
    fprintf(script_r,"set dgrid3d %d,%d\n",161,161);
    fprintf(script_r,"set pm3d\n");
    fprintf(script_r,"set hidden3d\n");
    fprintf(script_r,"set view map\n");
    fprintf(script_r,"unset key\n");
    fprintf(script_r,"set size square\n");
    fprintf(script_r,"set xl 'x (fm)'\n");
    fprintf(script_r,"set yl 'z (fm)'\n");
    fprintf(script_r,"set xr [-%.1f:%.1f]\n",graph_size,graph_size);
    fprintf(script_r,"set yr [-%.1f:%.1f]\n",graph_size,graph_size);
    fprintf(script_r,"set border 1+2+4\n");
    fprintf(script_r,"set lmargin at screen 0.0\n");
    fprintf(script_r,"set rmargin at screen 1.0\n");
    fprintf(script_r,"set bmargin at screen 0.2\n");
    fprintf(script_r,"set tmargin at screen 0.9\n");
    fflush(script_t);
  }
  
  count=0;
  for(t=t0;t<=tf;t+=dt){
    
    if( count == 0 || (((count)%5)==0) || count==20 ){
      
      if(argc>1){
        fprintf(script_r,"set cbl '$\\varepsilon$ ($GeV/fm^3$)'\n");
        if(titlef==0)
          fprintf(script_r,"unset title\n");
        else if(titlef==1)
          fprintf(script_r,"set title 't=%.1lf'\n",t);
        else
          fprintf(script_r,"set title 'Energy Density (GeV/fm^3) - Reaction Plane - (t=%.1lf)'\n",t);
          
        if(namef==0)
          fprintf(script_r,"splot 'e_react-%s-(t=%.1lf).dat' using 1:3:4\n",argv[1],t);
        else
          fprintf(script_r,"splot 'energy-reactplane-(t=%lf).dat'using 1:3:4\n",t);
        
        if(gfcf==0){
          fprintf(script_r,"set term epslatex standalone color colortext 12\n");
          fprintf(script_r,"set out 'e_react-%s-(t=%.1lf).tex'\n",argv[1],t);}
        else{
          fprintf(script_r,"set term jpeg\n");
          fprintf(script_r,"set out 'e_react-%s-(t=%.1lf).jpg'\n",argv[1],t);}
                
        fprintf(script_r,"replot\n");
        if(secf==0)
          fprintf(script_r,"set term wxt\n");
        else
          fprintf(script_r,"set term qt\n");
        
        
        fprintf(script_r,"set cbl '$T$ ($GeV$)'\n");
        if(titlef==0)
          fprintf(script_r,"unset title\n");
        else if(titlef==1)
          fprintf(script_r,"set title 't=%.1lf'\n",t);
        else
          fprintf(script_r,"set title 'Temperature (GeV) - Reaction Plane - (t=%.1lf)'\n",t);
        
        if(namef==0)
          fprintf(script_r,"splot 'e_react-%s-(t=%.1lf).dat' using 1:3:5\n",argv[1],t);
        else
          fprintf(script_r,"splot 'energy-reactplane-(t=%lf).dat'using 1:3:5\n",t);
        
        if(gfcf==0){
          fprintf(script_r,"set term epslatex standalone color colortext 12\n");
          fprintf(script_r,"set out 'T_react-%s-(t=%.1lf).tex'\n",argv[1],t);}
        else{
          fprintf(script_r,"set term jpeg\n");
          fprintf(script_r,"set out 'T_react-%s-(t=%.1lf).jpg'\n",argv[1],t);}
          
        fprintf(script_r,"replot\n");
        if(secf==0)
          fprintf(script_r,"set term wxt\n");
        else
          fprintf(script_r,"set term qt\n");
          
        fflush(script_r);
      }
            
      if(argc>1){
        fprintf(script_t,"set cbl '$\\varepsilon$ ($GeV/fm^3$)'\n");
        if(titlef==0)
          fprintf(script_t,"unset title\n");
        else if(titlef==1)
          fprintf(script_t,"set title 't=%.1lf'\n",t);
        else
          fprintf(script_t,"set title 'Energy Density (GeV/fm^3) - Transverse Plane - (t=%.1lf)'\n",t);
        
        if(namef==0)
          fprintf(script_t,"splot 'e_trans-%s-(t=%.1lf).dat' using 1:2:4\n",argv[1],t);
        else
          fprintf(script_t,"splot 'energy-trans-(t=%lf).dat' using 1:2:4\n",t);
        
        if(gfcf==0){
          fprintf(script_t,"set term epslatex standalone color colortext 12\n");
          fprintf(script_t,"set out 'e_trans-%s-(t=%.1lf).tex'\n",argv[1],t);}
        else{
          fprintf(script_t,"set term jpeg\n");
          fprintf(script_t,"set out 'e_trans-%s-(t=%.1lf).jpg'\n",argv[1],t);}
        
        fprintf(script_t,"replot\n");
        if(secf==0)
          fprintf(script_t,"set term wxt\n");
        else
          fprintf(script_t,"set term qt\n");
        
        fprintf(script_t,"set cbl '$T$ ($GeV$)'\n");
        if(titlef==0)
          fprintf(script_t,"unset title\n");
        else if(titlef==1)
          fprintf(script_t,"set title 't=%.1lf'\n",t);
        else
          fprintf(script_t,"set title 'Temperature (GeV) - Transverse Plane - (t=%.1lf)'\n",t);
        
        if(namef==0)
          fprintf(script_t,"splot 'e_trans-%s-(t=%.1lf).dat' using 1:2:5\n",argv[1],t);
        else
          fprintf(script_t,"splot 'energy-trans-(t=%lf).dat' using 1:2:5\n",t);
        
        if(gfcf==0){
          fprintf(script_t,"set term epslatex standalone color colortext 12\n");
          fprintf(script_t,"set out 'T_trans-%s-(t=%.1lf).tex'\n",argv[1],t);}
        else{
          fprintf(script_t,"set term jpeg\n");
          fprintf(script_t,"set out 'T_trans-%s-(t=%.1lf).jpg'\n",argv[1],t);}
        
        fprintf(script_t,"replot\n");
        if(secf==0)
          fprintf(script_t,"set term wxt\n");
        else
          fprintf(script_t,"set term qt\n");
          
        fflush(script_t);
      }
        
      printf("script - t=%.12f\n",t);
    }
    
    
    count+=1;
  }
  
  fclose(script_t);
  fclose(script_r);
  
  return 0;
}
