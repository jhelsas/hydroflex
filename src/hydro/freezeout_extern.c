#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "SPHtypes.h"
#include "SPH-lista.h"
#include "SPH-state.h"
#include "SPH-ODEsolver.h"
#include "kernel.h"
#include "utilities.h"

#define EoS EoS_qg

int main(int argc,char **argv){
  int i,l,D,k,err,N_sph,Nspecies,*N;
  double h,kh;
  FILE *dadossph;
  Box *lbox;
  SPHeq_list *sph_eq;
  SPHneq_list **sph_neq;
  
  if(argc==1)
    return (-1);
    
  D=3;
  N_sph=20540;
  h=0.5;
  kh=1.0;
  Nspecies=0;
  
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
    lbox->xmin[l]= -20.0;
    lbox->xmax[l]= 20.0;
    lbox->dx[l] = 0.2;
    
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
  err=sph_neqMalloc(Nspecies,N,D,&sph_neq);
  if(err!=0)
    return (40);
  
  
  if(N_sph!=0){
    dadossph=fopen(argv[1],"r");
    if(dadossph==NULL){
      printf("Falhou em abrir o arquivo de input\n");
      return (30+6);}
  }
  for(i=0;i<N_sph;i+=1){
    err=fscanf(dadossph,"%lf %lf %lf",&(sph_eq[i].p.ni),&(sph_eq[i].p.S),&(sph_eq[i].p.q));
    sph_eq[i].p.q=1.0;
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
  }
  if(dadossph!=NULL){fclose(dadossph);dadossph=NULL;}
  
  printf("Fazendo calculo de setup_sph\n");
  
  err=setup_sph(D,h,kh,N_sph,sph_eq,Nspecies,N,sph_neq,lbox,w_bspline,EoS);
  if(err!=0){
    printf("problemas no setup_sph err=%d\n",err);
    scanf("%d",&err);
  }
  
  printf("Fazendo o calculo do freezeout a tempo constante\n");
  if(argc==2)
    err=freezeout_const_time(D,N_sph,sph_eq,20,10,90,0,"","");
  else if(argc==3)
    err=freezeout_const_time(D,N_sph,sph_eq,20,10,90,1,argv[2],"");
  else
    err=freezeout_const_time(D,N_sph,sph_eq,20,10,90,2,argv[2],argv[3]);
    
  if(err!=0)
    printf("D nao e 3\n");

  return 0;
}
