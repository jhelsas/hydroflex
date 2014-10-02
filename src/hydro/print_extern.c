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
  int Npoints,Nrp,p;
  int i,l,D,k,err,N_sph,Nspecies,*N;
  double h,kh;
  double **xp,**xv,**xrp,**xvrp;
  double *x, *displ;
  char filename[140+1];
  FILE *dadossph;
  Box *lbox;
  funDDtD w,Dw;
  SPHeq_list *sph_eq;
  SPHneq_list **sph_neq;
  
  if(argc==1)
    return (-1);
    
  D=3;
  N_sph=atoi(argv[1]);
  h=0.5;
  kh=1.0;
  Nspecies=0;
  
  w.f=w_bspline;
  Dw.f=Dw_bspline;
  
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
    dadossph=fopen(argv[2],"r");
    if(dadossph==NULL){
      printf("Falhou em abrir o arquivo de input\n");
      return (30+6);}
  }
  for(i=0;i<N_sph;i+=1){
	/*
    err=fscanf(dadossph,"%lf %lf %lf",&(sph_eq[i].p.ni),&(sph_eq[i].p.S),&(sph_eq[i].p.Nc));
    if(sph_eq[i].p.Nc>0.0)
      sph_eq[i].p.q=1.0;
    else if(sph_eq[i].p.Nc<0.0)
      sph_eq[i].p.q=-1.0;
    else
      sph_eq[i].p.q=0.0;
    */
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
  }
  if(dadossph!=NULL){fclose(dadossph);dadossph=NULL;}
  
  
  x=(double*)malloc((D+1)*sizeof(double));
  if(x==NULL)
    printf("problemas alocacao x\n");
  
  displ=(double*)malloc((D+1)*sizeof(double));
  if(displ==NULL)
    printf("problemas alocacao displ\n");
  
  /*  
  for(l=0;l<=D;l+=1)
    displ[l]=0.0;
  */
  
  sprintf(filename,"gfc/points-trans.gfcp");
  Npoints=161*161;
  if(init_printing_grid(Npoints,D,filename,&xp)!=0)
    printf("Erro na Inicializacao de xp\n");
      
  sprintf(filename,"gfc/points-reactplane.gfcp");
  Nrp=161*161;
  if(init_printing_grid(Npoints,D,filename,&xrp)!=0)
    printf("Erro na Inicializacao de xrp\n");
  /*
  sprintf(filename,"gfc/points-4vel.gfcp");
  Nvel=81*81;
  if(init_printing_grid(Nvel,D,filename,&xv)!=0)
    printf("Erro na Inicializacao de xv\n");
    
  sprintf(filename,"gfc/points-4velrp.gfcp");
  Nvel=81*81;
  if(init_printing_grid(Nvel,D,filename,&xvrp)!=0)
    printf("Erro na Inicializacao de xvrp\n"); 
  */ 
  
  printf("Fazendo calculo de setup_sph\n");
  
  err=setup_sph(D,h,kh,N_sph,sph_eq,Nspecies,N,sph_neq,lbox,w_bspline,EoS);
  if(err!=0){
    printf("problemas no setup_sph err=%d\n",err);
    scanf("%d",&err);
  }
  
  for(p=0;p<argc;p+=1)
    printf("argv[%d]=%s\n",p,argv[p]);
  
  printf("Calculando densidade de energia no plano de reação\n");
  
  if(argc>=4){
   if(argc==5)
     sprintf(filename,"gfc/%s/e_react-%s.dat",argv[4],argv[3]);
   else
     sprintf(filename,"gfc/e_react-%s.dat",argv[3]);
  }
  else
   sprintf(filename,"gfc/energy-reactplane.dat");
        
  displ[0]=0.0;
  displ[1]=0.0;
  displ[2]=-1.7;
  displ[3]=0.0;
  
  if(print_proper_energy(D,Nspecies,h,kh,lbox,w.f,filename,Nrp,xrp,displ)!=0)
    printf("Erro na impressao da densidade de energia plano de reacao\n");
      
  printf("Calculando densidade de energia no plano transversal\n");
  
  if(argc>=4){
    if(argc>=5)
      sprintf(filename,"gfc/%s/e_trans-%s.dat",argv[4],argv[3]);
    else
      sprintf(filename,"gfc/e_trans-%s.dat",argv[3]);
  }
  else
    sprintf(filename,"gfc/energy-trans.dat");
    
  displ[0]=0.0;
  displ[1]=0.0;
  displ[2]=0.0;
  displ[3]=0.0;
  
  if(print_proper_energy(D,Nspecies,h,kh,lbox,w.f,filename,Npoints,xp,displ)!=0)
    printf("Erro na impressao da densidade de energia transversal\n");
    
  printf("Calculando densidade de carga no plano de reação\n");
  
  if(argc>=4){
   if(argc==5)
     sprintf(filename,"gfc/%s/c_react-%s.dat",argv[4],argv[3]);
   else
     sprintf(filename,"gfc/c_react-%s.dat",argv[3]);
  }
  else
   sprintf(filename,"gfc/charge-react.dat");
        
  displ[0]=0.0;
  displ[1]=0.0;
  displ[2]=-1.7;
  displ[3]=0.0;
  
  if(print_charge_density(D,Nspecies,h,kh,lbox,w.f,filename,Nrp,xrp,displ)!=0)
    printf("Erro na impressao da densidade de energia plano de reacao\n");
      
  printf("Calculando densidade de carga no plano transversal\n");
  
  if(argc>=4){
    if(argc>=5)
      sprintf(filename,"gfc/%s/c_trans-%s.dat",argv[4],argv[3]);
    else
      sprintf(filename,"gfc/c_trans-%s.dat",argv[3]);
  }
  else
    sprintf(filename,"gfc/charge-trans.dat");
    
  displ[0]=0.0;
  displ[1]=0.0;
  displ[2]=0.0;
  displ[3]=0.0;
  
  if(print_charge_density(D,Nspecies,h,kh,lbox,w.f,filename,Npoints,xp,displ)!=0)
    printf("Erro na impressao da densidade de energia transversal\n");
    
  for(p=0;p<Npoints;p+=1)
    free(xp[p]);
  free(xp);
  for(p=0;p<Nrp;p+=1)
    free(xrp[p]);
  free(xrp);
  /*
  for(p=0;p<Nvel;p+=1)
    free(xv[p]);
  free(xv);
  for(p=0;p<Nvelrp;p+=1)
    free(xvrp[p]);
  free(xvrp);*/
  free(N);
  free(displ);
  free(x);
    
  return 0;
}
