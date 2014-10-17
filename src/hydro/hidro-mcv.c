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

int main(int argc,char **argv)
{
  int count,i,l,p,k,D,N_sph,Nspecies,Nneq,err,Npoints,Nlong,Nvel,Nrp,Nvelrp,Kf,Idf,checkF,gfcf;
  int *N;
  const int scriptsize=6;
  char **scriptinfo;
  double t,t0,tf,dt,h,kh,dist,s,ni,S,a,H,Hp,E_T,u0sqr;
  double uTF,gammaF,wF,uzF,BF,phi0,sgn;
  double *x,*displ;
  double *Pc,step,graph_size;
  double **xp,**xv,**xrp,**xvrp,*xi,*ui,*xl,*xu,*dx,*xP,*xG,*xP2;
  SPHeq_list *sph_eq,*sph_eqTemp,*f0_eq,*f1_eq,*f2_eq,*f3_eq,*aux_eq,*inter_eq;
  double omega,L,lambda,xt;
  SPHneq_list **sph_neq,**sph_neqTemp,**f0_neq,**f1_neq,**f2_neq,**f3_neq,*aux_neq,**inter_neq;
  Box *lbox;
  char config_file[100+1],filename[100+1];
  funDDtD w,Dw;
  FILE *dadosin,*dadosout,*dadosH,*dadosF,*script_t,*script_r;
  struct stat file_stat;
  
  omega=0.1;
  L=1.0;
  
  scriptinfo=(char**)malloc((scriptsize+1)*sizeof(char*));
  for(i=0;i<=scriptsize;i+=1)
    scriptinfo[i]=(char*)malloc((100+1)*sizeof(char));
  
  checkF=0;
  gfcf=1;
  
  printf("Inicializacao do Programa\n");
  
  w.f=w_bspline;
  Dw.f=Dw_bspline;
  
  if(argc>=2)
    strcpy(config_file,argv[1]);
  else
    sprintf(config_file,"cfg/hidro-gold_RHIC.cfg");
  
  if(argc>=3)
    printf("hashname = %s\n",argv[2]);
  
  if(argc>=4){
    printf("directory name = %s\n",argv[3]);
    sprintf(filename,"gfc/%s",argv[3]);
    if(mkdir(filename, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)!=0)
      printf("Diretorio de graficos ja existe\n");
    else
      printf("Diretorio de particulas criado com sucesso\n");
      
    sprintf(filename,"data/%s",argv[3]);
    if(mkdir(filename, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)!=0)
      printf("Diretorio de particulas ja existe\n");
    else
      printf("Diretorio de particulas criado com sucesso\n");
  }
  
  err=init_simulation(config_file,&D,&N_sph,&Nspecies,&N,&h,&kh,&t0,&tf,&dt,&sph_eq,&sph_eqTemp,&f0_eq,&f1_eq,&f2_eq,&f3_eq,&sph_neq,&sph_neqTemp,&f0_neq,&f1_neq,&f2_neq,&f3_neq,&lbox);
  if(err!=0){
    printf("Initialization Error - %d\n",err);
    return err;
  }
    
  printf("D=%d\nNspecies=%d N_sph=%d\nh=%lf kh=%lf\nt0=%.12lf ,tf=%.12lf ,dt=%.12lf\n",D,Nspecies,N_sph,h,kh,t0,tf,dt);
    
  Pc=(double*)malloc((D+1)*sizeof(double));
  
  /*
  err=adapsize(D,N_sph,Nspecies,N,h,kh,t,&dt,sph_eq,sph_eqTemp,f0_eq,f1_eq,f2_eq,f3_eq,sph_neq,sph_neqTemp,f0_neq,f1_neq,f2_neq,f3_neq,lbox,w.f,Dw.f,EoS_qg,Deriv_MCV);
  if(err!=0)
    printf("Erro no adapsize - %d\n",err);
  printf("Terminada Adaptacao, dt=%lf\n",dt);
  */
  
  printf("oi-1\n");
    
  inter_neq=(SPHneq_list**)malloc(Nspecies*sizeof(SPHneq_list*));
  if(inter_neq==NULL)
    printf("problemas alocacao inter_neq\n");
    
  x=(double*)malloc((D+1)*sizeof(double));
  if(x==NULL)
    printf("problemas alocacao x\n");
  
  displ=(double*)malloc((D+1)*sizeof(double));
  if(displ==NULL)
    printf("problemas alocacao displ\n");  
  
  xl=(double *)malloc((D+1)*sizeof(double));
  xu=(double*)malloc((D+1)*sizeof(double));
  dx=(double*)malloc((D+1)*sizeof(double));
  for(l=0;l<=D;l+=1){xl[l]=-9.0;dx[l]=0.02;xu[l]=9.0+1.01*dx[l];}
  printf("antes\n");
  err=create_grid(D,&xG,xl,xu,dx,&Npoints);if(err!=0){ return err;}
  err=err=create_X_line(D,&xP,xl,xu,dx,&Npoints);if(err!=0){ return err;};
  err=err=create_Y_line(D,&xP2,xl,xu,dx,&Npoints);if(err!=0){ return err;};
  printf("logo depois\n");
  printf("Npoints=%d\n",Npoints);
  free(xl); free(xu); free(dx);
  
  /*
  sprintf(scriptinfo[0],argv[0]);
  sprintf(scriptinfo[1],argv[2]);
  sprintf(scriptinfo[2],"1");
  sprintf(scriptinfo[3],"0");
  sprintf(scriptinfo[4],"0");
  sprintf(scriptinfo[5],"0");
  sprintf(scriptinfo[6],argv[3]);
    
  err=scripting(scriptsize+1,scriptinfo,t0,tf,dt);
  */   
  
  printf("Inicio do Loop Temporal\n");
  
  dadosH=fopen("gfc/Deposit_check.dat","w");
  count=0;
  E_T=0.0;
  for(l=0;l<=D;l+=1)
    Pc[l]=0.;
  H=0.;
  Hp=0.;
  for(t=t0;t<=tf;t+=dt){
    /*printf("loop - t=%.12f\n",t);*/
    if(count%4==0)
      printf("tau = %.2f\n",t);
  
    err=ssph_2p1bi(D,t,h,kh,N_sph,sph_eq,Nspecies,N,sph_neq,lbox,w_bspline,EoS);
    if(err!=0){
      printf("problemas no setup_sph - t=%f, err=%d\n",t,err);
      scanf("%d",&err);
      if(err==99)
        break;
    }
        
    H=0.0;
    for(i=0;i<N_sph;i+=1){
      u0sqr=(sph_eq[i].p.u[0])*(sph_eq[i].p.u[0]);
      H += ((sph_eq[i].p.ni)/(sph_eq[i].p.rho))*((sph_eq[i].p.e_p+sph_eq[i].p.p_p)*u0sqr-sph_eq[i].p.p_p);
    }
    if(count > 1)        
      fprintf(dadosH,"%lf %lf %lf %lf %lf %lf\n",t,t*H,H+Pc[0],H+E_T,H, (H-Hp)/dt + (Hp/(t-dt)+H/t)/2.);
    fflush(dadosH);
    E_T += H*(1./t)*dt;
    
    Hp=H;
    
    if( (count%40)==0 || 
        count == 0 || count == 8 || count == 16 || 
        count == 32 || count == 48 || count == 56){
      printf("print in - t=%.12f\n",t);
      if(argc==3)
        sprintf(filename,"data/eqp-%s-(t=%.1lf).dat",argv[2],t);
      else if(argc==4)
        sprintf(filename,"data/%s/eqp-%s-(t=%.1lf).dat",argv[3],argv[2],t);
      else
        sprintf(filename,"data/sph_eqparticles-(t=%lf).dat",t);
      dadosout=fopen(filename,"w");
      if(dadosout==NULL){
        printf("Problemas em abrir %s\n",filename);
        return (-1);
      }
      for(i=0;i<N_sph;i+=1){        
        fprintf(dadosout,"%f %f %f ",sph_eq[i].p.ni,sph_eq[i].p.q,sph_eq[i].p.S);
        for(l=1;l<=D;l+=1)
          fprintf(dadosout,"%f ",sph_eq[i].p.x[l]);
        for(l=1;l<=D;l+=1)
          fprintf(dadosout,"%f ",sph_eq[i].p.u[l]);
          
        fprintf(dadosout,"%lf ",sph_eq[i].p.e_p);
         
        fprintf(dadosout,"\n");
      }
      fclose(dadosout);
      
      for(l=0;l<=D;l+=1)
        displ[l]=0.0;
      
      if(argc>=3){
        if(argc==4)
          sprintf(filename,"gfc/%s/Xaxis_dens-%s-(t=%.1lf).dat",argv[3],argv[2],t);
        else
          sprintf(filename,"gfc/Xaxis_dens-%s-(t=%.1lf).dat",argv[2],t);
      }
      else
        sprintf(filename,"gfc/X-edens-(t=%lf).dat",t);
        
      if(print_new(D,Nspecies,t,h,kh,lbox,w.f,filename,Npoints,xP,displ)!=0)
        printf("Erro na impressao da densidade de energia no eixo X\n");
        
      if(argc>=3){
        if(argc==4)
          sprintf(filename,"gfc/%s/Yaxis_dens-%s-(t=%.1lf).dat",argv[3],argv[2],t);
        else
          sprintf(filename,"gfc/Yaxis_dens-%s-(t=%.1lf).dat",argv[2],t);
      }
      else
        sprintf(filename,"gfc/Y-edens-(t=%lf).dat",t);
        
      if(print_new(D,Nspecies,t,h,kh,lbox,w.f,filename,Npoints,xP2,displ)!=0)
        printf("Erro na impressao da densidade de energia no eixo Y\n");
        
      if(argc>=3){
        if(argc==4)
          sprintf(filename,"gfc/%s/edens-%s-(t=%.1lf).dat",argv[3],argv[2],t);
        else
          sprintf(filename,"gfc/edens-%s-(t=%.1lf).dat",argv[2],t);
      }
      else
        sprintf(filename,"gfc/grid-edens-(t=%lf).dat",t);
        
      if(print_new(D,Nspecies,t,h,kh,lbox,w.f,filename,Npoints,xG,displ)!=0)
        printf("Erro na impressao da densidade de energia na malha\n");
        
      printf("print out - t=%.12f\n",t);
      
    }
    
    /*
    if(count==20){
      printf("Fazendo o calculo do freezeout a tempo constante\n");
      if(argc==2)
        err=freezeout_const_time(D,N_sph,sph_eq,20,10,90,0,"","");
      else if(argc==3)
        err=freezeout_const_time(D,N_sph,sph_eq,20,10,90,1,argv[2],"");
      else{
        sprintf(filename,"data/%s",argv[3]);
        err=freezeout_const_time(D,N_sph,sph_eq,20,10,90,2,argv[2],filename);
      }
      if(err!=0)
        printf("D nao e 3\n");
    }*/
    
    
    /*
    for(l=0;l<=D;l+=1)
      Pc[l]=0.0;
      
    for(i=0;i<N_sph;i+=1){
      u0sqr=(sph_eq[i].p.u[0])*(sph_eq[i].p.u[0]);
      Pc[0]+=((sph_eq[i].p.ni)/(sph_eq[i].p.rho))*((sph_eq[i].p.e_p+sph_eq[i].p.p_p)*u0sqr-sph_eq[i].p.p_p);
      for(l=1;l<=D;l+=1)
        Pc[l] += ((sph_eq[i].p.ni)/(sph_eq[i].p.rho))*(sph_eq[i].p.e_p+sph_eq[i].p.p_p)*(sph_eq[i].p.u[0])*(sph_eq[i].p.u[l]);
    }
      
    fprintf(dadosH,"%.12lf ",t);
    for(l=0;l<=D;l+=1)
      fprintf(dadosH,"%.12lf ",Pc[l]);
    */
    
    /*
    fprintf(dadosH,"%.12lf %.12lf ",sph_neq[0][0].p.N_p,sph_eq[0].p.S);
    fprintf(dadosH,"%.12lf %.12lf %.12lf %.12lf",sph_neq[0][0].p.x[3],sph_eq[0].p.x[3],sph_neq[0][0].p.u[3],sph_eq[0].p.u[3]);
    
    fprintf(dadosH,"\n");*/
    
    /*err=RK4(D,t,dt,h,kh,N_sph,Nspecies,N,sph_eq,sph_eqTemp,f0_eq,f1_eq,f2_eq,f3_eq,sph_neq,sph_neqTemp,f0_neq,f1_neq,f2_neq,f3_neq,lbox,w.f,Dw.f,EoS,ssph_2p1bi,Drv_2p1bi);*/
    /*err=RK2(D,t,dt,h,kh,N_sph,Nspecies,N,sph_eq,sph_eqTemp,f0_eq,f1_eq,sph_neq,sph_neqTemp,f0_neq,f1_neq,lbox,w.f,Dw.f,EoS,ssph_2p1bi,Drv_2p1bi);*/
    err=HE2(D,t,dt,h,kh,N_sph,Nspecies,N,sph_eq,sph_eqTemp,f0_eq,f1_eq,sph_neq,sph_neqTemp,f0_neq,f1_neq,lbox,w.f,Dw.f,EoS,ssph_2p1bi,Drv_2p1bi,Pc);
  
    if(err!=0){
      printf("problemas no Integrador - t=%f err=%d\n",t,err);return err;}
    count+=1;
  }
  
  /*
  if(argc>=3){
    fclose(script_t);
    fclose(script_r);
  }*/
  
  fclose(dadosH);dadosH=NULL;
  if(checkF!=0)
    fclose(dadosF);
  dadosF=NULL;
  printf("Fim do Loop Tempoaral\n");
  
  
  printf("Inicio das desalocacoes\n");
  
  sph_eqFree(N_sph,&sph_eq);
  sph_eqFree(N_sph,&sph_eqTemp);
  sph_eqFree(N_sph,&f0_eq);
  sph_eqFree(N_sph,&f1_eq);
  sph_eqFree(N_sph,&f2_eq);
  sph_eqFree(N_sph,&f3_eq);
  
  sph_neqFree(Nspecies,N,&sph_neq);
  sph_neqFree(Nspecies,N,&sph_neqTemp);
  sph_neqFree(Nspecies,N,&f0_neq);
  sph_neqFree(Nspecies,N,&f1_neq);
  sph_neqFree(Nspecies,N,&f2_neq);
  sph_neqFree(Nspecies,N,&f3_neq);
  
  for(k=0;k<Nspecies;k+=1)
    free(lbox->neqbox[k]);
  free(lbox->neqbox);lbox->neqbox=NULL;
  free(lbox->box);lbox->box=NULL;
  free(lbox->k);lbox->k=NULL;
  free(lbox->ind);lbox->ind=NULL;
  free(lbox->xmin);lbox->xmin=NULL;
  free(lbox->xmax);lbox->xmax=NULL;
  free(lbox->Lbox);lbox->Lbox=NULL;
  free(lbox->dx);lbox->dx=NULL;
  free(lbox);lbox=NULL;
  free(inter_neq);
  
  /*
  for(p=0;p<Npoints;p+=1)
    free(xp[p]);
  free(xp);
  for(p=0;p<Nrp;p+=1)
    free(xrp[p]);
  free(xrp);
  for(p=0;p<Nvel;p+=1)
    free(xv[p]);
  free(xv);
  for(p=0;p<Nvelrp;p+=1)
    free(xvrp[p]);
  free(xvrp);
  free(N);*/
  free(displ);
  free(x);
  free(Pc);
    
  printf("Fim do Programa\n");
    
  return 0;
}

