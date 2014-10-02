/*
 * Arquivo com as diferentes equacoes de estado que podem ser utilizadas
 * para o calculo das simulacoes.
 */

#include <stdio.h>
#include <math.h>
#include "SPHtypes.h"
#include "SPH-state.h"

int EoS_landau(SPHeq_particle *par)
{
	par->p_p=C_pi*pow(par->s_p,4.0/3.0);
	par->e_p=3.0*(par->p_p);
	par->h_p=(par->e_p)+(par->p_p);
	par->hsh_p=(4.0/3.0)*(par->p_p);
  par->T=((par->h_p)/(par->s_p));
	
	return 0;
}

int EoS_qg(SPHeq_particle *par)
{
	par->p_p=C_qg*pow(par->s_p,4.0/3.0);
	par->e_p=3.0*(par->p_p);
	par->h_p=(par->e_p)+(par->p_p);
	par->hsh_p=(4.0/3.0)*(par->p_p);
  par->T=((par->h_p)/(par->s_p));/*Isso tem que ser verificado -- eu to supondo todo mundo como sendo positivamente carregado*/
  par->nc = (par->q)*EtaCharg_qg*(par->s_p); 
  par->Nc = (par->nc)/(par->rho);
	
	return 0;
}

#define m_ref (0.139)
#define EtaCharg_dust 1.0

int EoS_dust(SPHeq_particle *par)
{
	par->p_p=0.0;
	par->e_p=m_ref*(par->S)*(par->rho_p);
	par->h_p=(par->e_p)+(par->p_p);
	par->hsh_p=(m_ref)*(par->s_p)-par->h_p;
  par->T=0;
  par->nc = /*EtaCharg_qg* */(par->s_p); 
  par->Nc = (par->nc)/(par->rho_p);
	
	return 0;
}

int freezeout_const_time(int D,int N_sph,SPHeq_list *sph_eq,int Ny,int NpT,int Nphi,
                         int namingcase,char *hash,char *filepath){
  int i,k,l,m,j,spinstat;
  const double MY_PI=3.1415926;
  double m_pion,expy;
  double ymin,ymax,y,dy;
  double pT,dpT,pTmin,pTmax;
  double phi, dphi;
  double P[D+1],g;
  double dNd3p[Nphi+1][NpT+1][Ny+1],dNdphidpT[Nphi+1][NpT+1],dNdphi[Nphi+1],dNpTdpT[NpT+1];
  char filename[150+1];
  FILE *phi_dist,*d3p_dist;
      
  if(namingcase==0){
    phi_dist=fopen("angular_distribution.dat","w");
    d3p_dist=fopen("full_p3_distribution.dat","w");
  }
  else if(namingcase ==1){
    printf("hash=%s\n",hash);
    sprintf(filename,"angular_distribution-%s.dat",hash);
    phi_dist=fopen(filename,"w");
    if(phi_dist==0)
      printf("problemas ao abrir phi_dist\n");
    sprintf(filename,"full_p3_distribution-%s.dat",hash);
    d3p_dist=fopen(filename,"w");    
    if(phi_dist==0)
      printf("problemas ao abrir d3p_dist\n");
  }
  else if(namingcase ==2){
    sprintf(filename,"%s/angular_distribution-%s.dat",filepath,hash);
    phi_dist=fopen(filename,"w");
    if(phi_dist==0)
      printf("problemas ao abrir phi_dist\n");
    sprintf(filename,"%s/full_p3_distribution-%s.dat",filepath,hash);
    d3p_dist=fopen(filename,"w");  
    if(phi_dist==0)
      printf("problemas ao abrir d3p_dist\n");
  }
  spinstat=0;
  
  if(D!=3)
    return 1;
  
  g=3.0;
  
  ymin=-0.5;
  ymax=0.5;
  pTmin=0.0;
  pTmax=5.0;
    
  dphi=(2.0*MY_PI)/((double)(Nphi));
  dpT=(pTmax-pTmin)/((double)(NpT));
  dy=(ymax-ymin)/((double)Ny);
    
  for(k=0;k<=Nphi;k+=1){
    dNdphi[k]=0.0;
    phi=((double)k)*dphi;
    
    printf("phi=%lf\n",phi);
    for(l=0;l<=NpT;l+=1){
      dNdphidpT[k][l]=0.0;
      pT=pTmin+((double)l)*dpT;
      
      for(m=0;m<=Ny;m+=1){     
        dNd3p[k][l][m]=0.0;
        y=ymin+((double)m)*dy;
        
        P[0]=sqrt(m_pi*m_pi+pT*pT)*cosh(y);
        P[1]=pT*cos(phi);
        P[2]=pT*sin(phi);
        P[3]=sqrt(m_pi*m_pi+pT*pT)*sinh(y);
        
        for(i=0;i<N_sph;i+=1){
          expy=P[0]*(sph_eq[i].p.u[0]);
          for(j=1;j<=D;j+=1)
            expy -= P[j]*(sph_eq[i].p.u[j]);
          
          if(spinstat==0)
            dNd3p[k][l][m]+= ((sph_eq[i].p.ni)/(sph_eq[i].p.rho))/(exp(expy/(sph_eq[i].p.T))-1);
          else
            dNd3p[k][l][m]+= ((sph_eq[i].p.ni)/(sph_eq[i].p.rho))/(exp(expy/(sph_eq[i].p.T))+1);
        }
        dNd3p[k][l][m]*=(g*P[0])/((2.0*MY_PI)*(2.0*MY_PI)*(2.0*MY_PI));
        
        fprintf(d3p_dist,"%lf %lf %lf %.12lf\n",phi,pT,y,dNd3p[k][l][m]);
        
        if(m==0 || m==Ny)
          dNdphidpT[k][l]+=(dNd3p[k][l][m]);
        else
          dNdphidpT[k][l]+= 4.0*dNd3p[k][l][m];
      }
      dNdphidpT[k][l]*=dy/6.0;
      
      if(l==0 || m==NpT)
        dNdphi[k]+=(pT*dNdphidpT[k][l]);
      else
        dNdphi[k]+= 4.0*pT*dNdphidpT[k][l];
    }
    dNdphi[k]*=dpT/6.0;
    
    fprintf(phi_dist,"%lf %.12lf\n",phi,dNdphi[k]);
    fflush(phi_dist);
  }
  
  /*
  for(l=0;l<=NpT;l+=1){
    for(k=0;k<=Nphi;k+=1){
      for(m=0;m<=Ny;m+=1){
        
      }
    }
  }*/
    
  fclose(phi_dist);
  fclose(d3p_dist);
  
  return 0;
}
