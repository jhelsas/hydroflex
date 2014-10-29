/*
 * Implementacao do metodo de runge-kutta para usar
 * particulas SPH.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SPHtypes.h"
#include "SPH-lista.h"
#include "SPH-state.h"

/* Checar a unidade de B depois // [du/dt ] = 1/fm = (GeV^2 /(GeV fm))/ GeV*/
/* Usando eBy = m_pi ^2 = 0.019321 ==> q_k By/(m_k hbarc) ~ 0.7055 */

/*
 * Artigo do Hirano Sensei:
 * http://arxiv.org/abs/1309.2823
 */

int eq_force(double t,int D,SPHeq_particle *sph_eqp,double *u){
  int l,m;
  double tau,E0,B0,E[D+1],B[D+1];
  if(D!=3)
    return 1;
    
  tau=1.0;
  
  B0=((0.019321*4)/hbarc)*exp(-t/tau);
  E0=B0/(2.0*tau);
  
  /*
  B0=((0.019321*4)/hbarc);
  E0=0.0;
  */
  /*
  if(sph_eqp->id==0)
    printf("Nc=%lf\nB0=%lf\nE0=%lf\nNc*B0=%lf\nrho_p=%lf\ns_p=%lf\n",(sph_eqp->Nc),B0,E0,B0*(sph_eqp->Nc),(sph_eqp->rho_p),(sph_eqp->s_p));
  */
  E[1]=-E0*(sph_eqp->x[3]);B[1]=0.0;
  E[2]=0.0;B[2]=B0;
  E[3]=E0*(sph_eqp->x[1]);B[3]=0.0;
  
  u[1] += (sph_eqp->Nc)*(E[1] + ((sph_eqp->u[2])*B[3]-(sph_eqp->u[3])*B[2])/(sph_eqp->u[0]));
  u[2] += (sph_eqp->Nc)*(E[2] + ((sph_eqp->u[3])*B[1]-(sph_eqp->u[1])*B[3])/(sph_eqp->u[0]));
  u[3] += (sph_eqp->Nc)*(E[3] + ((sph_eqp->u[1])*B[2]-(sph_eqp->u[2])*B[1])/(sph_eqp->u[0]));
  
  /*
  u[1] += (sph_eqp->Nb/50.0)*(E[1] + ((sph_eqp->u[2])*B[3]-(sph_eqp->u[3])*B[2])/(sph_eqp->u[0]));
  u[2] += (sph_eqp->Nb/50.0)*(E[2] + ((sph_eqp->u[3])*B[1]-(sph_eqp->u[1])*B[3])/(sph_eqp->u[0]));
  u[3] += (sph_eqp->Nb/50.0)*(E[3] + ((sph_eqp->u[1])*B[2]-(sph_eqp->u[2])*B[1])/(sph_eqp->u[0]));
  */
  return 0;
}

int neq_force(double t,int D,SPHneq_particle *sph_neqp,double *nu){
  int l,m;
  double tau,E0,B0,E[D+1],B[D+1];
  if(D!=3)
    return 1;
    
  tau=1.0;
  
  B0=((0.019321*0)/hbarc)*exp(-t/tau);
  E0=B0/tau;
  
  E[1]=-E0*(sph_neqp->x[3]);B[1]=0.0;
  E[2]=0.0;B[2]=B0;
  E[3]=E0*(sph_neqp->x[1]);B[3]=0.0;
  
  nu[1] += (sph_neqp->q)*(sph_neqp->N_p/(sph_neqp->m))*(E[1] + ((sph_neqp->u[2])*B[3]-(sph_neqp->u[3])*B[2])/(sph_neqp->u[0]));
  nu[2] += (sph_neqp->q)*(sph_neqp->N_p/(sph_neqp->m))*(E[2] + ((sph_neqp->u[3])*B[1]-(sph_neqp->u[1])*B[3])/(sph_neqp->u[0]));
  nu[3] += (sph_neqp->q)*(sph_neqp->N_p/(sph_neqp->m))*(E[3] + ((sph_neqp->u[1])*B[2]-(sph_neqp->u[2])*B[1])/(sph_neqp->u[0]));
  /*
  if(D==3){
    u[1] -= ((sph_eq[i].p.Nc)/((sph_eq[i].p.u[0])*hbarc))*By*(sph_eq[i].p.u[3]);
    u[2] += 0.0;
    u[3] += ((sph_eq[i].p.Nc)/((sph_eq[i].p.u[0])*hbarc))*By*(sph_eq[i].p.u[1]);
  }*/
	return 0;
}

int Deriv_MCV(double t,int D,int N_sph,int Nspecies,int *N,double h,double kh,
              SPHeq_list *sph_eq,SPHneq_list **sph_neq,Box *lbox,
              double (*w)(int,double,double),double (*Dw)(int,double,double),
              SPHeq_list *f_eq,SPHneq_list **f_neq
             )
{
  int i,k,l,lp,err;
  double dist,divv,x[D+1],u[D+1],nu[D+1],a,b,X,pforce,udotu,dS,dN,dNc,F[D+1][D+1],By;
  SPHeq_list *inter_eq,*aux_eq;
  SPHneq_list **inter_neq,*aux_neq;
  
  By=0.019321*0;
  
  inter_neq=(SPHneq_list**)malloc(Nspecies*sizeof(SPHneq_list*));
  
  for(k=0;k<Nspecies;k+=1){
    for(i=0;i<N[k];i+=1){
      
      dNc=0.0;
      f_neq[k][i].p.N_p=0.0;
      
      for(l=1;l<=D;l+=1){
        f_neq[k][i].p.x[l]=(sph_neq[k][i].p.u[l])/(sph_neq[k][i].p.u[0]);
        f_neq[k][i].p.u[l]=0.0;
        nu[l]=0.0;
      }
      
      for(l=0;l<=D;l+=1)
        x[l]=sph_neq[k][i].p.x[l];
    
      err=fetch_inter_eq(x,lbox,kh,&inter_eq);
      if(err!=0){
        printf("Problemas no fetch_inter_eq - Deriv_MCV\n");
        return err;
      }
      dN=0.0;
      while(inter_eq!=NULL){
        /*
        if((inter_eq->p).fo!=0){
          aux_eq=inter_eq;
          inter_eq=inter_eq->inext;
          aux_eq->inext=NULL;
          
          continue;
        }
        */
        
        dist=0.0;
        for(l=1;l<=D;l+=1)
          dist+=(sph_neq[k][i].p.x[l]-(inter_eq->p).x[l])*(sph_neq[k][i].p.x[l]-(inter_eq->p).x[l]);
        dist=sqrt(dist);
        
        pforce = (sph_neq[k][i].p.u[0])*((inter_eq->p).u[0]);
        
        for(l=1;l<=D;l+=1)
          pforce -= (sph_neq[k][i].p.u[l])*((inter_eq->p).u[l]);
        
        pforce = pforce*invtau_R /*+(sph_neq[k][i].p.n_p)*inveta0tau_R2*/;
        pforce *= (w(D,dist,h)*(((inter_eq->p).ni)/((inter_eq->p).rho)));
        dN+=pforce;
        
        aux_eq=inter_eq;
        inter_eq=inter_eq->inext;
        aux_eq->inext=NULL;
      }
      dN *= -((sph_neq[k][i].p.n_p)/(sph_neq[k][i].p.rho));
      f_neq[k][i].p.N_p = dN;
      
      /*
      err=neq_force(t,D,&(sph_neq[k][i].p),nu);
      if(err!=0)
        return err;
      for(l=1;l<=D;l+=1)
        f_neq[k][i].p.u[l] = nu[l];
      */
    }
  }
  /*
   ******************************************
   */
  
  for(i=0;i<N_sph;i+=1){
    /* 
    if(sph_eq[i].p.fo!=0)
     continue;
     */
     
    dS=0.0;
    pforce=0.0;
    divv=0.0;
    for(l=1;l<=D;l+=1){
      u[l] = 0.0;
      f_eq[i].p.x[l] = sph_eq[i].p.v[l];
      f_eq[i].p.u[l] = 0.0;
    }
    
    for(l=0;l<=D;l+=1)
      x[l]=sph_eq[i].p.x[l];
    
    err=fetch_inter_eq(x,lbox,kh,&inter_eq);
    if(err!=0){
      printf("Problemas no fetch_inter_eq - Deriv_MCV\n");
      return err;
      /*
      scanf("%d",&err);
      if(err!=0)
        return err;
      */
    }
    while(inter_eq!=NULL){
      /*
      if((inter_eq->p).fo!=0){
        aux_eq=inter_eq;
        inter_eq=inter_eq->inext;
        aux_eq->inext=NULL;
          
        continue;
      }
     */
      if((inter_eq->p).id==sph_eq[i].p.id){
        aux_eq=inter_eq;
        inter_eq=inter_eq->inext;
        aux_eq->inext=NULL;
        continue;
      }
      dist=0.0;
      for(l=1;l<=D;l+=1)
        dist+=(x[l]-(inter_eq->p).x[l])*(x[l]-(inter_eq->p).x[l]);
      dist=sqrt(dist);
      if(dist > 0){
        pforce=((inter_eq->p).ni)*((sph_eq[i].p.p_p/((sph_eq[i].p.rho)*(sph_eq[i].p.rho)))+((inter_eq->p).p_p/(((inter_eq->p).rho)*((inter_eq->p).rho))));
        for(l=1;l<=D;l+=1){
          divv += ((inter_eq->p).ni)*(sph_eq[i].p.v[l]-(inter_eq->p).v[l])*(Dw(D,dist,h)*((sph_eq[i].p.x[l]-(inter_eq->p).x[l])/dist));
          u[l] -= pforce*(Dw(D,dist,h))*((sph_eq[i].p.x[l]-(inter_eq->p).x[l])/dist);
        }
      }
      
      aux_eq=inter_eq;
      inter_eq=inter_eq->inext;
      aux_eq->inext=NULL;
    }
    
    for(k=0;k<Nspecies;k+=1){
      err=fetch_inter_neq(x,lbox,kh,k,&inter_neq);
      if(err!=0){
        printf("Problemas no fetch_inter_neq - %d - Deriv_MCV\n",k);
        return err;
      }
      while(inter_neq[k]!=0){
        dist=0.0;
        for(l=1;l<=D;l+=1)
          dist+=(x[l]-(inter_neq[k]->p).x[l])*(x[l]-(inter_neq[k]->p).x[l]);
        dist=sqrt(dist);
        pforce=(sph_eq[i].p.u[0])*((inter_neq[k]->p).u[0]);
        for(l=1;l<=D;l+=1)
          pforce -= (sph_eq[i].p.u[l])*((inter_neq[k]->p).u[l]);
        /*
        for(l=0;l<=D;l+=1)
          printf("%lf %lf %f\n",sph_eq[i].p.u[l],(inter_neq[k]->p).u[l],-(1./3.)*(sph_eq[i].p.u[l]));
        printf("udotu=%lf\n\n",pforce);
        */
        udotu=pforce;
        pforce = (pforce*invtau_R)*((inter_neq[k]->p).n_p)/* +(((inter_neq[k]->p).n_p)*((inter_neq[k]->p).n_p))*inveta0tau_R2 */;
        pforce*= (w(D,dist,h)/((inter_neq[k]->p).rho));
        
        dS+=pforce*((inter_neq[k]->p).m)*(((inter_neq[k]->p).ni))*udotu;        
        /* GAMBIARRA ALERT */
        /*
        for(l=1;l<=D;l+=1)
          u[l] += (4.0/3.0)*pforce*((inter_neq[k]->p).m)*(((inter_neq[k]->p).ni)/(sph_eq[i].p.rho))*((inter_neq[k]->p).u[l]); 
        */
        
        for(l=1;l<=D;l+=1)
          u[l] += pforce*((inter_neq[k]->p).m)*(((inter_neq[k]->p).ni)/(sph_eq[i].p.rho))*((inter_neq[k]->p).u[l]);
        
        /*Possibilidade:*/
        /*
        for(l=1;l<=D;l+=1)
          f_neq[k][(inter_neq[k]->p).id] -= *** (inter_neq[k]->p).ni ***  pforce*((inter_neq[k]->p).m)*((sph_eq[i].p.ni)/(sph_eq[i].p.rho))*((inter_neq[k]->p).u[l]);
        */
        /*Nao ignorar*/
        
        aux_neq=inter_neq[k];
        inter_neq[k]=inter_neq[k]->inext;
        aux_neq->inext=NULL;
      }
    }
    
    /*
     * Campo Eletro-Magnetico
     */
    /*
    if(D==3){
      u[3] += ((sph_eq[i].p.Nc)/((sph_eq[i].p.u[0])*hbarc))*By*(sph_eq[i].p.u[1]);
      u[1] -= ((sph_eq[i].p.Nc)/((sph_eq[i].p.u[0])*hbarc))*By*(sph_eq[i].p.u[3]);
      u[2] += 0.0;
    }
    */
    
    
    err=eq_force(t,D,&(sph_eq[i].p),u);
    if(err!=0)
      if(D!=1)
        return err;
    
    
    /*
     * Campo Eletro-Magnetico
     */
       
    f_eq[i].p.Nc=dNc;
    if(sph_eq[i].p.T>0)
      dS /= (sph_eq[i].p.T)*(sph_eq[i].p.rho);    /* Isso aqui tem que checar */
    /*if(dS>0)
      printf("Producao magica de entropia\n");*/
    f_eq[i].p.S=dS;
    divv=(-divv)/(sph_eq[i].p.rho);
    for(l=1;l<=D;l+=1)
      u[l] -= (dS*((sph_eq[i].p.hsh_p+sph_eq[i].p.h_p)/sph_eq[i].p.s_p) - divv*((sph_eq[i].p.hsh_p)/sph_eq[i].p.rho_p))*sph_eq[i].p.u[l];
    a=(sph_eq[i].p.h_p)/(sph_eq[i].p.rho_p);
    b=((-1.0)/((sph_eq[i].p.u[0])*(sph_eq[i].p.u[0])*(sph_eq[i].p.rho_p)))*(sph_eq[i].p.hsh_p);
    X=(-b/a)/(a+b*((sph_eq[i].p.u[0])*(sph_eq[i].p.u[0])-1.0));
    for(l=1;l<=D;l+=1){
      f_eq[i].p.u[l] = (u[l]/a);
      divv=0.0;
      for(lp=1;lp<=D;lp+=1)
        divv+=(sph_eq[i].p.u[lp])*u[lp];
      f_eq[i].p.u[l] += X*divv*(sph_eq[i].p.u[l]);
    }
  }
  
  free(inter_neq);
  return 0;
}

/*
 * Derivate of 2+1 D Boost Invariant hydrodynamics
 */

int Drv_2p1bi(double t,int D,int N_sph,int Nspecies,int *N,double h,double kh,
              SPHeq_list *sph_eq,SPHneq_list **sph_neq,Box *lbox,
              double (*w)(int,double,double),double (*Dw)(int,double,double),
              SPHeq_list *f_eq,SPHneq_list **f_neq
             )
{
  int i,k,l,lp,err;
  double dist,divv,x[D+1],u[D+1],nu[D+1],a,b,X,pforce,udotu,dS,dN,dNc,F[D+1][D+1],By;
  SPHeq_list *inter_eq,*aux_eq;
     
  for(i=0;i<N_sph;i+=1){
     
    if(sph_eq[i].p.fo!=0)
     continue;
     
    dS=0.0;
    pforce=0.0;
    divv=0.0;
    dNc=0.;
    for(l=1;l<=D;l+=1){
      u[l] = 0.0;
      f_eq[i].p.x[l] = sph_eq[i].p.v[l];
      f_eq[i].p.u[l] = 0.0;
    }
    
    for(l=0;l<=D;l+=1)
      x[l]=sph_eq[i].p.x[l];
    
    err=fetch_inter_eq(x,lbox,kh,&inter_eq);
    if(err!=0){
      printf("Problemas no fetch_inter_eq - Deriv_MCV\n");
      return err;
    }
    while(inter_eq!=NULL){
      
      if((inter_eq->p).fo!=0){
        aux_eq=inter_eq;
        inter_eq=inter_eq->inext;
        aux_eq->inext=NULL;
          
        continue;
      }
     
      if((inter_eq->p).id==sph_eq[i].p.id){
        aux_eq=inter_eq;
        inter_eq=inter_eq->inext;
        aux_eq->inext=NULL;
        continue;
      }
      dist=0.0;
      for(l=1;l<=D;l+=1)
        dist+=(x[l]-(inter_eq->p).x[l])*(x[l]-(inter_eq->p).x[l]);
      dist=sqrt(dist);
      if(dist > 0){
        pforce=((inter_eq->p).ni)*((sph_eq[i].p.p_p/((sph_eq[i].p.rho)*(sph_eq[i].p.rho)))+((inter_eq->p).p_p/(((inter_eq->p).rho)*((inter_eq->p).rho))));
        for(l=1;l<=D;l+=1){
          divv += ((inter_eq->p).ni)*(sph_eq[i].p.v[l]-(inter_eq->p).v[l])*(Dw(D,dist,h)*((sph_eq[i].p.x[l]-(inter_eq->p).x[l])/dist));
          u[l] -= t*pforce*(Dw(D,dist,h))*((sph_eq[i].p.x[l]-(inter_eq->p).x[l])/dist);
        }
      }
      
      aux_eq=inter_eq;
      inter_eq=inter_eq->inext;
      aux_eq->inext=NULL;
    }
              
    f_eq[i].p.Nc=dNc;
    if(sph_eq[i].p.T>0)
      dS /= (sph_eq[i].p.T)*(sph_eq[i].p.rho);    /* Isso aqui tem que checar */
    if(dS>0)
      printf("Producao magica de entropia\n");
    f_eq[i].p.S=dS;
    divv=(-divv)/(sph_eq[i].p.rho)+1.0/t;
    for(l=1;l<=D;l+=1)
      u[l] -= (dS*((sph_eq[i].p.hsh_p+sph_eq[i].p.h_p)/sph_eq[i].p.s_p) - divv*((sph_eq[i].p.hsh_p)/sph_eq[i].p.rho_p))*sph_eq[i].p.u[l];
    a=(sph_eq[i].p.h_p)/(sph_eq[i].p.rho_p);
    b=((-1.0)/((sph_eq[i].p.u[0])*(sph_eq[i].p.u[0])*(sph_eq[i].p.rho_p)))*(sph_eq[i].p.hsh_p);
    X=(-b/a)/(a+b*((sph_eq[i].p.u[0])*(sph_eq[i].p.u[0])-1.0));
    for(l=1;l<=D;l+=1){
      f_eq[i].p.u[l] = (u[l]/a);
      divv=0.0;
      for(lp=1;lp<=D;lp+=1)
        divv+=(sph_eq[i].p.u[lp])*u[lp];
      f_eq[i].p.u[l] += X*divv*(sph_eq[i].p.u[l]);
    }
  }
  
  return 0;
}

int HE2(int D,double t,double dt,double h,double kh,
        int N_sph, int Nspecies, int *N,
        SPHeq_list *sph_eq,SPHeq_list *sph_eqTemp,
        SPHeq_list *f0_eq,SPHeq_list *f1_eq,
        SPHneq_list **sph_neq,SPHneq_list **sph_neqTemp,
        SPHneq_list **f0_neq,SPHneq_list **f1_neq,
        Box *lbox,           
        double (*w)(int,double,double),double (*Dw)(int,double,double),
        int (*EoS)(SPHeq_particle*),
        int setup(int,double,double,double,int,SPHeq_list *,
                  int, int *, SPHneq_list **, Box *,
                  double (*)(int,double,double),int (*)(SPHeq_particle*)),
        int Deriv(double ,int ,int ,int ,int *,double ,double ,
              SPHeq_list *,SPHneq_list **,Box *,
              double (*)(int,double,double),double (*)(int,double,double),
              SPHeq_list *,SPHneq_list **),
        double *Pc
       )
{
  int i,j,k,l,err,pt0[D+1],pt1[D+1];
      
  err=Deriv(t,D,N_sph,Nspecies,N,h,kh,sph_eq,sph_neq,lbox,w,Dw,f0_eq,f0_neq);
  if(err!=0)
    return 2;
  
  /*
  for(i=0;i<N_sph;i+=1){
    sph_eq[i].p.sigma=(sph_eq[i].p.rho)*(f0_eq[i].p.S); 
    for(l=1;l<=D;l+=1)
      sph_eq[i].p.dudt[l]=f0_eq[i].p.x[l]; 
  }
  */

  for(i=0;i<N_sph;i+=1){
    for(l=1;l<=D;l+=1){
      sph_eqTemp[i].p.x[l] = sph_eq[i].p.x[l]+(dt)*(f0_eq[i].p.x[l]);
      sph_eqTemp[i].p.u[l] = sph_eq[i].p.u[l]+(dt)*(f0_eq[i].p.u[l]);
            
      f0_eq[i].p.x[l] *= dt; f0_eq[i].p.x[l]+=sph_eq[i].p.x[l];
      f0_eq[i].p.u[l] *= dt; f0_eq[i].p.u[l]+=sph_eq[i].p.u[l];
    }
    sph_eqTemp[i].p.S = sph_eq[i].p.S+(dt)*(f0_eq[i].p.S);
    f0_eq[i].p.S *=dt; f0_eq[i].p.S += sph_eq[i].p.S;
  }
  for(k=0;k<Nspecies;k+=1)
    for(i=0;i<N[k];i+=1){
      for(l=1;l<=D;l+=1){
	    sph_neqTemp[k][i].p.x[l] = sph_neq[k][i].p.x[l]+(dt)*(f0_neq[k][i].p.x[l]);
        sph_neqTemp[k][i].p.u[l] = sph_neq[k][i].p.u[l]+(dt)*(f0_neq[k][i].p.x[l]);
    
        f0_neq[k][i].p.x[l] *= dt; f0_neq[k][i].p.x[l]+=sph_neq[k][i].p.x[l];
        f0_neq[k][i].p.u[l] *= dt; f0_neq[k][i].p.u[l]+=sph_neq[k][i].p.u[l];
      }
      sph_neqTemp[k][i].p.N_p = sph_neq[k][i].p.N_p +(dt)*(f0_neq[k][i].p.N_p);
      f0_neq[k][i].p.N_p *= dt; f0_neq[k][i].p.N_p+=sph_neq[k][i].p.N_p;
	}
	
  err=setup(D,t+dt,h,kh,N_sph,sph_eqTemp,Nspecies,N,sph_neqTemp,lbox,w,EoS);
  if(err!=0){
    printf("erro no setup \nerr=%d\nD=%d\n",err,D);
    return 3;
  }
  
  err=Deriv(t+dt,D,N_sph,Nspecies,N,h,kh,sph_eqTemp,sph_neqTemp,lbox,w,Dw,f1_eq,f1_neq);
  if(err!=0)
    return 4;
  
  for(i=0;i<N_sph;i+=1){
    
    sph_eq[i].p.Sa = sph_eq[i].p.S; /* freezeout variables*/
    for(l=1;l<=D;l+=1){
      sph_eq[i].p.xa[l] = sph_eq[i].p.x[l]; /* freezeout variables */
      sph_eq[i].p.ua[l] = sph_eq[i].p.ua[l]; /* freezeout variables */
    }
    
    for(l=1;l<=D;l+=1){
      f1_eq[i].p.x[l] *= dt; f1_eq[i].p.x[l]+=sph_eq[i].p.x[l];
      f1_eq[i].p.u[l] *= dt; f1_eq[i].p.u[l]+=sph_eq[i].p.u[l];
      sph_eq[i].p.x[l] = (f0_eq[i].p.x[l]+f1_eq[i].p.x[l])/2.0;
      sph_eq[i].p.u[l] = (f0_eq[i].p.u[l]+f1_eq[i].p.u[l])/2.0;
    }
    f1_eq[i].p.S *= dt; f1_eq[i].p.S+=sph_eq[i].p.S;
    sph_eq[i].p.S = (f0_eq[i].p.S+f1_eq[i].p.S)/2.0;
  }
  
  for(k=0;k<Nspecies;k+=1)
    for(i=0;i<N[k];i+=1){
      for(l=1;l<=D;l+=1){
        f1_neq[k][i].p.x[l] *= dt; f1_neq[k][i].p.x[l]+=sph_neq[k][i].p.x[l];
        f1_neq[k][i].p.u[l] *= dt; f1_neq[k][i].p.u[l]+=sph_neq[k][i].p.u[l];
        sph_neq[k][i].p.x[l] = (f0_neq[k][i].p.x[l]+f1_neq[k][i].p.x[l])/2.0;
        sph_neq[k][i].p.u[l] = (f0_neq[k][i].p.u[l]+f1_neq[k][i].p.u[l])/2.0;
      }
      f1_neq[k][i].p.N_p *= dt; f1_neq[k][i].p.N_p+=sph_neq[k][i].p.N_p;
      sph_neq[k][i].p.N_p = (f0_neq[k][i].p.N_p+f1_neq[k][i].p.N_p)/2.0;
    }
  
  /*
   * 4-momentum Conservation Check
   */
      
  for(l=0;l<=D;l+=1){
    pt0[l] = 0.; pt1[l] = 0.;
  }
  
  for(i=0;i<N_sph;i+=1){
    pt0[0] += ((sph_eq[i].p.ni)/(sph_eq[i].p.rho))*((sph_eq[i].p.e_p+sph_eq[i].p.p_p)*(sph_eq[i].p.u[0])*(sph_eq[i].p.u[0])-sph_eq[i].p.p_p);
    pt1[0] += ((sph_eqTemp[i].p.ni)/(sph_eqTemp[i].p.rho))*((sph_eqTemp[i].p.e_p+sph_eqTemp[i].p.p_p)*(sph_eqTemp[i].p.u[0])*(sph_eqTemp[i].p.u[0])-sph_eqTemp[i].p.p_p);
    for(l=1;l<=D;l+=1){
      pt0[l] += ((sph_eq[i].p.ni)/(sph_eq[i].p.rho))*(sph_eq[i].p.e_p+sph_eq[i].p.p_p)*(sph_eq[i].p.u[0])*(sph_eq[i].p.u[l]);
      pt1[l] += ((sph_eqTemp[i].p.ni)/(sph_eqTemp[i].p.rho))*(sph_eqTemp[i].p.e_p+sph_eqTemp[i].p.p_p)*(sph_eqTemp[i].p.u[0])*(sph_eqTemp[i].p.u[l]);
    }
  }
    
  for(l=0;l<=D;l+=1)
    Pc[l] += (1./2.)*( pt0[l]/t + pt1[l]/(t+dt) )*dt;
  
  return 0;
}

int RK2(int D,double t,double dt,double h,double kh,
        int N_sph, int Nspecies, int *N,
        SPHeq_list *sph_eq,SPHeq_list *sph_eqTemp,
        SPHeq_list *f0_eq,SPHeq_list *f1_eq,
        SPHneq_list **sph_neq,SPHneq_list **sph_neqTemp,
        SPHneq_list **f0_neq,SPHneq_list **f1_neq,
        Box *lbox,           
        double (*w)(double,double),double (*Dw)(double,double),
        int (*EoS)(SPHeq_particle*),
        int setup(int,double,double,double,int,SPHeq_list *,
                  int, int *, SPHneq_list **, Box *,
                  double (*)(double,double),int (*)(SPHeq_particle*)),
        int Deriv(double ,int ,int ,int ,int *,double ,double ,
              SPHeq_list *,SPHneq_list **,Box *,
              double (*)(double,double),double (*)(double,double),
              SPHeq_list *,SPHneq_list **)
       )
{
  int i,j,k,l,err;
  
  /*printf("Calculo com sph_eq\n");*/
  err=Deriv(t,D,N_sph,Nspecies,N,h,kh,sph_eq,sph_neq,lbox,w,Dw,f0_eq,f0_neq);
  if(err!=0)
    return 2;
  
  for(i=0;i<N_sph;i+=1){
    sph_eq[i].p.sigma=(sph_eq[i].p.rho)*(f0_eq[i].p.S); /* freezeout variables */
    for(l=1;l<=D;l+=1)
      sph_eq[i].p.dudt[l]=f0_eq[i].p.x[l]; /* freezeout variables */
  }
    
  for(i=0;i<N_sph;i+=1){
    for(l=1;l<=D;l+=1){
      sph_eqTemp[i].p.x[l] = sph_eq[i].p.x[l]+(dt/2.0)*(f0_eq[i].p.x[l]);
      sph_eqTemp[i].p.u[l] = sph_eq[i].p.u[l]+(dt/2.0)*(f0_eq[i].p.u[l]);
            
      f0_eq[i].p.x[l] *= dt; f0_eq[i].p.x[l]+=sph_eq[i].p.x[l];
      f0_eq[i].p.u[l] *= dt; f0_eq[i].p.u[l]+=sph_eq[i].p.u[l];
    }
    /*Outros Updates*/
    sph_eqTemp[i].p.S = sph_eq[i].p.S+(dt/2.0)*(f0_eq[i].p.S);
    f0_eq[i].p.S *=dt; f0_eq[i].p.S += sph_eq[i].p.S;    
    /*sph_eqTemp[i].p.Nc=sph_eq[i].p.Nc+(dt/2.0)*(f0_eq[i].p.Nc);
    f0_eq[i].p.Nc *=dt; f0_eq[i].p.Nc += sph_eq[i].p.Nc;*/
    
    /*
    sph_eqTemp[i].p.N = sph_eq[i].p.N+(dt/2.0)*(f0_eq[i].p.N);
    f0_eq[i].p.N *=dt; f0_eq[i].p.N += sph_eq[i].p.N;
    */
  }
  for(k=0;k<Nspecies;k+=1)
    for(i=0;i<N[k];i+=1){
	    for(l=1;l<=D;l+=1){
	      sph_neqTemp[k][i].p.x[l] = sph_neq[k][i].p.x[l]+(dt/2.0)*(f0_neq[k][i].p.x[l]);
        sph_neqTemp[k][i].p.u[l] = sph_neq[k][i].p.u[l]+(dt/2.0)*(f0_neq[k][i].p.u[l]);
        f0_neq[k][i].p.x[l] *= dt; f0_neq[k][i].p.x[l]+=sph_neq[k][i].p.x[l];
        f0_neq[k][i].p.u[l] *= dt; f0_neq[k][i].p.u[l]+=sph_neq[k][i].p.u[l];
      }
      /*Outros updates*/
      sph_neqTemp[k][i].p.N_p = sph_neq[k][i].p.N_p +(dt/2.0)*(f0_neq[k][i].p.N_p);
      f0_neq[k][i].p.N_p *= dt; f0_neq[k][i].p.N_p+=sph_neq[k][i].p.N_p;
	  }
  
  /*printf("Calculo com sph_eqTemp\n");*/
  err=0;
  err=setup(D,t+dt/2.0,h,kh,N_sph,sph_eqTemp,Nspecies,N,sph_neqTemp,lbox,w,EoS);
  if(err!=0)
    return 3;
  err=Deriv(t+dt/2.0,D,N_sph,Nspecies,N,h,kh,sph_eqTemp,sph_neqTemp,lbox,w,Dw,f1_eq,f1_neq);
  if(err!=0)
    return 4;
  for(i=0;i<N_sph;i+=1){
    
    sph_eq[i].p.Sa = sph_eq[i].p.S; /* freezeout variables*/
    for(l=1;l<=D;l+=1){
      sph_eq[i].p.xa[l] = sph_eq[i].p.x[l]; /* freezeout variables */
      sph_eq[i].p.ua[l] = sph_eq[i].p.ua[l]; /* freezeout variables */
    }
    
    for(l=1;l<=D;l+=1){
      f1_eq[i].p.x[l] *= dt; f1_eq[i].p.x[l]+=sph_eq[i].p.x[l];
      f1_eq[i].p.u[l] *= dt; f1_eq[i].p.u[l]+=sph_eq[i].p.u[l];
      sph_eq[i].p.x[l] = f1_eq[i].p.x[l];
      sph_eq[i].p.u[l] = f1_eq[i].p.u[l];
    }
    /*Outros updates*/
    f1_eq[i].p.S *= dt; f1_eq[i].p.S+=sph_eq[i].p.S;
    sph_eq[i].p.S = f1_eq[i].p.S;    
    /*f1_eq[i].p.Nc *= dt; f1_eq[i].p.Nc += sph_eq[i].p.Nc;
    sph_eq[i].p.Nc = f1_eq[i].p.Nc;*/
    /*
    f1_eq[i].p.N *= dt; f1_eq[i].p.N+=sph_eq[i].p.N;
    sph_eq[i].p.N = (f0_eq[i].p.N+f1[i].p.N)/2.0;
     */
  }
  
  for(k=0;k<Nspecies;k+=1)
    for(i=0;i<N[k];i+=1){
      for(l=1;l<=D;l+=1){
        f1_neq[k][i].p.x[l] *= dt; f1_neq[k][i].p.x[l]+=sph_neq[k][i].p.x[l];
        f1_neq[k][i].p.u[l] *= dt; f1_neq[k][i].p.u[l]+=sph_neq[k][i].p.u[l];
        /*sph_neq[k][i].p.x[l] = (f0_neq[k][i].p.x[l]+f1_neq[k][i].p.x[l])/2.0;
        sph_neq[k][i].p.u[l] = (f0_neq[k][i].p.u[l]+f1_neq[k][i].p.u[l])/2.0;*/
        sph_neq[k][i].p.x[l] = f1_neq[k][i].p.x[l];
        sph_neq[k][i].p.u[l] = f1_neq[k][i].p.u[l];
      }
      /*Outros updates*/
      f1_neq[k][i].p.N_p *= dt; f1_neq[k][i].p.N_p+=sph_neq[k][i].p.N_p;
      sph_neq[k][i].p.N_p = f1_neq[k][i].p.N_p;
    }
  return 0;
}

int RK4(int D,double t,double dt,double h,double kh,
        int N_sph, int Nspecies, int *N,
        SPHeq_list *sph_eq,SPHeq_list *sph_eqTemp,
        SPHeq_list *f0_eq,SPHeq_list *f1_eq,
        SPHeq_list *f2_eq,SPHeq_list *f3_eq,
        SPHneq_list **sph_neq,SPHneq_list **sph_neqTemp,
        SPHneq_list **f0_neq,SPHneq_list **f1_neq,
        SPHneq_list **f2_neq,SPHneq_list **f3_neq,
        Box *lbox,           
        double (*w)(double,double),double (*Dw)(double,double),
        int (*EoS)(SPHeq_particle*),
        int setup(int,double,double,double,int,SPHeq_list *,
                  int, int *, SPHneq_list **, Box *,
                  double (*)(double,double),int (*)(SPHeq_particle*)),
        int Deriv(double ,int ,int ,int ,int *,double ,double ,
              SPHeq_list *,SPHneq_list **,Box *,
              double (*)(double,double),double (*)(double,double),
              SPHeq_list *,SPHneq_list **)
       )
{
  int i,j,k,l,err;
    
  /*err=setup_sph(D,h,kh,N_sph,sph_eq,Nspecies,N,sph_neq,lbox,w,EoS);
  if(err!=0)
    return 1;*/
    
  err=Deriv(t,D,N_sph,Nspecies,N,h,kh,sph_eq,sph_neq,lbox,w,Dw,f0_eq,f0_neq);
  if(err!=0)
    return 2;
    
  for(i=0;i<N_sph;i+=1){
    sph_eq[i].p.sigma=(sph_eq[i].p.rho)*(f0_eq[i].p.S); /* freezeout variables */
    for(l=1;l<=D;l+=1)
      sph_eq[i].p.dudt[l]=f0_eq[i].p.x[l]; /* freezeout variables */
  } 
   
  for(i=0;i<N_sph;i+=1){
    for(l=1;l<=D;l+=1){
      sph_eqTemp[i].p.x[l] = sph_eq[i].p.x[l]+(dt/2.0)*(f0_eq[i].p.x[l]);
      sph_eqTemp[i].p.u[l] = sph_eq[i].p.u[l]+(dt/2.0)*(f0_eq[i].p.u[l]);
            
      f0_eq[i].p.x[l] *= dt; f0_eq[i].p.x[l]+=sph_eq[i].p.x[l];
      f0_eq[i].p.u[l] *= dt; f0_eq[i].p.u[l]+=sph_eq[i].p.u[l];
    }
    /*Outros Updates*/
    sph_eqTemp[i].p.S = sph_eq[i].p.S+(dt/2.0)*(f0_eq[i].p.S);
    f0_eq[i].p.S *=dt; f0_eq[i].p.S += sph_eq[i].p.S;
  }
  for(k=0;k<Nspecies;k+=1)
    for(i=0;i<N[k];i+=1){
		  for(l=1;l<=D;l+=1){
	      sph_neqTemp[k][i].p.x[l] = sph_neq[k][i].p.x[l]+(dt/2.0)*(f0_neq[k][i].p.x[l]);
        sph_neqTemp[k][i].p.u[l] = sph_neq[k][i].p.u[l]+(dt/2.0)*(f0_neq[k][i].p.u[l]);
        f0_neq[k][i].p.x[l] *= dt; f0_neq[k][i].p.x[l]+=sph_neq[k][i].p.x[l];
        f0_neq[k][i].p.u[l] *= dt; f0_neq[k][i].p.u[l]+=sph_neq[k][i].p.u[l];
      }
      /*Outros updates*/
      sph_neqTemp[k][i].p.N_p = sph_neq[k][i].p.N_p +(dt/2.0)*(f0_neq[k][i].p.N_p);
      f0_neq[k][i].p.N_p *= dt; f0_neq[k][i].p.N_p+=sph_neq[k][i].p.N_p;
	  }
    
  err=setup(D,t+dt/2.0,h,kh,N_sph,sph_eqTemp,Nspecies,N,sph_neqTemp,lbox,w,EoS);
  if(err!=0)
    return 3;
  err=Deriv(t+dt/2.0,D,N_sph,Nspecies,N,h,kh,sph_eqTemp,sph_neqTemp,lbox,w,Dw,f1_eq,f1_neq);
  if(err!=0)
    return 4;
  for(i=0;i<N_sph;i+=1){
    for(l=1;l<=D;l+=1){
      sph_eqTemp[i].p.x[l] = sph_eq[i].p.x[l]+(dt/2.0)*(f1_eq[i].p.x[l]);
      sph_eqTemp[i].p.u[l] = sph_eq[i].p.u[l]+(dt/2.0)*(f1_eq[i].p.u[l]);
      f1_eq[i].p.x[l] *= dt; f1_eq[i].p.x[l]+=sph_eq[i].p.x[l];
      f1_eq[i].p.u[l] *= dt; f1_eq[i].p.u[l]+=sph_eq[i].p.u[l];
    }
    /*Outros updates*/
    sph_eqTemp[i].p.S = sph_eq[i].p.S+(dt/2.0)*(f1_eq[i].p.S);
    f1_eq[i].p.S *= dt; f1_eq[i].p.S+=sph_eq[i].p.S;
  }
  
  for(k=0;k<Nspecies;k+=1)
    for(i=0;i<N[k];i+=1){
      for(l=1;l<=D;l+=1){
	      sph_neqTemp[k][i].p.x[l] = sph_neq[k][i].p.x[l]+(dt/2.0)*(f1_neq[k][i].p.x[l]);
        sph_neqTemp[k][i].p.u[l] = sph_neq[k][i].p.u[l]+(dt/2.0)*(f1_neq[k][i].p.u[l]);
        f1_neq[k][i].p.x[l] *= dt; f1_neq[k][i].p.x[l]+=sph_neq[k][i].p.x[l];
        f1_neq[k][i].p.u[l] *= dt; f1_neq[k][i].p.u[l]+=sph_neq[k][i].p.u[l];
      }
      /*Outros updates*/
      sph_neqTemp[k][i].p.N_p = sph_neq[k][i].p.N_p +(dt/2.0)*(f1_neq[k][i].p.N_p);
      f1_neq[k][i].p.N_p *= dt; f1_neq[k][i].p.N_p+=sph_neq[k][i].p.N_p;
    }
  
    
  err=setup(D,t+dt/2.0,h,kh,N_sph,sph_eqTemp,Nspecies,N,sph_neqTemp,lbox,w,EoS);
  if(err!=0)
    return 3;
  err=Deriv(t+dt/2.0,D,N_sph,Nspecies,N,h,kh,sph_eqTemp,sph_neqTemp,lbox,w,Dw,f2_eq,f2_neq);
  if(err!=0)
    return 4;
  for(i=0;i<N_sph;i+=1){
    for(l=1;l<=D;l+=1){
      sph_eqTemp[i].p.x[l] = sph_eq[i].p.x[l]+(dt)*(f2_eq[i].p.x[l]);
      sph_eqTemp[i].p.u[l] = sph_eq[i].p.u[l]+(dt)*(f2_eq[i].p.u[l]);
      f2_eq[i].p.x[l] *= dt; f2_eq[i].p.x[l]+=sph_eq[i].p.x[l];
      f2_eq[i].p.u[l] *= dt; f2_eq[i].p.u[l]+=sph_eq[i].p.u[l];
    }
    /*Outros updates*/
    sph_eqTemp[i].p.S = sph_eq[i].p.S+(dt)*(f2_eq[i].p.S);
    f2_eq[i].p.S *= dt; f2_eq[i].p.S+=sph_eq[i].p.S;
  }
  
  for(k=0;k<Nspecies;k+=1)
    for(i=0;i<N[k];i+=1){
      for(l=1;l<=D;l+=1){
	      sph_neqTemp[k][i].p.x[l] = sph_neq[k][i].p.x[l]+(dt)*(f2_neq[k][i].p.x[l]);
        sph_neqTemp[k][i].p.u[l] = sph_neq[k][i].p.u[l]+(dt)*(f2_neq[k][i].p.u[l]);
        f2_neq[k][i].p.x[l] *= dt; f2_neq[k][i].p.x[l]+=sph_neq[k][i].p.x[l];
        f2_neq[k][i].p.u[l] *= dt; f2_neq[k][i].p.u[l]+=sph_neq[k][i].p.u[l];
      }
      /*Outros updates*/
      sph_neqTemp[k][i].p.N_p = sph_neq[k][i].p.N_p +(dt)*(f2_neq[k][i].p.N_p);
      f2_neq[k][i].p.N_p *= dt; f2_neq[k][i].p.N_p+=sph_neq[k][i].p.N_p;
    }  
    
  err=setup(D,t+dt,h,kh,N_sph,sph_eqTemp,Nspecies,N,sph_neqTemp,lbox,w,EoS);
  if(err!=0)
    return 3;
  err=Deriv(t+dt,D,N_sph,Nspecies,N,h,kh,sph_eqTemp,sph_neqTemp,lbox,w,Dw,f3_eq,f3_neq);
  if(err!=0)
    return 4;
  for(i=0;i<N_sph;i+=1){
    
    sph_eq[i].p.Sa = sph_eq[i].p.S; /* freezeout variables*/
    for(l=1;l<=D;l+=1){
      sph_eq[i].p.xa[l] = sph_eq[i].p.x[l]; /* freezeout variables */
      sph_eq[i].p.ua[l] = sph_eq[i].p.ua[l]; /* freezeout variables */
    }
    
    for(l=1;l<=D;l+=1){
      f3_eq[i].p.x[l] *= dt; f3_eq[i].p.x[l]+=sph_eq[i].p.x[l];
      f3_eq[i].p.u[l] *= dt; f3_eq[i].p.u[l]+=sph_eq[i].p.u[l];
      sph_eq[i].p.x[l] = (f0_eq[i].p.x[l]+2.*f1_eq[i].p.x[l]+2.*f2_eq[i].p.x[l]+f3_eq[i].p.x[l])/6.0;
      sph_eq[i].p.u[l] = (f0_eq[i].p.u[l]+2.*f1_eq[i].p.u[l]+2.*f2_eq[i].p.u[l]+f3_eq[i].p.u[l])/6.0;
    }
    /*Outros updates*/
    f3_eq[i].p.S *= dt; f3_eq[i].p.S+=sph_eq[i].p.S;
    sph_eq[i].p.S = (f0_eq[i].p.S+2.*f1_eq[i].p.S+2.*f2_eq[i].p.S+f3_eq[i].p.S)/6.0;
  }
  
  for(k=0;k<Nspecies;k+=1)
    for(i=0;i<N[k];i+=1){
      for(l=1;l<=D;l+=1){
        f3_neq[k][i].p.x[l] *= dt; f3_neq[k][i].p.x[l]+=sph_neq[k][i].p.x[l];
        f3_neq[k][i].p.u[l] *= dt; f3_neq[k][i].p.u[l]+=sph_neq[k][i].p.u[l];
        sph_neq[k][i].p.x[l] = (f0_neq[k][i].p.x[l]+2.*f1_neq[k][i].p.x[l]+2.*f2_neq[k][i].p.x[l]+f3_neq[k][i].p.x[l])/6.0;
        sph_neq[k][i].p.u[l] = (f0_neq[k][i].p.u[l]+2.*f1_neq[k][i].p.u[l]+2.*f2_neq[k][i].p.u[l]+f3_neq[k][i].p.u[l])/6.0;
      }
      /*Outros updates*/
      f3_neq[k][i].p.N_p *= dt; f3_neq[k][i].p.N_p+=sph_neq[k][i].p.N_p;
      sph_neq[k][i].p.N_p = (f0_neq[k][i].p.N_p+2.*f1_neq[k][i].p.N_p+2.*f2_neq[k][i].p.N_p+f3_neq[k][i].p.N_p)/6.0;
    }
  return 0;
}

int adapsize(int D,int N_sph,int Nspecies,int *N,
             double h,double kh,double t,double *dt,
             SPHeq_list *sph_eq,SPHeq_list *sph_eqTemp,
             SPHeq_list *f0_eq,SPHeq_list *f1_eq,
             SPHeq_list *f2_eq,SPHeq_list *f3_eq,
             SPHneq_list **sph_neq,SPHneq_list **sph_neqTemp,
             SPHneq_list **f0_neq,SPHneq_list **f1_neq,
             SPHneq_list **f2_neq,SPHneq_list **f3_neq,
             Box *lbox,           
             double (*w)(double,double),double (*Dw)(double,double),
             int (*EoS)(SPHeq_particle*),
             int setup(int,double,double,double,int,SPHeq_list *,
                       int, int *, SPHneq_list **, Box *,
                       double (*)(double,double),int (*)(SPHeq_particle*)),
             int Deriv(double ,int ,int ,int ,int *,double ,double ,
                       SPHeq_list *,SPHneq_list **,Box *,
                       double (*)(double,double),double (*)(double,double),
                       SPHeq_list *,SPHneq_list **)
             )
{
  int i,k,l,Nneq,err;
  double dta,dtb,errAdap,eqAdap,neqAdap,modAdap,normAdap,tolabs,tolrel;
  SPHeq_list *sph_eqAdap1,*sph_eqAdap2;
  SPHneq_list **sph_neqAdap1,**sph_neqAdap2;
 
  tolabs=0.000001;
  tolrel=0.001;
  
  err=sph_eqMalloc(N_sph,D,&sph_eqAdap1);
  if(err!=0)
    return (-30);
  
  err=sph_neqMalloc(Nspecies,N,D,&sph_neqAdap1);
  if(err!=0)
    return (-40);
  
  err=sph_eqMalloc(N_sph,D,&sph_eqAdap2);
  if(err!=0)
    return (-30);
  
  err=sph_neqMalloc(Nspecies,N,D,&sph_neqAdap2);
  if(err!=0)
    return (-40);
    
  for(i=0;i<N_sph;i+=1){
    sph_eqAdap1[i].p.ni= sph_eq[i].p.ni;
    sph_eqAdap1[i].p.S = sph_eq[i].p.S;
    for(l=0;l<=D;l+=1){
      sph_eqAdap1[i].p.x[l]=sph_eq[i].p.x[l];
      sph_eqAdap1[i].p.u[l]=sph_eq[i].p.u[l];
    }
  }
  for(i=0;i<N_sph;i+=1){
    sph_eqAdap2[i].p.ni= sph_eq[i].p.ni;
    sph_eqAdap2[i].p.S = sph_eq[i].p.S;
    for(l=0;l<=D;l+=1){
      sph_eqAdap2[i].p.x[l]=sph_eq[i].p.x[l];
      sph_eqAdap2[i].p.u[l]=sph_eq[i].p.u[l];
    }
  }
  
  for(k=0;k<Nspecies;k+=1)
    for(i=0;i<N[k];i+=1){
      sph_neqAdap1[k][i].p.m=sph_neq[k][i].p.m;
      sph_neqAdap1[k][i].p.q=sph_neq[k][i].p.q;
      sph_neqAdap1[k][i].p.ni=sph_neq[k][i].p.ni;
      sph_neqAdap1[k][i].p.N_p=sph_neq[k][i].p.N_p;
      for(l=0;l<=D;l+=1){
        sph_neqAdap1[k][i].p.x[l]=sph_neq[k][i].p.x[l];
        sph_neqAdap1[k][i].p.u[l]=sph_neq[k][i].p.u[l];
      }
    }
    
  for(k=0;k<Nspecies;k+=1)
    for(i=0;i<N[k];i+=1){
      sph_neqAdap2[k][i].p.m=sph_neq[k][i].p.m;
      sph_neqAdap2[k][i].p.q=sph_neq[k][i].p.q;
      sph_neqAdap2[k][i].p.ni=sph_neq[k][i].p.ni;
      sph_neqAdap2[k][i].p.N_p=sph_neq[k][i].p.N_p;
      for(l=0;l<=D;l+=1){
        sph_neqAdap2[k][i].p.x[l]=sph_neq[k][i].p.x[l];
        sph_neqAdap2[k][i].p.u[l]=sph_neq[k][i].p.u[l];
      }
    }
    
  Nneq=0;
  for(k=0;k<Nspecies;k+=1)
    Nneq+=N[k];
    
  dta=*dt;
  do{
	
	printf("dta=%lf\n",dta);
	  
    for(i=0;i<N_sph;i+=1){
      sph_eqAdap1[i].p.S = sph_eq[i].p.S;
      for(l=0;l<=D;l+=1){
        sph_eqAdap1[i].p.x[l]=sph_eq[i].p.x[l];
        sph_eqAdap1[i].p.u[l]=sph_eq[i].p.u[l];
      }
    }

    for(i=0;i<N_sph;i+=1){
      sph_eqAdap2[i].p.S = sph_eq[i].p.S;
      for(l=0;l<=D;l+=1){
        sph_eqAdap2[i].p.x[l]=sph_eq[i].p.x[l];
        sph_eqAdap2[i].p.u[l]=sph_eq[i].p.u[l];
      }
    }
  
    for(k=0;k<Nspecies;k+=1)
      for(i=0;i<N[k];i+=1){
        sph_neqAdap1[k][i].p.N_p=sph_neq[k][i].p.N_p;
        for(l=0;l<=D;l+=1){
          sph_neqAdap1[k][i].p.x[l]=sph_neq[k][i].p.x[l];
          sph_neqAdap1[k][i].p.u[l]=sph_neq[k][i].p.u[l];
        }
      }
      
    for(k=0;k<Nspecies;k+=1)
      for(i=0;i<N[k];i+=1){
        sph_neqAdap2[k][i].p.N_p=sph_neq[k][i].p.N_p;
        for(l=0;l<=D;l+=1){
          sph_neqAdap2[k][i].p.x[l]=sph_neq[k][i].p.x[l];
          sph_neqAdap2[k][i].p.u[l]=sph_neq[k][i].p.u[l];
        }
      }
     
    printf("Calculando Adap1\n");
     
    err=setup(D,t,h,kh,N_sph,sph_eqAdap1,Nspecies,N,sph_neqAdap1,lbox,w,EoS);
    if(err!=0)
      return err;
    err=RK2(D,t,dta,h,kh,N_sph,Nspecies,N,sph_eqAdap1,sph_eqTemp,f0_eq,f1_eq,sph_neqAdap1,sph_neqTemp,f0_neq,f1_neq,lbox,w,Dw,EoS,setup,Deriv);
    if(err!=0)
      return err;
  
    printf("Calculando Adap2 - 1\n");
    dtb=dta/2;
    err=setup(D,t,h,kh,N_sph,sph_eqAdap2,Nspecies,N,sph_neqAdap2,lbox,w,EoS);
    if(err!=0)
      return err;
    err=RK2(D,t,dtb,h,kh,N_sph,Nspecies,N,sph_eqAdap2,sph_eqTemp,f0_eq,f1_eq,sph_neqAdap2,sph_neqTemp,f0_neq,f1_neq,lbox,w,Dw,EoS,setup,Deriv);
    if(err!=0)
      return err;
    
    printf("Calculando Adap2 - 2\n");
    err=setup(D,t+dtb,h,kh,N_sph,sph_eqAdap2,Nspecies,N,sph_neqAdap2,lbox,w,EoS);
    if(err!=0)
      return err;
    err=RK2(D,t+dtb,dtb,h,kh,N_sph,Nspecies,N,sph_eqAdap2,sph_eqTemp,f0_eq,f1_eq,sph_neqAdap2,sph_neqTemp,f0_neq,f1_neq,lbox,w,Dw,EoS,setup,Deriv);
    if(err!=0)
      return err;
    
    printf("Checando Erro\n");
    
    eqAdap=0.0;
    for(i=0;i<N_sph;i+=1){
      eqAdap+=fabs(sph_eqAdap2[i].p.S-sph_eqAdap1[i].p.S);
      
      modAdap=0.0;
      for(l=1;l<=D;l+=1)
        modAdap+=(sph_eqAdap2[i].p.u[l]-sph_eqAdap1[i].p.u[l])*(sph_eqAdap2[i].p.u[l]-sph_eqAdap1[i].p.u[l]);
      eqAdap+=sqrt(modAdap);      
      
      modAdap=0.0;
      for(l=1;l<=D;l+=1)
        modAdap+=(sph_eqAdap2[i].p.x[l]-sph_eqAdap1[i].p.x[l])*(sph_eqAdap2[i].p.x[l]-sph_eqAdap1[i].p.x[l]);
      eqAdap+=(m_pi/hbarc)*sqrt(modAdap);
    }
    if(N_sph>0)
      eqAdap/=(double)N_sph;
    
    neqAdap=0.0;
    for(k=0;k<Nspecies;k+=1)
      for(i=0;i<N[k];i+=1){
	    neqAdap+=fabs(sph_neqAdap2[k][i].p.N_p-sph_neqAdap1[k][i].p.N_p);
         
        modAdap=0.0;
        for(l=1;l<=D;l+=1)
          modAdap+=(sph_neqAdap2[k][i].p.u[l]-sph_neqAdap1[k][i].p.u[l])*(sph_neqAdap2[k][i].p.u[l]-sph_neqAdap1[k][i].p.u[l]);
        neqAdap+=sqrt(modAdap);
        
        modAdap=0.0;
        for(l=1;l<=D;l+=1)
          modAdap+=(sph_neqAdap2[k][i].p.x[l]-sph_neqAdap1[k][i].p.x[l])*(sph_neqAdap2[k][i].p.x[l]-sph_neqAdap1[k][i].p.x[l]);
        neqAdap+=(m_pi/hbarc)*sqrt(modAdap);
      }
    if(Nneq>0)
      neqAdap/=(double)Nneq;
    
    errAdap=eqAdap+neqAdap;
    
    eqAdap=0.0;
    for(i=0;i<N_sph;i+=1){
      eqAdap+=fabs(sph_eq[i].p.S);
      
      modAdap=0.0;
      for(l=1;l<=D;l+=1)
        modAdap+=(sph_eq[i].p.u[l])*(sph_eq[i].p.u[l]);
      eqAdap+=sqrt(modAdap);
      
      modAdap=0.0;
      for(l=1;l<=D;l+=1)
        modAdap+=(sph_eq[i].p.x[l])*(sph_eq[i].p.x[l]);
      eqAdap+=(m_pi/hbarc)*sqrt(modAdap);
    }
    if(N_sph>0)
      eqAdap/=(double)N_sph;
    
    neqAdap=0.0;
    for(k=0;k<Nspecies;k+=1)
      for(i=0;i<N[k];i+=1){
        neqAdap+=fabs(sph_neq[k][i].p.N_p);
        
        modAdap=0.0;
        for(l=1;l<=D;l+=1)
          modAdap+=(sph_neq[k][i].p.u[l])*(sph_neq[k][i].p.u[l]);
        neqAdap=sqrt(modAdap);
        
        modAdap=0.0;
        for(l=1;l<=D;l+=1)
          modAdap+=(sph_neq[k][i].p.x[l])*(sph_neq[k][i].p.x[l]);
        neqAdap+=(m_pi/hbarc)*modAdap;
      }
    if(N_sph>0)  
      neqAdap/=(double)Nneq;    
    normAdap=eqAdap+neqAdap;
        
    dtb=dta;
    dta=dta/2;
  }while(errAdap>tolabs+tolrel*normAdap);
    
  sph_eqFree(N_sph,&sph_eqAdap1);
  sph_neqFree(Nspecies,N,&sph_neqAdap1);
  sph_eqFree(N_sph,&sph_eqAdap2);
  sph_neqFree(Nspecies,N,&sph_neqAdap2);
  
  *dt=dtb;
  
  return 0;
}
