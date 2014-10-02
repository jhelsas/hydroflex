#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "SPHtypes.h"
#include "SPH-state.h"
 
double maxD(double *a,int D){
  double max;
  int i;
  max=a[0];
  for(i=0;i<D;i+=1)
    if(a[i]>max)
      max=a[i];
  return max;
}

double minD(double *a,int D){
  double min;
  int i;
  min=a[0];
  for(i=0;i<D;i+=1)
    if(a[i]<min)
      min=a[i];
  return min;
}

/************************************************/

int SPHeq_list_blen(SPHeq_list *lis){
  int len=0;
  while(lis!=NULL){
    len+=1;
    lis=lis->bnext;
  }
  return len;
}

int SPHeq_list_ilen(SPHeq_list *lis){
  int len=0;
  while(lis!=NULL){
    len+=1;
    lis=lis->inext;
  }
  return len;
}

int SPHeq_append_bhead(SPHeq_list **lis,SPHeq_list *el){
  SPHeq_list *aux;
  if(el==NULL)
    return 1;
  aux=el;
  while(aux->bnext != NULL)
    aux=aux->bnext;
  aux->bnext = *lis;
  *lis = el;
  return 0;
}

int SPHeq_append_ihead(SPHeq_list **lis,SPHeq_list *el){
  SPHeq_list *aux;
  if(el==NULL)
    return 1;
  aux=el;
  while(aux->inext != NULL)
    aux=aux->inext;
  aux->inext = *lis;
  *lis = el;
  return 0;
}

int SPHeq_clean_list(int N_sph,SPHeq_list *lsph){
  int i;
  
  if(lsph==NULL && N_sph>0)
    return 1;
  for(i=0;i<N_sph;i+=1){
    lsph[i].bnext=NULL;
    lsph[i].inext=NULL;
  }
  return 0;
}

/****************************************************/

int SPHneq_list_blen(SPHneq_list *lis){
  int len=0;
  while(lis!=NULL){
    len+=1;
    lis=lis->bnext;
  }
  return len;
}

int SPHneq_list_ilen(SPHneq_list *lis){
  int len=0;
  while(lis!=NULL){
    len+=1;
    lis=lis->inext;
  }
  return len;
}

int SPHneq_list_len(SPHneq_list *lis){
  int len=0;
  while(lis!=NULL){
    len+=1;
    lis=lis->next;
  }
  return len;
}

int SPHneq_append_bhead(SPHneq_list **lis,SPHneq_list *el){
  SPHneq_list *aux;
  if(el==NULL)
    return 1;
  aux=el;
  while(aux->bnext != NULL)
    aux=aux->bnext;
  aux->bnext = *lis;
  *lis = el;
  return 0;
}

int SPHneq_append_ihead(SPHneq_list **lis,SPHneq_list *el){
  SPHneq_list *aux;
  if(el==NULL)
    return 1;
  aux=el;
  while(aux->inext != NULL)
    aux=aux->inext;
  aux->inext = *lis;
  *lis = el;
  return 0;
}

int SPHneq_clean_list(int Nspecies,int *N,SPHneq_list **sph_neq){
  int k,i;
  if(sph_neq==NULL)
    return 2;
  for(k=0;k<Nspecies;k+=1){
    if(sph_neq[k]==NULL && N[k]>0)
      return 1;
	for(i=0;i<N[k];i+=1){
      sph_neq[k][i].bnext=NULL;
      sph_neq[k][i].inext=NULL; 
	}
  }
  return 0;
}

/****************************************************/

int clean_box(Box *lbox){
  int k,i;
  
  if(lbox->Nbox<=0)
    return 1;
  for(i=0;i<lbox->Nbox;i+=1)
    lbox->box[i]=NULL;    
  for(k=0;k<lbox->Nspecies;k+=1)
    for(i=0;i<lbox->Nbox;i+=1)
      lbox->neqbox[k][i]=NULL;
  return 0;
}

int null_box(Box *lbox){/*testar depois*/
  int k,i,nullbox;

  nullbox=0;
  if(lbox->Nbox<=0)
    return 2;
  for(i=0;i<lbox->Nbox;i+=1)
    if(lbox->box[i]!=NULL)
      nullbox=1;

  for(k=0;k<lbox->Nspecies;k+=1)
    for(i=0;i<lbox->Nbox;i+=1)
      if(lbox->neqbox[k][i]!=NULL)
        nullbox=1;
  
  return nullbox;
}

int linear_index(Box *lbox)
{
	long int i,lindex;
	lindex=lbox->k[(lbox->D)-1];
	for(i=lbox->D-2;i>=0;i-=1)
		lindex=lbox->k[i]+lindex*(lbox->Lbox[i]);
	return lindex;
}

int set_box(int N_sph,SPHeq_list *sph_eq,int Nspecies,int *N,SPHneq_list **sph_neq,Box *lbox){
  int i,k,l,err=0;
  long int lindex;
  double xm;
  for(i=0;i<N_sph;i+=1){
    sph_eq[i].bnext=NULL;
    for(l=1;l<=lbox->D;l+=1){
      xm=((sph_eq[i].p).x[l] - lbox->xmin[l-1])/(lbox->dx[l-1]);
      lbox->k[l-1]=(int) xm;
    }
    lindex=linear_index(lbox);
    if(lindex >= lbox->Nbox || lindex < 0){
      printf("particula %d fora das caixas: ",i);
      for(l=0;l<=lbox->D;l+=1)
        printf("%lf ",sph_eq[i].p.x[l]);
      printf("\n");
      return 1;
    }
    err=SPHeq_append_bhead(&(lbox->box[lindex]),&(sph_eq[i]));
  }
	
  for(k=0;k<Nspecies;k+=1)
    for(i=0;i<N[k];i+=1){
      for(l=1;l<=lbox->D;l+=1){
	    xm = (sph_neq[k][i].p.x[l]-lbox->xmin[l-1])/(lbox->dx[l-1]);
	    lbox->k[l-1]=(int) xm;
	  }
      lindex = linear_index(lbox);
      if(lindex >= lbox->Nbox || lindex < 0){
        printf("particula (%d,%d) fora das caixas\n",k,i);
        for(l=0;l<=lbox->D;l+=1)
          printf("%lf ",sph_neq[k][i].p.x[l]);
        printf("\n");
	      return -1;
	    }
	  err=SPHneq_append_bhead(&(lbox->neqbox[k][lindex]),&(sph_neq[k][i]));
    }
	
	return 0;
}

int fetch_inter(double *x, Box *lbox,double kh,SPHeq_list **inter_eq,SPHneq_list ***inter_neq)
{
  int k,l,range,lindex,Nboxes,box,boxaux,isGoodBox,err;
  double xm;
  SPHeq_list *aux_eq;
  SPHneq_list *aux_neq;
  
  xm=kh/minD(lbox->dx,lbox->D);
  range = (int)xm;
  Nboxes=1;
  
  (*inter_eq)=NULL;
  for(k=0;k<lbox->Nspecies;k+=1)
    (*inter_neq)[k]=NULL;
    
  for(l=1;l<=lbox->D;l+=1){
    xm=(x[l] - lbox->xmin[l-1])/(lbox->dx[l-1]);
    lbox->ind[l-1]=(int)xm;
    Nboxes*=2*range+1;
  }
  for(box=0;box<Nboxes;box+=1){
    isGoodBox=0;
    boxaux=box;
    for(l=0;l<lbox->D;l+=1){
      lbox->k[l] = lbox->ind[l] + (boxaux%(2*range+1))-range;
      boxaux = boxaux/(2*range+1);
      if(lbox->k[l]<0 || lbox->k[l]>=lbox->Lbox[l])
	      isGoodBox=1;
	  }
    if(isGoodBox!=0)
      continue;
    lindex=linear_index(lbox);

    aux_eq=lbox->box[lindex];
    while(aux_eq!=NULL){
      err=SPHeq_append_ihead(inter_eq,aux_eq);
      if(err!=0){
		    err=0;  
		    printf("problemas - SPHneq_append_ihead\n");
      }
      aux_eq=aux_eq->bnext;
    }

    for(k=0;k<lbox->Nspecies;k+=1){
      aux_neq=lbox->neqbox[k][lindex];
      while(aux_neq!=NULL){
        err=SPHneq_append_ihead(&((*inter_neq)[k]),aux_neq);/*QUE PORRA DE SINTAXE E ESSA?????????????*/
        if(err!=0){
          err=0;
          printf("problemas - SPHneq_append_ihead\n");
        }
        aux_neq=aux_neq->bnext;
      }
    }
  }

  return 0;
}

int fetch_inter_eq(double *x, Box *lbox,double kh,SPHeq_list **inter_eq)
{
  int l,range,lindex,Nboxes,box,boxaux,isGoodBox,err;
  double xm;
  SPHeq_list *aux_eq;
  
  xm=kh/minD(lbox->dx,lbox->D);
  range = (int)xm;
  Nboxes=1;
  
  (*inter_eq)=NULL;
    
  for(l=1;l<=lbox->D;l+=1){
    xm=(x[l] - lbox->xmin[l-1])/(lbox->dx[l-1]);
    lbox->ind[l-1]=(int)xm;
    Nboxes*=2*range+1;
  }
  for(box=0;box<Nboxes;box+=1){
    isGoodBox=0;
    boxaux=box;
    for(l=0;l<lbox->D;l+=1){
      lbox->k[l] = lbox->ind[l] + (boxaux%(2*range+1))-range;
      boxaux = boxaux/(2*range+1);
      if(lbox->k[l]<0 || lbox->k[l]>=lbox->Lbox[l])
	      isGoodBox=1;
	  }
    if(isGoodBox!=0)
      continue;
    lindex=linear_index(lbox);

    aux_eq=lbox->box[lindex];
    while(aux_eq!=NULL){
      err=SPHeq_append_ihead(inter_eq,aux_eq);
      if(err!=0){
		    err=0;  
		    printf("problemas - SPHneq_append_ihead\n");
      }
      aux_eq=aux_eq->bnext;
    }
  }

  return 0;
}

int fetch_inter_neq(double *x, Box *lbox,double kh,int k,SPHneq_list ***inter_neq)
{
  int l,range,lindex,Nboxes,box,boxaux,isGoodBox,err;
  double xm;
  SPHneq_list *aux_neq;
  
  xm=kh/minD(lbox->dx,lbox->D);
  range = (int)xm;
  Nboxes=1;
  
  (*inter_neq)[k]=NULL;
    
  for(l=1;l<=lbox->D;l+=1){
    xm=(x[l] - lbox->xmin[l-1])/(lbox->dx[l-1]);
    lbox->ind[l-1]=(int)xm;
    Nboxes*=2*range+1;
  }
  for(box=0;box<Nboxes;box+=1){
    isGoodBox=0;
    boxaux=box;
    for(l=0;l<lbox->D;l+=1){
      lbox->k[l] = lbox->ind[l] + (boxaux%(2*range+1))-range;
      boxaux = boxaux/(2*range+1);
      if(lbox->k[l]<0 || lbox->k[l]>=lbox->Lbox[l])
	      isGoodBox=1;
	  }
    if(isGoodBox!=0)
      continue;
    lindex=linear_index(lbox);

    aux_neq=lbox->neqbox[k][lindex];
    while(aux_neq!=NULL){
      err=SPHneq_append_ihead(&((*inter_neq)[k]),aux_neq);/*QUE PORRA DE SINTAXE E ESSA?????????????*/
      if(err!=0){
        err=0;
        printf("problemas - SPHneq_append_ihead\n");
      }
      aux_neq=aux_neq->bnext;
    }
  }

  return 0;
}

/*
 * Esse Setup_SPH foi feita utilizando o artigo do BJP com densidade de 
 * entropia como quantidade conservada.
 * 
 * Utiliza tambem poeira fora de equilibrio
 */
 
int setup_sph(int D,double t,double h,double kh,int N_sph,SPHeq_list *sph_eq,
              int Nspecies, int *N, SPHneq_list **sph_neq,
              Box *lbox,double (*w)(double,double),int (*EoS)(SPHeq_particle*))
{
  int i,k,l,err;
  double x[D+1],dist,us,p_can,sgn;
  SPHeq_list *inter_eq,*aux_eq;
  SPHneq_list **inter_neq,*aux_neq;
  
  inter_neq=(SPHneq_list **)malloc(Nspecies*sizeof(SPHneq_list*));
  
  err=clean_box(lbox);
  if(err!=0)
    printf("erro na clean_box\n");
  err=0;
  err=SPHeq_clean_list(N_sph,sph_eq);
  if(err!=0)
    printf("erro na SPHeq_clean_list\n");
  err=0;
  err=SPHneq_clean_list(Nspecies,N,sph_neq);
  if(err!=0)
    printf("erro na SPHneq_clean_list\n");
  err=set_box(N_sph,sph_eq,Nspecies,N,sph_neq,lbox);
  if(err!=0){
    printf("erro na set_box\n");
    return err;
  }
  err=null_box(lbox);
  if(err!=1)
    printf("caixas vazias-%d\n",err);
  
  for(i=0;i<N_sph;i+=1){
    /* 
     if(sph_eq[i].p.fo!=0)
       continue;
     */
    
    sph_eq[i].p.u[0]=1.0;
    for(l=1;l<=D;l+=1)
      sph_eq[i].p.u[0] += (sph_eq[i].p.u[l])*(sph_eq[i].p.u[l]);
    sph_eq[i].p.u[0]=sqrt(sph_eq[i].p.u[0]);
    for(l=0;l<=D;l+=1)
      sph_eq[i].p.v[l]=(sph_eq[i].p.u[l])/(sph_eq[i].p.u[0]);
    
    for(l=0;l<=D;l+=1)
      x[l]=sph_eq[i].p.x[l];
    err=0;
    err=fetch_inter_eq(x,lbox,kh,&inter_eq);
    if(err!=0)
      printf("erro no fetch_inter_eq\n");
    sph_eq[i].p.rho=0.0;
    while(inter_eq!=NULL){
      dist=0.0;
      for(l=1;l<=D;l+=1)
        dist+=(x[l]-(inter_eq->p).x[l])*(x[l]-(inter_eq->p).x[l]);
      dist=sqrt(dist);
      
      sph_eq[i].p.rho += ((inter_eq->p).ni)*w(dist,h);
      
      aux_eq=inter_eq;
      inter_eq=inter_eq->inext;
      aux_eq->inext=NULL;
    }
    
    sph_eq[i].p.rho_pa = sph_eq[i].p.rho_p; /* Freezout Variable */
    sph_eq[i].p.Ta = sph_eq[i].p.T; /* Freezout Variable */
    
    sph_eq[i].p.rho_p=(sph_eq[i].p.rho)/(sph_eq[i].p.u[0]);
    sph_eq[i].p.s_p = (sph_eq[i].p.rho_p)*(sph_eq[i].p.S);
    sph_eq[i].p.s = (sph_eq[i].p.s_p)*(sph_eq[i].p.u[0]);
    
    sph_eq[i].p.nc = ((sph_eq[i].p.nc)*EtaCharg_qg)/(sph_eq[i].p.u[0]);
    sph_eq[i].p.Nc = (sph_eq[i].p.nc)/(sph_eq[i].p.rho_p);

    err=EoS(&(sph_eq[i].p));
    if(err!=0)
      printf("Erro na Eq. de Estado");
    /*
      if(sph_eq[i].p.T < T_fo)
        sph_eq[i].p.fo=1;
    */
  }
  
  /*Calculando efeito de carga via carga media*/
  /*for(i=0;i<N_sph;i+=1){ 
    for(l=0;l<=D;l+=1)
      x[l]=sph_eq[i].p.x[l];
    err=0;
    err=fetch_inter_eq(x,lbox,kh,&inter_eq);
    if(err!=0)
      printf("erro no fetch_inter_eq\n");
    sph_eq[i].p.nb=0.0;
    while(inter_eq!=NULL){
      dist=0.0;
      for(l=1;l<=D;l+=1)
        dist+=(x[l]-(inter_eq->p).x[l])*(x[l]-(inter_eq->p).x[l]);
      dist=sqrt(dist);
      
      sph_eq[i].p.nb += ((inter_eq->p).ni)*((inter_eq->p).Nc)*w(dist,h);
      
      aux_eq=inter_eq;
      inter_eq=inter_eq->inext;
      aux_eq->inext=NULL;
    }
    sph_eq[i].p.Nb = (sph_eq[i].p.nb)/(sph_eq[i].p.rho);
  }*/
  
  for(k=0;k<Nspecies;k+=1){
    for(i=0;i<N[k];i+=1){
      sph_neq[k][i].p.u[0]=1.0;
      for(l=1;l<=D;l+=1)
        sph_neq[k][i].p.u[0]+=(sph_neq[k][i].p.u[l])*(sph_neq[k][i].p.u[l]);
      sph_neq[k][i].p.u[0]=sqrt(sph_neq[k][i].p.u[0]);
      for(l=0;l<=D;l+=1)
        sph_neq[k][i].p.v[l]=(sph_neq[k][i].p.u[l])/(sph_neq[k][i].p.u[0]);
      
      for(l=0;l<=D;l+=1)
        x[l]=sph_neq[k][i].p.x[l];
      err=0;
      err=fetch_inter_neq(x,lbox,kh,k,&(inter_neq));
      if(err!=0)
        printf("Erro no fetch_inter_neq - %d \n",k);
      sph_neq[k][i].p.rho=0.0;
      while(inter_neq[k]!=NULL){
        dist=0.0;
        for(l=1;l<=D;l+=1)
          dist+=(x[l]-(inter_neq[k]->p).x[l])*(x[l]-(inter_neq[k]->p).x[l]);
        dist=sqrt(dist);
        
        sph_neq[k][i].p.rho += ((inter_neq[k]->p).ni)*w(dist,h);
        
        aux_neq=inter_neq[k];
        inter_neq[k]=inter_neq[k]->inext;
        aux_neq->inext=NULL;
      }
      sph_neq[k][i].p.rho_p=(sph_neq[k][i].p.rho)/(sph_neq[k][i].p.u[0]);
      sph_neq[k][i].p.n = (sph_neq[k][i].p.rho)*(sph_neq[k][i].p.N_p);
      sph_neq[k][i].p.n_p = (sph_neq[k][i].p.rho_p)*(sph_neq[k][i].p.N_p);
    }
  }
  
  free(inter_neq);
  return 0;
}

/*
 * 2+1 D Boost Invariant Setup SPH
 */

int ssph_2p1bi(int D,double t,double h,double kh,int N_sph,SPHeq_list *sph_eq,
               int Nspecies, int *N, SPHneq_list **sph_neq,
               Box *lbox,double (*w)(double,double),int (*EoS)(SPHeq_particle*))
{
  int i,k,l,err;
  double x[D+1],dist,us,p_can,sgn;
  SPHeq_list *inter_eq,*aux_eq;
  SPHneq_list **inter_neq,*aux_neq;

  if(D!=2){
    printf("inside ssph_2p1bi D=%d\n",D);return (-10);}  
  
  inter_neq=(SPHneq_list **)malloc(Nspecies*sizeof(SPHneq_list*));
  
  err=clean_box(lbox); if(err!=0){ printf("erro na clean_box\n"); return err;}
  err=SPHeq_clean_list(N_sph,sph_eq); if(err!=0){ printf("erro na SPHeq_clean_list\n"); return err;}
  err=SPHneq_clean_list(Nspecies,N,sph_neq); if(err!=0){printf("erro na SPHneq_clean_list\n");return err;}
  err=set_box(N_sph,sph_eq,Nspecies,N,sph_neq,lbox); if(err!=0){ printf("erro na set_box\n");return err;}
  err=null_box(lbox); if(err!=1){ printf("caixas vazias-%d\n",err);}
  
  for(i=0;i<N_sph;i+=1){
    /* 
     if(sph_eq[i].p.fo!=0)
       continue;
     */
    
    sph_eq[i].p.u[0]=1.0;
    for(l=1;l<=D;l+=1)
      sph_eq[i].p.u[0] += (sph_eq[i].p.u[l])*(sph_eq[i].p.u[l]);
    sph_eq[i].p.u[0]=sqrt(sph_eq[i].p.u[0]);
    for(l=0;l<=D;l+=1)
      sph_eq[i].p.v[l]=(sph_eq[i].p.u[l])/(sph_eq[i].p.u[0]);
    
    for(l=0;l<=D;l+=1)
      x[l]=sph_eq[i].p.x[l];
    err=0;
    err=fetch_inter_eq(x,lbox,kh,&inter_eq);
    if(err!=0)
      printf("erro no fetch_inter_eq\n");
    sph_eq[i].p.rho=0.0;
    /*sph_eq[i].p.s = 0.0;*/
    while(inter_eq!=NULL){
      dist=0.0;
      for(l=1;l<=D;l+=1)
        dist+=(x[l]-(inter_eq->p).x[l])*(x[l]-(inter_eq->p).x[l]);
      dist=sqrt(dist);
      
      sph_eq[i].p.rho += ((inter_eq->p).ni)*w(dist,h);
      /*sph_eq[i].p.s   += ((inter_eq->p).ni)*((inter_eq->p).S)*w(dist,h);*/
      
      aux_eq=inter_eq;
      inter_eq=inter_eq->inext;
      aux_eq->inext=NULL;
    }
    
    /*sph_eq[i].p.rho_pa = sph_eq[i].p.rho_p;  /*Freezout Variable */
    /*sph_eq[i].p.Ta = sph_eq[i].p.T;         /* Freezout Variable */
    
    sph_eq[i].p.rho_p=(sph_eq[i].p.rho)/((sph_eq[i].p.u[0])*t); /* scaling com tempo */
    /*sph_eq[i].p.s_p = (sph_eq[i].p.s)/((sph_eq[i].p.u[0])*t);*/
    sph_eq[i].p.s_p = (sph_eq[i].p.rho_p)*(sph_eq[i].p.S);
    sph_eq[i].p.s = (sph_eq[i].p.rho)*(sph_eq[i].p.S);

    sph_eq[i].p.nc = 0.0;
    sph_eq[i].p.Nc = 0.0;

    /*
     * sph_eq[i].p.s_p = (sph_eq[i].p.rho_p)*(sph_eq[i].p.S);
       sph_eq[i].p.s = (sph_eq[i].p.rho)*(sph_eq[i].p.S);  */

    err=EoS(&(sph_eq[i].p));
    if(err!=0)
      printf("Erro na Eq. de Estado");
    /*
      if(sph_eq[i].p.T < T_fo)
        sph_eq[i].p.fo=1;
    */
  }
  
  free(inter_neq);
  return 0;
}
