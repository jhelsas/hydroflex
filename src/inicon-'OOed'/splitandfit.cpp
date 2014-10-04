#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "splitandfit.h"
#include <string.h>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

using namespace std;

#define _USE_MATH_DEFINES

int Domain :: D=0;
vector <Domain> Domain :: dms;

/*
 * Arquivo de implementação das rotinas de divisão e fit
 * 
 * Leia splitandfit.h para quais funções usar e quais funções 
 * são internas da biblioteca
 * 
 */

/*
 * Public Functions
 * 
 * domain():
 *   Default Constructor
 * 
 * domain(int,int,int):
 *   Constructor
 * 
 * init(int,int,int):
 *  Initializer Fuction
 *  Usage: Starting From the Default Constructor,
 *    Input the dimension, Number of vertex and type of domain 
 */

Domain :: Domain(): xv(0){
  Nv=0;good=0;type=-1;
  S=0.;
}

Domain :: Domain(int dim, int Nvertex, int Type): xv(dim*Nvertex){
  Domain :: D=dim; Nv=Nvertex; good=0; type=Type;
  S=0.;
}

Domain :: Domain (int dim,int Ntri,int Type,double *Xv): xv(dim*Ntri*(D+1)){
  int err;
  Domain :: D=dim; Nv=(D+1)*Ntri; good=0; type=Type;
  S=0.;
  err=init_triangle((Domain :: D),Ntri,Xv);
  if(err!=0)
    cout << "Problemas construtor de triangulos\n";
}
 
Domain :: Domain (int dim,int Nvertex,int Type,double *xi,double *xf): xv(dim*Nvertex){
  int err;
  Domain :: D=dim; Nv=Nvertex; good=0; type=Type;
  S=0.;
  err=init_cube((Domain :: D),xi,xf);
  if(err!=0)
    cout << "Problemas construtor de cubo\n";
}

int Domain ::init(int dim, int Nvertex, int Type){
  D=dim; Nv=Nvertex; good=0; type=Type;
  xv.resize(dim*Nvertex);
  S=0.;
  return 0;
}

Domain ::~Domain(){}

/*
 * Public Functions
 * 
 * domain_count: 
 *   Output the number of domains on dom
 * 
 * domain_check:
 *   check if all the domains are good, or if there are spurious domains
 * 
 * el_print:
 *   Print one domain coordinates 
 * 
 * print_domain_n:
 *   Debug function - print the n-th element of dom's vertex coordinates
 * 
 * print_domain:
 *   print the coordinates of all vertexes of the domain list
 * 
 * check_inside:
 *   check if a point is inside a given domain
 *   to use with functions as input for the code
 */
 
int Domain ::domain_count(){
  int i=0;
  vector <Domain> :: iterator el;
  for(el=(Domain :: dms).begin();el!=(Domain :: dms).end();el++)
    if(el->good==0)
      i+=1;
  return i;
}

int Domain ::domain_check(){
  vector <Domain> :: iterator el;
  for(el=(Domain :: dms).begin();el!=(Domain :: dms).end();el++)
    if(el->good!=0)
      return 1;  
  return 0;
}

int Domain ::el_print(){
  int j,l;
  for(j=0;j<Nv;j+=1){
    cout << "( ";
    for(l=0;l<Domain :: D;l+=1)
      cout << xv[j*(Domain :: D)+l] << " ";
    cout << ") ";
  }
  cout << endl;
  return 0;
}

int Domain ::print_domain_n(int n){
  int i,l;
 
  cout << "domain " << n << ":\n";
  for(i=0;i<(Domain :: dms)[n].Nv;i+=1){
    cout << "  vertex " << i << ": ";
    for(l=0;l<(Domain :: D);l+=1)
      cout << (Domain :: dms)[n].xv[i*(Domain :: D)+l] << " ";
    cout << "\n";
  }
  return 0;
}

int Domain ::print_domain(){
  int n,i,l;
  for(n=0;n<(Domain :: dms).size();n+=1){
    cout << "domain " << n << ":\n";
    for(i=0;i<(Domain :: dms)[n].Nv;i+=1){
      cout << "  vertex " << i << ": ";
      for(l=0;l<(Domain :: D);l+=1)
        cout << (Domain :: dms)[n].xv[i*(Domain :: D)+l] << " ";
      cout << "\n";
    }
    cout << "aspect ratio: " << (Domain :: dms)[n].aspect_ratio() << "\n";
  }
  cout << "\n";
  return 0;
}

// voltar aqui
int Domain ::check_inside(double *x){ 
  int i,l,err;
  double lmb[(Domain :: D)],min=100000.,max=-1000000.;
  
  if(type == 0){
    return 0;
  }
  else if(type == 1){
    err=bc_coord(x,lmb);if(err!=0){printf("problemas gsl_wrapper\n");return 0;}
    err=bc_check((Domain :: D),lmb);
    if(err==0)
      return 0;
    else
      return 1;
  }
  else
    return 1;
  return 1;
}

/*
 * Public Function
 * 
 * init_cube:
 *   Initalize the domain with a (d-)cube domain aligned with 
 *   artesian axis, the lower and upper limits of the cube are given by
 *   xl[] and xu[] 
 * 
 *   e.g.:
 *   2-cube = square (regular) or paralelogram (irregular)
 *   3-cube = cube (regular) or paralelepiped (irregular)
 * 
 * init_triangle:
 *   
 *   Initalize the domain with a (d-)tiangles/simplex form the number of
 *   triangles and from the list of the triangles positions
 * 
 *   e.g.:
 *   2-triangle/simplex = triangles (regular or irregular)
 *   3-triangle/simples = tetrahedron (regular or irregular)
 */

int Domain ::init_cube(int dim,double xl[],double xu[]){
  int i,l,err,ind[dim];
  Domain :: D = dim;
  if((Domain :: dms).size()>0)
    return 1;
  err=this->init((Domain :: D),1<<(Domain :: D),0); if(err!=0) return err;
  for(i=0;i<Nv;i+=1){
    err=rollout(ind,&i,(Domain :: D),2);
    for(l=0;l<(Domain :: D);l+=1){
      if(ind[l]==0)
        xv[i*(Domain :: D)+l]=xl[l]; 
      else
        xv[i*(Domain :: D)+l]=xu[l];
    }
  }
  (Domain :: dms).push_back(*this);
  return 0;
}

int Domain ::init_triangle(int dim,double Ntri,double *xv){
  int i,n,l,err;
  Domain :: D=dim;
  if((Domain :: dms).size()>0)
    return 1;
  
  for(n=0;n<Ntri;n+=1){
    Domain mdel;
    err=mdel.init((Domain :: D),(Domain :: D)+1,1); if(err!=0) return err;
    for(i=0;i<mdel.Nv;i+=1)
      for(l=0;l<(Domain :: D);l+=1)
        mdel.xv[i*(Domain :: D)+l] = xv[n*(Domain :: D)*((Domain :: D)+1)+i*(Domain :: D)+l];
    (Domain :: dms).push_back(mdel);
  }

  return 0;
}

/*
 * Public Function
 * 
 * domain_split:
 *   Main function of the code
 *   
 *   Takes a domain list dom, and splits til the weight within it 
 *   non-breaked domain is less than a given cutoff
 * 
 *   The weight is calculated as the integral of the function F
 *   on the associated domain. The breaking is done automatically, and the
 *   function returns either sucess or fail.
 * 
 *   The result is stored on the input domain list dom. It stills contains
 *   spurious (breaked), thus, it's necessary to call clean_domain before
 *   turning it into SPH particles
 */
 
int Domain ::domain_split(double cutoff,gsl_monte_function F){
  unsigned int id;
  int i,err,l;
  double xl[(Domain :: D)],xu[(Domain :: D)],S,erd;
  size_t tcalls=3000,calls = 500000;
  const gsl_rng_type *T;
  wparams *lpar;
  gsl_rng *r;
  gsl_monte_miser_state *s_m=gsl_monte_miser_alloc (D);
    
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  lpar=(wparams*)(F.params);
    
  id=0;
  while(id<(Domain :: dms).size()){
	
    lpar->mdel=&((Domain :: dms)[id]);
    
    S=0.0;
    if((Domain :: dms)[id].type==0){
      for(l=0;l<(Domain :: D);l+=1){
        xl[l]=(Domain :: dms)[id].xv[(Domain :: D)*0+l];
        xu[l]=(Domain :: dms)[id].xv[(Domain :: D)*((1<<(Domain :: D))-1)+l]; 
      }
    }
    else{
      for(l=0;l<(Domain :: D);l+=1){
        xl[l]=100000.;
        xu[l]=-10000.;
        for(i=0;i<((Domain :: dms)[id].Nv);i+=1){
          if((Domain :: dms)[id].xv[(Domain :: D)*i+l] < xl[l])
            xl[l] = (Domain :: dms)[id].xv[(Domain :: D)*i+l];
          if((Domain :: dms)[id].xv[(Domain :: D)*i+l] > xu[l])
            xu[l] = (Domain :: dms)[id].xv[(Domain :: D)*i+l];
        }
      }
    }
    
    gsl_monte_miser_integrate(&F,xl,xu,(Domain :: D),calls,r,s_m,&S,&erd);
    (Domain :: dms)[id].S=S;
    
    if(S > cutoff){
      (Domain :: dms)[id].good=1;
      if((Domain :: dms)[id].type==0){
        err=(Domain :: dms)[id].cubic_split(id);if(err!=0){return err;}}
      else if((Domain :: dms)[id].type==1){
        if(D==2){
          err=(Domain :: dms)[id].triangle_midpoint_split(id);if(err!=0){ cout << "in: " << err << endl;return err;}}
        else{
          err=(Domain :: dms)[id].bc_simplex_split(id);if(err!=0) return err;}
      }
    }
    else{
      gsl_monte_miser_integrate(&F,xl,xu,(Domain :: D),calls,r,s_m,&S,&erd);
      (Domain :: dms)[id].S=S;
    }
    
    id++;
  }
  
  gsl_monte_miser_free (s_m);
  gsl_rng_free (r);
  
  return 0;
}

/*
 * Public Functions
 * 
 * clean_domain:
 *   Use after domain_split
 *   Clean Spurious (already split) Domains
 * 
 * print_sph:
 *   Use after clean_domain
 *   Converts the domains weights and positions to produce 
 *   SPH particles
 */

int Domain ::clean_domain(){
  unsigned int i;
  vector <Domain> :: iterator el;  
  for(i=0;i<(Domain :: dms).size();i+=1)
    if((Domain :: dms)[i].good!=0){
      (Domain :: dms).erase((Domain :: dms).begin()+i);
      i-=1;}
  if(domain_check()!=0)
    cout<< "algum dominio esta ruim\n";
  
  return 0;
}

int Domain ::print_sph(const char *filename){
  int j,l,it=0,N;
  ofstream sphfile;
  vector <Domain> :: iterator el;
  if((Domain :: D)<=0 || filename==NULL || (Domain :: dms).size()==0)
    return 1;
  N=(Domain :: dms).size();
  sphfile.open(filename);
  sphfile << (Domain :: D) << "  " << N << endl;
  for(el=(Domain :: dms).begin();el!=(Domain :: dms).end();el++){
    if(el->type==0){
      for(l=0;l<(Domain :: D);l+=1) 
        sphfile << (el->xv[(Domain :: D)*0+l]+el->xv[(Domain :: D)*((1<<(Domain :: D))-1)+l])/2. << " ";
      for(l=0;l<D;l+=1) 
        sphfile << 0. << " ";
    }
    else{
      double x[l];
      for(l=0;l<(Domain :: D);l+=1){
        x[l]=0.;
        for(j=0;j<el->Nv;j+=1)
          x[l]+= el->xv[(Domain :: D)*j+l];
        x[l]=x[l]/(el->Nv);
        sphfile << x[l] << " ";
      }
      for(l=0;l<(Domain :: D);l+=1)
        sphfile << 0. << " ";
    }
    sphfile << el->S <<"\n";
    it++;
  }
  sphfile.close();
  return 0;
}

/* 
 * Internal Function
 * 
 * rollin , rollout , aspect_ratio
 */

int rollin(int *ind,int *indx,int D,int n){
  int i,index=0,tmp=1;
  if(D<=0 || n<=0)
    return 1;
  for(i=0;i<D;i+=1){
    index+=tmp*ind[i];tmp*=n;}
  *indx=index;
  return 0;
}

int rollout(int *ind,int *indx,int D,int n){
  int i,index=*indx;
  if(D<=0 || n<=0)
    return 1;
  for(i=0;i<D;i+=1){
    ind[i] = index%n;index=index/n;}
  return 0;
}

double Domain ::aspect_ratio(){
  int i,j,l;
  double dist,min=100000000.,max=-1000000000.;
  if(type!=1)
    return -1.;
  for(i=0;i<Nv;i+=1){
    for(j=0;j<Nv;j+=1){
      if(i==j)
        continue;
      dist=0.;
      for(l=0;l<(Domain :: D);l+=1)
        dist+= (xv[i*(Domain :: D)+l]-xv[j*(Domain :: D)+l])*(xv[i*(Domain :: D)+l]-xv[j*(Domain :: D)+l]);
      dist=sqrt(dist);
      if(dist>max)
        max=dist;
      if(dist<min)
        min=dist;
    }
  }
  return max/min;
}

/*
 * Internal Function
 * 
 * Spliting routine for (d-dimensional) cubic domains
 * Split 1 (d-)cube into 2^d smaller (d-)cubes
 * 
 * e.g.:
 * 2-cube = square (regular) or paralelogram (irregular)
 * 3-cube = cube (regular) or paralelepiped (irregular)
 */

int Domain ::cubic_split(int n){
  int i,j,k,l,err,indx,indi[D],ind[D],inda[D];
  double xv[(1<<(Domain :: D))*(Domain :: D)];
  
  if((Domain :: dms)[n].type!=0)
    return -1;
  
  for(i=0;i<(1<<(Domain :: D));i+=1)
    for(l=0;l<(Domain :: D);l+=1)
      xv[i*(Domain :: D)+l] = (Domain :: dms)[n].xv[i*(Domain :: D)+l];
    
  for(i=0;i<(1<<(Domain :: D));i+=1){
    Domain tdel((Domain :: D),1<<(Domain :: D),0); // eu ainda vou ter de entender pq isso está aí
    tdel.good=0;
    err=rollout(ind,&i,(Domain :: D),2); if(err!=0) return err;
    
    for(l=0;l<(Domain :: D);l+=1){
      
      for(k=0;k<(Domain :: D);k+=1)
        inda[k] = ind[k];
      inda[l]=1-ind[l];
      
      err=rollin(inda,&indx,(Domain :: D),2); if(err!=0) return err;
      
      for(j=0;j<(Domain :: dms)[n].Nv;j+=1){
        err=rollout(indi,&j,(Domain :: D),2); if(err!=0) return err;      
        if(indi[l]==ind[l])
          tdel.xv[(Domain :: D)*j+l]=xv[(Domain :: D)*i+l];
        else
          tdel.xv[(Domain :: D)*j+l]=(xv[(Domain :: D)*i+l]+xv[(Domain :: D)*indx+l])/2.;
      }
    }
    tdel.S=0.0;
    (Domain :: dms).push_back(tdel);
  }
  
  return 0;
}

/*
 * Internal function
 * 
 * Spliting routine for (d-dimensional) simplex into d+1 simplexes 
 * using the baricenter as a vertex. 
 * Note: Not optimal split
 * 
 * e.g.:
 * 2-triangle/simplex = triangles (regular or irregular)
 * 3-triangle/simples = tetrahedron (regular or irregular)
 */

int Domain ::bc_simplex_split(int n){
  int i,j,l,ind; 
  double bc[(Domain :: D)];
  
  if((Domain :: dms)[n].type!=1)
    return -1;
  for(l=0;l<(Domain :: D);l+=1){
    bc[l]=0.;
    for(i=0;i<(Domain :: dms)[n].Nv;i+=1)
      bc[l]+=(Domain :: dms)[n].xv[(Domain :: D)*i+l];
    bc[l]*= 1./((double)(Domain :: dms)[n].Nv);
  }
  
  for(i=0;i<(Domain :: dms)[n].Nv;i+=1){
    Domain tdel((Domain :: D),(Domain :: D)+1,1);
    tdel.good=0;
    for(j=0;j<=(Domain :: D);j+=1){
      for(l=0;l<(Domain :: D);l+=1){
        if(j==i)
          tdel.xv[(Domain :: D)*j+l]=bc[l];
        else
          tdel.xv[(Domain :: D)*j+l]=(Domain :: dms)[n].xv[(Domain :: D)*j+l];
      }
    }
    tdel.S=0.0;
    (Domain :: dms).push_back(tdel);
  }
  return 0;
}

/*
 * Internal functions
 * 
 * bc_coord:
 * Calculate the baricentric coordinates *lmb of a given point *r 
 * to a given (d-) triangular domain
 *  
 * bc_check
 * Check if a point is inside in a given (d-)triangular domain
 * using it's baricentric coordinates.
 * 
 * bc_test:
 * print the baricentric coordinates of *r and if they are inside
 * the domain or not
 */

int Domain ::bc_coord(double *r,double *lmb){
  int l,n,s;
  double T[(Domain :: D)*(Domain :: D)],B[(Domain :: D)];
  if(type!=1)
    return -1;
  
  for(l=0;l<D;l+=1){
    B[l] = r[l] - xv[(Domain :: D)*(Domain :: D)+l];
    for(n=0;n<(Domain :: D);n+=1)
      T[l*(Domain :: D)+n] = xv[n*(Domain :: D)+l] - xv[(Domain :: D)*(Domain :: D)+l];
  }
  
  gsl_matrix_view m = gsl_matrix_view_array (T, (Domain :: D), (Domain :: D));
  gsl_vector_view b = gsl_vector_view_array (B, (Domain :: D));
  gsl_vector *x = gsl_vector_alloc ((Domain :: D));  
  gsl_permutation * p = gsl_permutation_alloc ((Domain :: D));
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
  for(l=0;l<(Domain :: D);l+=1)
    lmb[l] = gsl_vector_get(x,l);
  gsl_permutation_free (p);
  gsl_vector_free (x);
  
  return 0;
}
 
int bc_check(int D,double *lmb){
  int l;
  double Lmb;
  if(D<=0)
    return -1;
  Lmb=0.;
  for(l=0;l< D;l+=1){
    Lmb+=lmb[l];
    if(lmb[l]<0 || lmb[l]>1)
      return 1;
  }
  if(Lmb <0 || Lmb > 1)
    return 1;
  return 0;
}

int Domain ::bc_test(double *r){
  int err,l;
  double lmb[(Domain :: D)],Lmb;
  err=bc_coord(r,lmb);if(err!=0) return err;
  cout << "bc coordinates: " ;
  Lmb=0.;
  for(l=0;l<(Domain :: D);l+=1){
    cout << lmb[l] << " " ;Lmb+= lmb[l];}
  cout << Lmb << "\n ";
    
  err= bc_check((Domain :: D),lmb);
  if(err==0)
    cout << "ponto no interior\n";
  else
    cout << "ponto no exterior\n";
  return 0;
}

/*
 * Internal function
 * 
 * Spliting routine for (2-dimensional) triangles into 4 triangles 
 * midpoints of each side as vertexes.
 * 
 */
 
int Domain ::triangle_midpoint_split(int n){
  int i,l,ind1,ind2;
  double midp[((Domain :: D)+1)*(Domain :: D)];
  
  if((Domain :: dms)[n].type!=1 || (Domain :: D)!=2)
    return -1;
  
  for(i=0;i<(Domain :: dms)[n].Nv;i+=1){
    ind1=(i+1)%((Domain :: dms)[n].Nv);
    ind2=(i+2)%((Domain :: dms)[n].Nv);
    for(l=0;l<(Domain :: D);l+=1)
      midp[i*(Domain :: D)+l] = ((Domain :: dms)[n].xv[ind1*D+l]+(Domain :: dms)[n].xv[ind2*D+l])/2.;
  }
  
  for(i=0;i<(Domain :: dms)[n].Nv;i+=1){
    Domain tdel((Domain :: D),(Domain :: D)+1,1);
    tdel.good=0;
    
    ind1=(i+1)%((Domain :: dms)[n].Nv);
    ind2=(i+2)%((Domain :: dms)[n].Nv);
    
    for(l=0;l<(Domain :: D);l+=1){
      tdel.xv[(Domain :: D)*i+l]=(Domain :: dms)[n].xv[(Domain :: D)*i+l];
      tdel.xv[(Domain :: D)*ind1+l]=midp[ind1*(Domain :: D)+l];
      tdel.xv[(Domain :: D)*ind2+l]=midp[ind2*(Domain :: D)+l];
    }
    tdel.S=0.0;
    (Domain :: dms).push_back(tdel);
  }
  
  {
    Domain tdel((Domain :: D),(Domain :: D)+1,1);
    tdel.good=0;
    
    for(i=0;i<(Domain :: dms)[n].Nv;i+=1)
      for(l=0;l<(Domain :: D);l+=1)
        tdel.xv[(Domain :: D)*i+l] = midp[i*D+l];
        
    tdel.S=0.0;
    (Domain :: dms).push_back(tdel);
  }
  
  return 0;
}

/*
 * Internal function
 * 
 * Spliting routine for (3-dimensional) tetrahedrons into 8 triangles 
 * midpoints of each side as vertexes.
 * 
 */
 
/*
int tetrahedron_midpoint_split(int D,vector <domain> & dom,int n){
  int i,j,l,ind[n];
  double midp[(D+1)*D*D];
  
  if(dom[n].type!=1 || D!=3)
    return -1;
    
  for(i=0;i<dom[n].Nv;i+=1){
    for(j=i;j<dom[n].Nv;j+=1)
      ind[i]=(i+j)%(dom[n].Nv);
    
    for(j=i;j<dom[n].){
      for(l=0;l<D;l+=1)
        midp[((D+1)*i+j)*D+l] = 
    }
  }
  
  return 0;
}*/


/*
 * Public Functions
 * 
 * unit_hexagon and unit2_hexagon
 *   Utility Functions. Create 6 triangles that form an hexagon
 *   to use as input for init_triangle with Ntri=6
 * 
 *   unit_hexagon creates an hexagon with side = 1.
 *   unit2_hexagon creates an hexagon with side = 2.
 */

int unit_hexagon(int Ntri,int D,double *xv){
  if(Ntri!=6 || D!=2)
    return -1;
    
  xv[0*(D+1)*D+0*D+0]=1.; xv[0*(D+1)*D+0*D+1]=0.;
  xv[0*(D+1)*D+1*D+0]=.5; xv[0*(D+1)*D+1*D+1]=sqrt(3.)/2.;
  xv[0*(D+1)*D+2*D+0]=0.; xv[0*(D+1)*D+2*D+1]=0.;
  
  xv[1*(D+1)*D+0*D+0]=.5; xv[1*(D+1)*D+0*D+1]=sqrt(3.)/2.;
  xv[1*(D+1)*D+1*D+0]=-.5;xv[1*(D+1)*D+1*D+1]=sqrt(3.)/2.;
  xv[1*(D+1)*D+2*D+0]=0.; xv[1*(D+1)*D+2*D+1]=0.;
  
  xv[2*(D+1)*D+0*D+0]=0.; xv[2*(D+1)*D+0*D+1]=0.;
  xv[2*(D+1)*D+1*D+0]=-.5;xv[2*(D+1)*D+1*D+1]=sqrt(3.)/2.;
  xv[2*(D+1)*D+2*D+0]=-1.;xv[2*(D+1)*D+2*D+1]=0.;
  
  xv[3*(D+1)*D+0*D+0]=-1.;xv[3*(D+1)*D+0*D+1]=0.;
  xv[3*(D+1)*D+1*D+0]=-.5;xv[3*(D+1)*D+1*D+1]=-sqrt(3.)/2.;
  xv[3*(D+1)*D+2*D+0]=0.; xv[3*(D+1)*D+2*D+1]=0.;
  
  xv[4*(D+1)*D+0*D+0]=-.5;xv[4*(D+1)*D+0*D+1]=-sqrt(3.)/2.;
  xv[4*(D+1)*D+1*D+0]=.5; xv[4*(D+1)*D+1*D+1]=-sqrt(3.)/2.;
  xv[4*(D+1)*D+2*D+0]=0.; xv[4*(D+1)*D+2*D+1]=0.;
  
  xv[5*(D+1)*D+0*D+0]=.5;xv[5*(D+1)*D+0*D+1]=-sqrt(3.)/2.;
  xv[5*(D+1)*D+1*D+0]=1.;xv[5*(D+1)*D+1*D+1]=0.;
  xv[5*(D+1)*D+2*D+0]=0.;xv[5*(D+1)*D+2*D+1]=0.;
  
  return 0;
}
 
int unit2_hexagon(int Ntri,int D,double *xv){
  if(Ntri!=6 || D!=2)
    return -1;
    
  xv[0*(D+1)*D+0*D+0]=2.; xv[0*(D+1)*D+0*D+1]=0.;
  xv[0*(D+1)*D+1*D+0]=1.; xv[0*(D+1)*D+1*D+1]=sqrt(3.);
  xv[0*(D+1)*D+2*D+0]=0.; xv[0*(D+1)*D+2*D+1]=0.;
  
  xv[1*(D+1)*D+0*D+0]=1.; xv[1*(D+1)*D+0*D+1]=sqrt(3.);
  xv[1*(D+1)*D+1*D+0]=-1.;xv[1*(D+1)*D+1*D+1]=sqrt(3.);
  xv[1*(D+1)*D+2*D+0]=0.; xv[1*(D+1)*D+2*D+1]=0.;
  
  xv[2*(D+1)*D+0*D+0]=0.; xv[2*(D+1)*D+0*D+1]=0.;
  xv[2*(D+1)*D+1*D+0]=-1.;xv[2*(D+1)*D+1*D+1]=sqrt(3.);
  xv[2*(D+1)*D+2*D+0]=-2.;xv[2*(D+1)*D+2*D+1]=0.;
  
  xv[3*(D+1)*D+0*D+0]=-2.;xv[3*(D+1)*D+0*D+1]=0.;
  xv[3*(D+1)*D+1*D+0]=-1.;xv[3*(D+1)*D+1*D+1]=-sqrt(3.);
  xv[3*(D+1)*D+2*D+0]=0.; xv[3*(D+1)*D+2*D+1]=0.;
  
  xv[4*(D+1)*D+0*D+0]=-1.;xv[4*(D+1)*D+0*D+1]=-sqrt(3.);
  xv[4*(D+1)*D+1*D+0]=1.; xv[4*(D+1)*D+1*D+1]=-sqrt(3.);
  xv[4*(D+1)*D+2*D+0]=0.; xv[4*(D+1)*D+2*D+1]=0.;
  
  xv[5*(D+1)*D+0*D+0]=1.;xv[5*(D+1)*D+0*D+1]=-sqrt(3.);
  xv[5*(D+1)*D+1*D+0]=2.;xv[5*(D+1)*D+1*D+1]=0.;
  xv[5*(D+1)*D+2*D+0]=0.;xv[5*(D+1)*D+2*D+1]=0.;
  
  return 0;
}


int create_grid(int D,double **xpo,double *xl,double *xu,
                double *dx,int *Np){
  double *xp;
  if(D<=0)
    return -1;
  else if(D==1){
    int i,Nx;
    Nx=(int)((xu[0]-xl[0])/dx[0]);
    xp=new double [Nx*D];
    for(i=0;i<Nx;i+=1)
        xp[i*D+0] = xl[0]+((double)i)*(dx[0]);
    *Np=Nx;
  }
  else if(D==2){
    int k,i,j,Nx,Ny;
    Nx=(int)((xu[0]-xl[0])/dx[0]);
    Ny=(int)((xu[1]-xl[1])/dx[1]);
    xp=new double [Nx*Ny*D];
        
    for(j=0;j<Ny;j+=1)
      for(i=0;i<Nx;i+=1){
        xp[(j*Nx+i)*D+0] = xl[0]+((double)i)*(dx[0]);
        xp[(j*Nx+i)*D+1] = xl[1]+((double)j)*(dx[1]);
      }
    
    *Np=Nx*Ny;
  }
  else if(D==3){
    int i,j,k,Nx,Ny,Nz;
    Nx=(int)((xu[0]-xl[0])/dx[0]);
    Ny=(int)((xu[1]-xl[1])/dx[1]);
    Nz=(int)((xu[2]-xl[2])/dx[2]);
    xp=new double [Nx*Ny*Nz*D];
    for(k=0;k<Nz;k+=1)
      for(j=0;j<Ny;j+=1)
        for(i=0;i<Nx;i+=1){
          xp[((k*Ny+j)*Nx+i)*D+0] = xl[0]+((double)i)*(dx[0]);
          xp[((k*Ny+j)*Nx+i)*D+1] = xl[1]+((double)j)*(dx[1]);
          xp[((k*Ny+j)*Nx+i)*D+2] = xl[2]+((double)k)*(dx[2]);
        }
      
    *Np=Nx*Ny*Nz;
  }
  *xpo=xp;
  return 0;
}

int sph_read(char* filename,int *Dout,int *Nout,double **xout,double **uout,double **Sout){
  int i,l,D,N;
  double *x,*u,*S;
  ifstream sphfile;
  
  sphfile.open(filename);
  sphfile >> D >> N;
  *Dout=D;
  *Nout=N;
  x=new double [N*D];
  u=new double [N*D];
  S=new double [N];
  for(i=0;i<N;i+=1){
    for(l=0;l<D;l+=1)
      sphfile >> x[i*D+l];
    for(l=0;l<D;l+=1)
      sphfile >> u[i*D+l];
    sphfile >> S[i];
  }
  sphfile.close();
  *xout=x;
  *uout=u;
  *Sout=S;
  
  return 0;
}

#define DIM 2

//#define A_d 1.  
 #define A_d (15.)/(7.*M_PI) 
// #define A_d (3.0)/(2.0*MYPI) 

double w_bspline(double r,double h)
{
	double R;
	if(h<=0.) exit(10);
	R=fabs(r)/h;
	if(R>=2.)
		return 0;
	else if((1.<=R)&&(R<2.))
		return (A_d)*(1./6.)*(2.-R)*(2.-R)*(2.-R)/pow(h,DIM);
	else
		return ((A_d)*((2./3.)-(R*R) + (R*R*R/2.0)))/pow(h,DIM);
}

int sph_dens(int D,int N,int Npoints,double *xp,double *x,
             double *S,double h,double xl[],double xu[],
             double (*tf)(double*,size_t,void*),char* filename,void *p){
  int k,i,l,err;
  double a,s,dist;
  Domain *mdel = new Domain();
  wparams par;void *downp;
  ofstream plotfile;
  if(D<=0)
    return -1;
  
  par.p=p;  
  err=mdel->init_cube(D,xl,xu);if(err!=0) return err; 
  par.mdel=mdel;
  
  plotfile.open(filename);
  for(k=0;k<Npoints;k+=1){
    s=0.;
    for(i=0;i<N;i+=1){
      dist=0.;
      for(l=0;l<D;l+=1)
        dist+=(xp[k*D+l]-x[i*D+l])*(xp[k*D+l]-x[i*D+l]);
      dist=sqrt(dist);
      
      s+=S[i]*w_bspline(dist,h);
    }
    a=tf(xp+k*D,D,&par);
    for(l=0;l<D;l+=1)
      plotfile << xp[k*D+l] << " ";
    plotfile << s << " " << a << "\n";
  }
  plotfile.close();
  
  return 0;
}

