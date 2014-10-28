typedef struct SPHeq_particle
{
	int D,id,*ind,fo;
	double ni,rho,q,S,s,s_p,sigma,Nb,Nc;
	double rho_p,e_p,p_p,h_p,hsh_p,T,nb,nc;
	double *x,*u,*v;
  double Ta,rho_pa,Sa,*xa,*ua,*dudt;
} SPHeq_particle;

typedef struct SPHeq_list
{
	SPHeq_particle p;
	struct SPHeq_list *bnext, *inext;
}SPHeq_list;

typedef struct SPHneq_particle
{
	int D,id,*ind;
	double m,ni,q,R;
	double rho,rho_p,n,n_p,N_p;
	double *x,*v,*u;
} SPHneq_particle;

typedef struct SPHneq_list
{
	SPHneq_particle p;
	struct SPHneq_list *next,*bnext, *inext;
}SPHneq_list;

typedef struct Box
{
	int Nbox,Nspecies,D,*ind,*k,*Lbox;
	double *xmin,*xmax,*dx;
	SPHeq_list **box;
	SPHneq_list ***neqbox;
} Box;
