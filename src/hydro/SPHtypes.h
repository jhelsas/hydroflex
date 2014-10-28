/*
 * Description of SPHeq_particle members:
 * 
 * D :    number of spatial dimentions
 * id:    unique identification tag for the sph particle
 * fi:    flag if particle is freezing but has not yet freezed-out
 * fo:    flag if particle if freezed-out
 * ni:    reference density SPH weight
 * rho:   lab frame reference density 
 * rho_p: local rest frame reference density
 * S:     entropy SPH weight - particle integrated entropy
 * s:     lab frame entropy density
 * s_p:   local rest frame entropy density
 * sigma: entropy souce -> \partial_\mu(s u^\mu) = \sigma
 * q:     electric charge sign
 * Nb:    Barionic density SPH weight
 * Nc:    Electric charge density SPH density
 * nb:    Barionic density
 * nc:    Electric charge density
 * T:     termodynamic temperature
 * e_p:   proper energy density/thermodynamic energy density
 * p_p:   (proper) pressure/thermodynamic pressure
 * h_p:   thermodynamic enthalpy -> h_p = e_p + p_p
 * hsh_p: s \frac{\partial h_p}{\partial s} - h
 * x:     SPH particle spacetime position
 * u:     SPH particle 4-velocity
 * v:     SPH particle velocity -> v[l] = u[l]/u[0]
 * 
 * ---
 * 
 * Freeze-out members:
 * 
 * T_, rhop_,sp_, x_, u_ : use for 3 point reconstruction, 
 *                         n-point in the general case
 *                         they accout for the value of the
 *                         associated quantities on the previous time
 *                         steps. 
 */
typedef struct SPHeq_particle
{
	int D,id,*ind,fo; 
	double ni,rho,q,S,s,s_p,sigma,Nb,Nc;
	double rho_p,e_p,p_p,h_p,hsh_p,T,nb,nc;
	double *x,*u,*v;
    double Ta,rhopa,spa,*xa,*ua,*gradsa;
    double Tb,rhopb,spb,*xb,*ub,*gradsb;
    double Tc,rhopc,spc,*xc,*uc,*gradsc;
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
