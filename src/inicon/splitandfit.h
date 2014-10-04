/*
 * type = 0 : D dimensional cubic domain
 * type = 1 : 2 dimensional triangular domain
 * type = 2 : 3 dimensional tetraedric domain
 */
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

using namespace std;
 
class domain{
  public:
    int D,Nv,type,good;
    vector <double> xv; 
    double S;
  
    int init(int,int,int);
    domain();
    domain(int,int,int);
};

typedef struct wparams{
  domain *mdel;
  void *p;
} wparams;


/*******************************************************
 * Domain Related Functions  
 *
 * Public Functions: to use as 
 *   interface of the library
 * 
 * Internal Functions: to use internally and not
 *   should be called outside the library, 
 *   except for debbuging purposes
 * 
 ********************************************************/

/*
 * Public Functions
 */
 
int domain_count(vector <domain> dom);

int domain_check(vector <domain> dom);

int el_print(domain el);

int print_domain_n(int D,vector <domain> dom,int n);

int print_domain(int D,vector <domain> dom);

int check_inside(double *x,size_t D,domain *mdel);

int init_cube(double xl[],double xu[],vector <domain> & dom,int D);

int init_triangle(int D,double Ntri,double *xv,vector <domain> & dom);

int domain_split(int D,double cutoff,vector <domain>& dom, gsl_monte_function F);

int clean_domain(vector <domain> &dom);

int print_sph(int D,const char *filename,vector <domain> &dom);

int print_moving_sph(int D,const char *filename,vector <domain> &dom,
                    int (*velocity)(double *,size_t,void *,double *),
                    void *par);

/*
 * Internal Functions
 */

int rollin(int *ind,int *indx,int D,int n);

int rollout(int *ind,int *indx,int D,int n);

double aspect_ratio(domain mdel);

int cubic_split(int D,vector <domain> & dom, int n);

int bc_simplex_split(int D,vector <domain> & dom,int n);

int bc_coord(int D,double *r,double *lmb,domain mdel);

int bc_check(int D,double *lmb);

int bc_test(int D,double *r,domain mdel);

int triangle_midpoint_split(int D,vector <domain> & dom,int n);

/**************************************************************
 *  Utility Functions  
 * 
 *  Usefull but not essential for the functioning of the
 *   library. Init and check/print functions
 * 
 **************************************************************/
 
int unit_hexagon(int Ntri,int D,double *xv);

int unit2_hexagon(int Ntri,int D,double *xv);

int create_grid(int D,double **xpo,double *xl,double *xu,
                double *dx,int *Np);
                
int sph_read(char* filename,int *Dout,int *Nout,double **x,double **u,double **S);

double w_bspline(double r,double h);

int sph_dens(int D,int N,int Npoints,double *xp,double *x,
             double *S,double h,double xl[],double xu[],
             double (*tf)(double*,size_t,void*),char* filename,void *p);
             
int sph_density_ploting(int D,double tau,int Npoints,double *xp,double xl[],double xu[],
                        int N,double *x,double *u,double *S,double h,
                        double (*ts)(double*,size_t,void*),
                        double (*ts_p)(double*,size_t,void*),
                        double (*te_p)(double*,size_t,void*),char* filename,void *p);
