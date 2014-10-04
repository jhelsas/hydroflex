#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

using namespace std;
 
/*
 * type = 0 : D dimensional cubic domain
 * type = 1 : 2 dimensional triangular domain
 * type = 2 : 3 dimensional tetraedric domain
 */
 
class Domain{
  public:
    static int D;
    static vector<Domain> dms; 
    int Nv,type,good;
    vector <double> xv;
    double S;
  
    int init(int,int,int);
    Domain();
    Domain(int,int,int);
    Domain (int,int,int,double*);
    Domain (int,int,int,double*,double*);
    ~Domain();
    int domain_count();
    int domain_check();
    int el_print();
    int print_domain_n(int n);
    int print_domain();
    int init_cube(int dim,double xl[],double xu[]);
    int init_triangle(int dim,double Ntri,double *xv);
    int check_inside(double *x);
    int domain_split(double cutoff,gsl_monte_function F);
    int clean_domain();
    int print_sph(const char *filename);
    double aspect_ratio();
    int cubic_split(int n);
    int bc_simplex_split(int n);
    int bc_coord(double *r,double *lmb);
    int bc_test(double *r);
    int triangle_midpoint_split(int n);
};

    
typedef struct wparams{
  Domain *mdel;
  int id;
  void *p;
} wparams;


/*******************************************************
 * Domain Related Functions  
 ***************************************************************/

int rollin(int *ind,int *indx,int D,int n);

int rollout(int *ind,int *indx,int D,int n);

int bc_check(int D,double *lmb);
 
int unit_hexagon(int Ntri,int D,double *xv);

int unit2_hexagon(int Ntri,int D,double *xv);

int create_grid(int D,double **xpo,double *xl,double *xu,
                double *dx,int *Np);
                
int sph_read(char* filename,int *Dout,int *Nout,double **x,double **u,double **S);

double w_bspline(double r,double h);

int sph_dens(int D,int N,int Npoints,double *xp,double *x,
             double *S,double h,double xl[],double xu[],
             double (*tf)(double*,size_t,void*),char* filename,void *p);
