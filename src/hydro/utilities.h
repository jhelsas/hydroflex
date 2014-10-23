#define PI (3.1415926)

typedef struct funDDtD
{
	double (*f)(int,double,double);
} funDDtD;

int sph_eqMalloc(int N_sph,int D,SPHeq_list **sph_eq);

int sph_neqMalloc(int Nspecies,int *N,int D,SPHneq_list ***sph_neq);

int sph_eqFree(int N_sph,SPHeq_list **sph_eq);

int sph_neqFree(int Nspecies,int *N,SPHneq_list ***sph_neq);

int init_simulation(char *config_file,int *D_out,
                    int *N_sph_out,int *Nspecies_out,int **N_out,
                    double *h_out,double *kh_out,
                    double *t0_out,double *tf_out,double *dt_out,
                    SPHeq_list **sph_eq_out,SPHeq_list **sph_eqTemp_out,
                    SPHeq_list **f0_eq_out,SPHeq_list **f1_eq_out,
                    SPHeq_list **f2_eq_out,SPHeq_list **f3_eq_out,
                    SPHneq_list ***sph_neq_out,SPHneq_list ***sph_neqTemp_out,
                    SPHneq_list ***f0_neq_out,SPHneq_list ***f1_neq_out,
                    SPHneq_list ***f2_neq_out,SPHneq_list ***f3_neq_out,
                    Box **lbox_out);

int init_printing_grid(int Npoints,int D,char *filename,double ***xp_out);

int create_grid(int D,double **xpo,double *xl,double *xu,
                double *dx,int *Np);

int create_X_line(int D,double **xpo,double *xl,double *xu,
                  double *dx,int *Np);

int create_Y_line(int D,double **xpo,double *xl,double *xu,
                  double *dx,int *Np);
                  
int print_proper_entropy(int D,double h,double kh,Box *lbox,double (*w)(int,double,double),
                         char *filename,int Npoints,double **xp);

int print_proper_energy(int D,int Nspecies,double h,double kh,Box *lbox,
                        double (*w)(int,double,double),char *filename,int Npoints,double **xp,double *displ);
                        
int print_new(int D,int Nspecies,double t,double h,double kh,Box *lbox,
              double (*w)(int,double,double),char *filename,
              int Npoints,double *xp,double *displ);
              
int print_4vel_profile(int D,int Nspecies,double h,double kh,Box *lbox,
                       double (*w)(int,double,double),char *filename,int Npoints,double **xp);
                       
int scripting(int argc,char **argv,
              double t0,double tf,double dt);
