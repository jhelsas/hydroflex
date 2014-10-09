int HE2(int D,double t,double dt,double h,double kh,
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
              SPHeq_list *,SPHneq_list **),
        double *Pt
       );
       
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
       );
       
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
       );
    
int Deriv_MCV(double t,int D,int N_sph,int Nspecies,int *N,double h,double kh,
              SPHeq_list *sph_eq,SPHneq_list **sph_neq,Box *lbox,
              double (*w)(double,double),double (*Dw)(double,double),
              SPHeq_list *f_eq,SPHneq_list **f_neq
             );
    
int Drv_2p1bi(double t,int D,int N_sph,int Nspecies,int *N,double h,double kh,
              SPHeq_list *sph_eq,SPHneq_list **sph_neq,Box *lbox,
              double (*w)(double,double),double (*Dw)(double,double),
              SPHeq_list *f_eq,SPHneq_list **f_neq
             );

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
             );
