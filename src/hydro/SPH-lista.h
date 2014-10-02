/*
 * Funcoes de Manipulacao de lista e de preparacao importadas do 
 * programa SSPH. apenas setup_sph foi modificada.
 */

double maxD(double *a,int D);

double minD(double *a,int D);

int SPHeq_list_blen(SPHeq_list *lis);

int SPHeq_list_ilen(SPHeq_list *lis);

int SPHeq_append_bhead(SPHeq_list **lis,SPHeq_list *el);

int SPHeq_append_ihead(SPHeq_list **lis,SPHeq_list *el);

int SPHeq_clean_list(int N_sph,SPHeq_list *lsph);

int SPHneq_list_blen(SPHneq_list *lis);

int SPHneq_list_ilen(SPHneq_list *lis);

int SPHneq_list_len(SPHneq_list *lis);

int SPHneq_append_bhead(SPHneq_list **lis,SPHneq_list *el);

int SPHneq_append_ihead(SPHneq_list **lis,SPHneq_list *el);

int SPHneq_clean_list(int Nspecies,int *N,SPHneq_list **sph_neq);

int clean_box(Box *lbox);

int null_box(Box *lbox);

int linear_index(Box *lbox);

int set_box(int N_sph,SPHeq_list *sph_eq,int Nspecies,int *N,SPHneq_list **sph_neq,Box *lbox);

int fetch_inter(double *x, Box *lbox,double kh,SPHeq_list **inter_eq,SPHneq_list ***inter_neq);

int fetch_inter_eq(double *x, Box *lbox,double kh,SPHeq_list **inter_eq);

int fetch_inter_neq(double *x, Box *lbox,double kh,int k,SPHneq_list ***inter_neq);

int setup_sph(int D,double t,double h,double kh,int N_sph,SPHeq_list *sph_eq,
              int Nspecies, int *N, SPHneq_list **sph_neq,
              Box *lbox,double (*w)(double,double),int (*EoS)(SPHeq_particle*));

int ssph_2p1bi(int D,double t,double h,double kh,int N_sph,SPHeq_list *sph_eq,
               int Nspecies, int *N, SPHneq_list **sph_neq,
               Box *lbox,double (*w)(double,double),int (*EoS)(SPHeq_particle*));
