double fd(double *x,size_t n,void *params);

double winicon(double x[],size_t dim,void *par);

double gauss_whs(double x[],size_t dim,void *par);

double woodsaxon(double x[],size_t dim,void *par);

double woodsaxon_whotspot(double x[],size_t dim,void *par);

double gubser_entropy(double *x,size_t dim,void *par);

double gubser_proper_entropy(double *x,size_t dim, void *par);

double gubser_proper_energy(double *x,size_t dim, void *par);

int null_velocity(double *x,size_t dim,void *par,double *u);

int gubser_velocity(double *x,size_t dim,void *par,double *u);
