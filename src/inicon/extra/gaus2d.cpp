double inicon(double x[],const int dim,void * par) {

  if(dim != 2) {
    fprintf(stderr, "error: dim != 2");
    abort();
  }  

  double xx = x[0];
  double yy = x[1];

  double *p = (double*)par;
  double    A = p[0];
  double   x0 = p[1];
  double   y0 = p[2];
  double sigx = p[3];
  double sigy = p[4];

  return A*exp(-0.5*(pow((xx-x0)/sigx,2.)+pow((yy-y0)/sigy,2.)));

}

