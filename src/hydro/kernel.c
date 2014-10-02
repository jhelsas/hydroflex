/*
 *      kernel.c
 *      
 *      Copyright 2009 José Hugo Elsas <jhelsas@jhelsas-desktop>
 *      

 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MYPI 3.1415926
/*#define A_d 1.*/  /*dim = 1*/
#define A_d (15.)/(7.*MYPI) /*dim=2*/
/*#define A_d (3.0)/(2.0*MYPI) */ /*dim=3*/

#define DIM 2

double signum(double z)
{
	if(z>=0) return 1.0;
	else return (-1.0);
}

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

double Dw_bspline(double r,double h)
{
	double R;
	if(h<=0.) exit(10);
	R=fabs(r)/h;
	if(R>=2.)
		return 0;
	else if((1.<=R)&&(R<2.))
		return (signum(r)*(-A_d*(2.0-R)*(2.0-R)))/(2.0*(h*pow(h,DIM)));
	else
		return (signum(r)*(A_d*(-2.0*R+(3./2.)*R*R)))/(h*pow(h,DIM));
	
}

double DDw_bspline(double r,double h)
{
	double R;
	if(h<=0.) exit(10);
	R=fabs(r)/h;
	if(R>=2.)
		return 0;
	else if((1.<=R)&&(R<2.))
		return ((A_d)*(2.0-R))/((h*h*pow(h,DIM)));
	else
		return (A_d*(-2.0+3.0*R))/(h*h*pow(h,DIM));
	
}


#define e_bspline 10

double w_bspline_calcado(double r,double h)
{
	return (w_bspline(r,h)+e_bspline)/(1+4.0*h*e_bspline);
}

double Dw_bspline_calcado(double r,double h)
{
	return Dw_bspline(r,h)/(1+4.0*h*e_bspline);
}

double DDw_bspline_calcado(double r,double h)
{
	return DDw_bspline(r,h)/(1+4.0*h*e_bspline);
}

/*********************/

double w_gauss(double r,double h)
{
	double R;
	if(h<=0.0) exit(10);
	R=fabs(r)/h;
	if(R>=2.0)
		return 0.0;
	else
		return (3.0*exp(-(4.5*R*R))/(sqrt(2*MYPI)*h));
}

double Dw_gauss(double r,double h)
{
	double R;
	if(h<=0.) exit(10);
	R=fabs(r)/h;
	if(R>=2.)
		return 0;
	else
		return -signum(r)*((27.0*R*exp(-(4.5*R*R)))/(sqrt(2*MYPI)*h*h));
}

double DDw_gauss(double r,double h)
{
	double R;
	if(h<=0.) exit(10);
	R=fabs(r)/h;
	if(R>=2.)
		return 0;
	else
		return ((27.0*(9.0*R*R-1.0)*exp(-(4.5*R*R)))/(sqrt(2*MYPI)*h*h*h));
}

double w_gauss_calcado(double r,double h)
{
	return (w_gauss(r,h)+e_bspline*w_gauss(0,h))/(1+4.0*h*e_bspline*w_gauss(0,h));
}

double Dw_gauss_calcado(double r,double h)
{
	return Dw_gauss(r,h)/(1+4.0*h*e_bspline*w_gauss(0,h));
}

double DDw_gauss_calcado(double r,double h)
{
	return DDw_gauss(r,h)/(1+4.0*h*e_bspline*w_gauss(0,h));
}

/************************************************/

double w_sgauss(double r,double h)
{
	double R;
	if(h<=0.0) exit(10);
	R=fabs(r)/h;
	if(R>=2.0)
		return 0.0;
	else
		return ((1.5-R*R)*exp(-R*R))/(sqrt(MYPI)*h);
}

double Dw_sgauss(double r,double h)
{
	double R;
	if(h<=0.) exit(10);
	R=fabs(r)/h;
	if(R>=2.)
		return 0;
	else
		return -signum(r)*(((2.0*(R*R*(R-1)))*exp(-R*R))/(sqrt(MYPI)*h*h));
}

double DDw_sgauss(double r,double h)
{
	double R;
	if(h<=0.) exit(10);
	R=fabs(r)/h;
	if(R>=2.)
		return 0;
	else
		return ((2.0*(3.0*R*R-2.0*R+2.0*R*R*R*R-R*R*R)*exp(-R*R))/(sqrt(MYPI)*h*h*h));
}
