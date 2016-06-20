/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#ifndef _AMI_POL_H_
#define _AMI_POL_H_

#include <cstdlib>
#include <math.h>

long double ami_horner(long double *pol,int degree,long double x,
                       long double *fp);
long double ami_root_bisection(long double *pol,int degree,long double a,
                               long double b,long double TOL);
int ami_polynomial_root(double *pol, int degree, double *root_r,
                        double *root_i);

#endif
