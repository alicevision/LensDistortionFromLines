/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


/**
 * \file ami_pol.c
 * \brief library to estimate real polynomial roots.
 * \author Luis Alvarez
*/

#include "ami_pol.h"
#include "stdlib.h"
#include "stdio.h"


/**
 * \fn long double ami_horner(long double *pol,int degree,long double x,
                              long double *fp)
 * \brief function to evaluate a polynomial and its derivate using Horner
          Algorithm
 * \return polynomial evaluation
 * \author Luis Alvarez
*/
long double ami_horner(
    long double *pol /** coefficient polynomial vector
                          pol[0]+pol[1]*x+pol[2]*x^2+... */,
    int degree /** polynomial degree */,
    long double x /** point to evaluate the polynomial */,
    long double *fp /** output polynomial derivate evaluation */)
{
    long double PPX=pol[degree],PX=pol[degree];
    int k;
    for(k=degree-1; k>0; k--) {
        PX=PX*x+pol[k];
        PPX=PPX*x+PX;
    }
    PX=PX*x+pol[0];
    *fp=PPX;
    return(PX);
}

/**
 * \fn long double ami_root_bisection(long double *pol,int degree,long double a,
       long double b,long double TOL)
 * \brief function to compute the polynomial root in an interval. We use a
          combination of bisection and Newton-Raphson techniques
 * \return polynomial root in the interval
 * \author Luis Alvarez
*/
long double ami_root_bisection(
    long double *pol /** coefficient polynomial vector
                          pol[0]+pol[1]*x+pol[2]*x^2+... */ ,
    int degree  /** polynomial degree */,
    long double a  /** left interval extreme */,
    long double b /** right interval extreme */,
    long double TOL /** convergence accuracy*/)
{
    long double a2=a,b2=b,fa,fb,c,c0,fc,fp=0.;
    int iter=0,maxiter=100000;

    /* we evaluate the polynomial at interval extrema */
    fa=ami_horner(pol,degree,a,&fp);
    fb=ami_horner(pol,degree,b,&fp);
    /* we check if there is a root in the interval */
    if(fa*fb>0) return(0.5*(a+b));
    /* we evaluate the interval middle point */
    c=(a2+b2)*0.5;
    c0=c+1e30;
    /* we perform Newton-Raphson or Bisection iterations */
    while( (b2-a2)>(TOL*fabs(b2)+1e-10) &&  fabs(c-c0)>(TOL*fabs(c0)+1e-10) &&
          (iter++)<maxiter) {
        fc=ami_horner(pol,degree,c,&fp);
        /* we try Newton-Raphson */
        c0=c;
        if(fp!=0.) {
            c=c-fc/fp;
            /* printf("c=%f\n",(double) c); */
            if( c>=a2 && c<=b2 ) continue;
        }
        /* if Newton-Raphson fails, we try bisection */
        c=(a2+b2)*0.5;
        fc=ami_horner(pol,degree,c,&fp);
        c0=c+1e30;
        /* we check if we have arrive to the root */
        if(fabs(fc)<1e-20) break;
        /* we perform bisection iteration */
        if((fa*fc)<0) {
            b2=c;
            fb=fc;
        }
        else {
            a2=c;
            fa=fc;
        }
        c=(a2+b2)*0.5;
    }
    if(iter>=maxiter){return((a2+b2)*0.5);}
    return(c);
}

/**
 * \fn int ami_polynomial_root(double *pol, int degree, double *root_r,
                               double *root_i)
 * \brief function to estimate real polynomial roots. We use a recursive
          procedure computing the derivative polynomial roots.
          This function only works for relative low polynomial degree
 * \return The number of real polynomial root
 * \author Luis Alvarez
*/
int ami_polynomial_root(
    double *pol /** coefficient polynomial vector pol[degree]+pol[degree-1]*
                x+pol[degree-2]*x^2+...
                FOR HISTORICAL REASONS THE ORDER POLYNOMIAL COEFICIENT IS GIVEN
                IS DIFFERENT THAT IN THE OTHERS POLYNOMIAL FUNCTIONS */,
    int degree /** polynomial degree */,
    double *root_r /** output real component of polynomial roots */,
    double *root_i /** output complex component of polynomial roots. Since this
                 function only real roots this complex part is fitted to 0.
                 FOR HISTORICAL REASONS WE KEEP THIS COMPLEX PART IN THE
                 FUNCTION*/)
{

    long double a,b,*p,*p_aux,*ap,*f,fp,*pol2,fa,fb,TOL=1e-14;
    long double max,n_factor;
    int m,k,l,Nr,N=degree;

    p=p_aux=ap=f=pol2=NULL;
    fp=0.0;
    a=b=fa=fb=max=n_factor=0.0; /*Added by Luis Gomez*/
    Nr=0;/*Added by Luis Gomez*/

    /* WE ALLOCATE MEMORY */
    p=(long double*)malloc(sizeof(long double)*(N+2));
    p_aux=(long double*)malloc(sizeof(long double)*(N+2));
    ap=(long double*)malloc(sizeof(long double)*(N+1));
    f=(long double*)malloc(sizeof(long double)*(N+1));
    pol2=(long double*)malloc(sizeof(long double)*(N+2));

    /* POLYNOMIAL NORMALIZATION */
    if(pol[0]!=1.) {
        for(k=0; k<=N; k++) pol[k]/=pol[0];
    }

    /* WE REORDER THE POLYNOM */
    for(k=0; k<=N; k++) pol2[k]=pol[N-k];
    pol2[N+1]=0.0; /*Added by Luis Gomez*/

    /* WE NORMALIZE POLINOMIAL COEFFICIENTS TO AVOID POLYNOMIAL EVALUATION IN
       LARGE NUMBERS*/
    n_factor=0;
    for(k=0; k<N; k++) {
        if(fabs(pol2[k])>10.) {
            max=log(fabs((double) pol2[k])/10.)/(N-k);
            if(max>n_factor) n_factor=max;
        }
    }

    n_factor=exp((double) n_factor);
    max=n_factor;
    for(k=N-1; k>=0; k--) {
        pol2[k]/=max;
        max*=n_factor;
    }

    /* WE COMPUTE FACTORIAL NUMBERS */
    f[0]=1.;
    for(k=1; k<=N; k++) f[k]=f[k-1]*k;

    /* WE COMPUTE THE INITIAL INTERVAL WHERE ALL ROOTS ARE INCLUDED */
    max=fabs(pol2[0]);
    for(k=1; k<N; k++) {
        if(fabs(pol2[k])>max)
            max=fabs(pol2[k]);
    }
    /* WE INITIALIZE THE INTERVALS */
    p[0]=-1.-max/fabs(pol2[N]);
    /* WE COMPUTE THE POLYNOMIAL DEGREE-1 DERIVATE ROOT */
    p[1]=-(pol2[N-1]*f[N-1])/(pol2[N]*f[N]);
    for(k=2; k<=(N+1); k++) p[k]=-p[0];
    for(k=0; k<=N; k++)p_aux[k]=p[k];
    p_aux[N+1]=0.0; /*Added by Luis Gomez*/
    /* WE COMPUTE THE DERIVATIVE POLINOMIAL ROOTS IN A RECURSIVE WAY */
    for(k=N-2; k>=0; k--) {
        /* we compute polynomial derivative coefficients */
        for(l=0; l<=(N-k); l++) {
            ap[l]=pol2[l+k]*(f[k+l]/f[l]);
        }
        /*we check the valid intervals to compute derivative polynomial roots*/
        m=1;
        fa=ami_horner(ap,N-k,p_aux[0],&fp);
        for(l=1; l<=(N-k); l++) {
            fb=ami_horner(ap,N-k,p_aux[l],&fp);
            if((fa*fb)<=0 || fabs(fb)<=TOL) {
                p[m++]=p_aux[l];
                fa=fb;
            }
        }
        for(l=m; l<N; l++) p[l]=p[N];
        /* we compute the derivative polynomial roots in each interval */
        for(l=1; l<=(N-k); l++) {
            p_aux[l]=ami_root_bisection(ap,N-k,p[l-1],p[l],TOL);
        }
    }

    /* WE STORE THE ROOTS */
    root_i[0]=0.; /* we fit the complex component of the root to 0 */
    Nr=0;
    for(k=1; k<=N; k++) {
        /* we check if the polynomial has a root in the point */
        a=(p_aux[k]+p_aux[k-1])*0.5;
        b=(p_aux[k]+p_aux[k+1])*0.5;
        fa=ami_horner(pol2,degree,a,&fp);
        fb=ami_horner(pol2,degree,b,&fp);
        if(fa*fb<0 || fabs(ami_horner(pol2,degree,p_aux[k],&fp))<TOL) {
            root_i[Nr]=0; /* we fit the complex component of the root to 0 */
            root_r[Nr++]=p_aux[k]*n_factor; /* we denormalize the root */
        }
    }
    /* we free the memory */
    free(p);
    free(ap);
    free(f);
    free(pol2);
    free(p_aux);
    /* we return the number of real roots estimated */
    return(Nr);

}
