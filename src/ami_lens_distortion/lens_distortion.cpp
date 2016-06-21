/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */

 #ifndef AMI_DLL_CPP
  #define AMI_DLL_CPP
#endif

/**
 * \file lens_distortion.cpp
 * \brief Functions for lens distortion model basic operations
 * \author Luis Alvarez \n \n
 */


#include <stdio.h>
#include <math.h>
#include "lens_distortion.h"
#include "../ami_utilities/utilities.h"
#include "../ami_pol/ami_pol.h"
/**
 * \fn ami_lens_distortion_polynomial_update_distance_2v(double *x, double *y,int Np,double *a,int Na,double x0,double y0,int k1,int k2,double **pol,double alfa)
 *  \brief Function to add the information of a line point sequence to the 4 degree polynomial to compute the lens distortion model
 *  \pre Any parameter can be null.
 *  \param[in]  x,y Distorted line coordinates
 *  \param[in] Np Number of points
 *  \param[in] a Polynomial defining the lens distortion model
 *  \param[in] Na Degree of polynomial model for lens distortion
 *  \param[in] x0,y0 Coordinates of the image center
 *  \param[in] k1 Coeficient 1 of the polynomial to be updated
 *  \param[in] k2  Coeficient 2 of the polynomial to be updated
 *  \param pol 4 degree 2 variable polynom to minimize
 *  \param[in] alfa  Weight of the distance in the polynom energy
 *  \return 0

 */
AMI_DLL_CPP int ami_lens_distortion_polynomial_update_distance_2v(
   double *x, double *y,  /* DISTORTED LINE COORDINATES (INPUT)*/
   int Np, /* NUMBER OF POINTS (INPUT)*/
   double *a, /* POLYNOMIAL DEFINING THE LENS DISTORTION MODEL (INPUT) */
   int Na,   /* DEGREE OF POLYNOMIAL MODEL FOR LENS DISTORTION */
   double x0,double y0, /* COORDINATES OF THE IMAGE CENTER */
   int k1,  /* COEFICIENT 1 OF THE POLYNOMIAL TO BE UPDATED*/
   int k2,  /* COEFICIENT 2 OF THE POLYNOMIAL TO BE UPDATED*/
   double **pol, /* 4 DEGREE 2 VARIABLE POLYNOM TO MINIMIZE (INPUT-OUTPUT) */
   double alfa) /* WEIGHT OF THE DISTANCE IN THE POLYNOM ENERGY */
{
  int i, j, k/*,l*/;
  double *A, *x2, *y2, *d1, *d2/*,s_yy,y_m,x2_m,s_max,x_m,s_xx,y2_m ,xA_m*/;
  //double /*xd1_m,yA_m,yd2_m,xd2_m,yd1_m*/;
  double paso, **pol1, **pol2, **p_xx, **p_xy, **p_yy;
  /* WE CHECK alfa VALUE */
  if(Np < 3)
    return(-1);
  if(alfa == 0)
    return(-2);

  /* WE ALLOCATE MEMORY */
  A=(double*)malloc( sizeof(double)*Np );
  x2=(double*)malloc( sizeof(double)*Np );
  y2=(double*)malloc( sizeof(double)*Np );
  d1=(double*)malloc( sizeof(double)*Np );
  d2=(double*)malloc( sizeof(double)*Np );
  ami_calloc2d(pol1,double,2,2);
  ami_calloc2d(pol2,double,2,2);
  ami_calloc2d(p_xx,double,3,3);
  ami_calloc2d(p_xy,double,3,3);
  ami_calloc2d(p_yy,double,3,3);

  /* WE COMPUTE THE DISTANCE TO THE IMAGE CENTER */
  for(i=0; i<Np; i++)
    d1[i]=sqrt( (x[i]-x0)*(x[i]-x0)+(y[i]-y0)*(y[i]-y0) );

  /* WE COMPUTE THE POINT TRANSFORMATION WITH THE CURRENT LENS DISTORTION MODEL */
  for(i=0; i<Np; i++)
  {
    A[i]=ami_polynomial_evaluation(a,Na,d1[i]);
    x2[i]=x0+(x[i]-x0)*A[i];
    y2[i]=y0+(y[i]-y0)*A[i];
  }


 /* WE COMPUTE THE POLYNOMS CORRESPONDING TO THE DISTANCE ERROR */

 for(i=0; i<=2; i++)
  for(j=0; j<=2; j++) 
    p_xx[i][j]=0.;
  for(i=0; i<Np; i++)
  {
    paso=0;
    for(k=1; k<=Na; k++) 
      if(k!=k1 && k!=k2)
        paso+=a[k]*pow(d1[i],(double) k);
    pol1[0][0]=paso*d1[i];
    pol1[1][0]=pow(d1[i],(double) k1+1);
    pol1[0][1]=pow(d1[i],(double) k2+1);
    ami_2v_polynom_multiplication(pol1,1,pol1,1,p_xx);
  }
  for(i=0; i<=2; i++)
    for(j=0; j<=2; j++)
      p_xx[i][j]=alfa*p_xx[i][j]/Np;

  /* WE UPDATE THE ERROR POLYNOM */
  ami_2v_polynom_multiplication(p_xx,2,p_xx,2,pol);

  /* WE FREE THE MEMORY */
  free(A); free(x2); free(y2); free(d1); free(d2);
  ami_free2d(p_xx); ami_free2d(p_xy);  ami_free2d(p_yy);
  ami_free2d(pol1); ami_free2d(pol2);

  return(0);
}



/**
 * \fn double ami_lens_distortion_estimation_2v(double **x,double **y,int Nl,int *Np,double x0,double y0,double *a,int Na,int k1,int k2,double alfa)
 * \brief Update of the lens distortion polynomial model for 2 variables. If alfa>0 we adapt a[0] to minimize the square distance beewten distorted and undistorted points and we add a term to the polynomial also minimizing such distance with weight alfa
 * \pre Any parameter can be null.
 * \param [out]  x,y Distorted line coordinates
 * \param [out] Np Number of points
 * \param [in] x0,y0 Coordinates of the image center
 * \param [in] a Polynomial defining the lens distortion model
 * \param [in] Na Degree of polynomial model for lens distortion
 * \param [in] k1 Coeficient 1 of the polynomial to be updated
 * \param [in] k2  Coeficient 2 of the polynomial to be updated
 * \param [in] alfa  Weight of the distance in the polynom energy
 * \param Nl  Not described.
 * \return Return the distance
 */
AMI_DLL_CPP double ami_lens_distortion_estimation_2v(double **x,double **y,int Nl,int *Np,
   double x0,double y0,double *a,int Na,int k1,int k2,double alfa)
{

  int i,/*j,*/k, m;
  double **pol_v2, d, suma_Ad, suma_dd,/*paso,*/A, Error=0;
  /* WE ALLOCATE MEMORY */
  ami_calloc2d(pol_v2,double,5,5);

  /* WE UPDATE a[0] BY MINIMIZING THE DISTANCE OF THE DISTORTED POINTS TO
  THE UNDISTORTED POINTS */
  if(alfa > 0)
  {
    suma_dd=suma_Ad=0;
    for(m=0; m<Nl; m++)
    {
      for(i=0; i<Np[m]; i++)
      {
        d=sqrt( (x[m][i]-x0)*(x[m][i]-x0)+(y[m][i]-y0)*(y[m][i]-y0) );
        A=0;
        for(k=1; k<=Na; k++)
          A+=a[k]*pow(d,(double) k+1);
        suma_dd+=d*d;
        suma_Ad+=A*d;
      }
    }
    a[0]=1-suma_Ad/suma_dd;
    //printf("\n a[0] ACTUALIZADO = %e\n",a[0]);
  }


  for(m=0; m<Nl; m++)
  {
    /* WE UPDATE DE POLYNOM TO MINIMIZE */
    ami_lens_distortion_polynomial_update_2v(x[m],y[m],Np[m],a,Na,x0,y0,k1,k2,pol_v2);
    ami_lens_distortion_polynomial_update_distance_2v(x[m],y[m],Np[m],a,Na,x0,y0,k1,k2,pol_v2,alfa);

  }

  /* WE UPDATE THE POLYNOMIAL LENS DISTORTION MODEL */
  ami_lens_distortion_model_update_2v(a,k1,k2,pol_v2);

  ami_free2d(pol_v2);


  for(m=0; m<Nl; m++)
    Error+=ami_LensDistortionEnergyError(x[m],y[m],Np[m],x0,y0,a,Na);
  return(Error/Nl);
}

/**
 * \fn int ami_lens_distortion_model_update_2v(double *a,int Na,int k1,int k2,double **pol)
 *  \brief Function to update the lens distortion model by minimizing a 4 degree  2 variable polynom
 *  \pre Any parameter can be null.
 *  \param a Polynomial defining the lens distortion model
 *  \param [in] Na Degree of polynomial model for lens distortion
 *  \param [in] k1 Coeficient 1 of the polynomial to be updated
 *  \param [in] k2  Coeficient 2 of the polynomial to be updated
 *  \param [in] pol  4 degree polynom to minimize
 *  \return 0
 */
AMI_DLL_CPP int ami_lens_distortion_model_update_2v(
   double *a, /* POLYNOMIAL DEFINING THE LENS DISTORTION MODEL (INPUT-OUTPUT) */
   int k1,  /* COEFICIENT 1 OF THE POLYNOMIAL TO BE UPDATED (INPUT)*/
   int k2,  /* COEFICIENT 2 OF THE POLYNOMIAL TO BE UPDATED (INPUT)*/
   double **pol) /* 4 DEGREE POLYNOM TO MINIMIZE (INPUT) */
{
  int j, i, M, Nr=0, m;
  double *x, **pol_x, **pol_y, *pol_r, xr, yr, Emin, *rx, *ry, *b2, *p3;
  double sx, sy,/*paso,*/ Energy; /* NORMALIZATION FACTORS */
  double p_r3[6][6][19];
  for(i=0; i<6; i++)
    for(j=0; j<6; j++)
      for(m=0; m<19; m++)
        p_r3[i][j][m]=0.;
  /* WE ALLOCATE MEMORY */
  x=(double*)malloc( sizeof(double)*3 );
  ami_calloc2d(pol_x,double,5,5);
  ami_calloc2d(pol_y,double,5,5);
  p3=(double*) malloc(sizeof(double)*5);

  //for(i=0;i<=4;i++) for(j=0;j<=4;j++) printf("pol_o[%d][%d]=%e\n",i,j,pol[i][j]);
  /* WE NORMALIZE POLYNOM COEFICIENT */
  sx=pow(pol[4][0],(double) 0.25);
  sy=pow(pol[0][4],(double) 0.25);
  for(i=0;i<=4;i++)
  {
    for(j=0; j<=4; j++)
    {
      if(i > 0)
        pol[i][j]/=pow(sx,(double) i);
      if(j > 0)
        pol[i][j]/=pow(sy,(double) j);
    }
  }
  //for(i=0;i<=4;i++) for(j=0;j<=4;j++) printf("pol[%d][%d]=%e\n",i,j,pol[i][j]);

  /* WE COMPUTE THE DERIVATIVES OF THE POLYNOM */
  ami_2v_polynom_derivatives(pol,4,pol_x, pol_y);
  /* WE FILL THE MATRIX TO COMPUTE THE DETERMINANT */
  for(i=0; i<=3; i++)
  {
    for(m=0; m<=4; m++)
    {
      p_r3[2][i+2][m]=p_r3[1][i+1][m]=p_r3[0][i][m]=pol_x[3-i][m];
      p_r3[5][i+2][m]=p_r3[4][i+1][m]=p_r3[3][i][m]=pol_y[3-i][m];
    }
  }
  /* WE COMPUTE THE RESOLVENT POLYNOM */
  pol_r=(double*) malloc(sizeof(double)*19);
  ami_polynom_determinant(p_r3,18,6,pol_r);
  /* WE COMPUTE THE RESOLVENT POLYNOM DEGREE */
  for(i=0; i<=18; i++)
  {
    if(pol_r[i] != 0)
      Nr=i;
  }
  /* WE COMPUTE THE ROOT OF THE RESOLVENT POLYNOM */
  rx=(double*) malloc(sizeof(double)*Nr);
  ry=(double*) malloc(sizeof(double)*Nr);
  b2=(double*) malloc(sizeof(double)*(Nr+1));
  
  for(i=0; i<=Nr; i++)
    b2[i] = pol_r[Nr-i];
  Nr=ami_polynomial_root(b2,Nr,rx,ry);
  /* WE COMPUTE THE X COMPONENT BY REPLACING THE ROOTS IN THE DERIVATIVES
  OF THE POLYNOM */
  xr=0; 
  yr=0; 
  Emin=10e90;
  
  for(i=0; i<Nr; i++)
  {
    if(fabs(ry[i])> 0.000000000000001)
      continue;
    ami_2v_polynom_to_1v_polynom(pol_x,4,p3,rx[i],1);
    M=ami_RootCubicPolynomial(p3,3,x);
    
    for(m=0; m<M; m++)
    {
     Energy=ami_2v_polynom_evaluation(pol,4,x[m],rx[i]);
     if(Energy < Emin)
     { 
       Emin=Energy; 
       xr=rx[i]; 
       yr=x[m];
     }
    }
    ami_2v_polynom_to_1v_polynom(pol_y,4,p3,rx[i],1);
    M=ami_RootCubicPolynomial(p3,3,x);
    
    for(m=0; m<M; m++)
    {
      Energy=ami_2v_polynom_evaluation(pol,4,x[m],rx[i]);
      if(Energy < Emin)
      { 
        Emin=Energy;
        xr=rx[i];
        yr=x[m];
      }
    }
  }
  /* WE UPDATE THE DISTORSION POLYNOMIAL MODEL */
  a[k1]+=(yr/sx);
  a[k2]+=(xr/sy);

  /* WE FREE THE MEMORY */
  free(x); ami_free2d(pol_x); ami_free2d(pol_y); free(pol_r);
  free(p3); free(rx); free(ry); free(b2);

  return(0);
}

/**
 * \fn int ami_lens_distortion_polynomial_update_2v(double *x, double *y,int Np,double *a,int Na,double x0,double y0,int k1,int k2,double **pol)
 *  \brief Function to add the information of a line point sequence to the 4 degree polynomial to compute the lens distortion model
 *  \pre Any parameter can be null.
 *  \param [in] x,y Distorted line coordinates
 *  \param [in] Np Number of points
 *  \param [in] a Polynomial defining the lens distortion model
 *  \param [in] Na Degree of polynomial model for lens distortion
 *  \param [in] x0,y0 Coordinates of the image center
 *  \param [in] k1 Coeficient 1 of the polynomial to be updated
 *  \param [in] k2  Coeficient 2 of the polynomial to be updated
 *  \param pol  4 degree 2 variable polynom to minimize
 *  \return 0
 */
AMI_DLL_CPP int ami_lens_distortion_polynomial_update_2v(
   double *x, double *y,  /* DISTORTED LINE COORDINATES (INPUT)*/
   int Np, /* NUMBER OF POINTS (INPUT)*/
   double *a, /* POLYNOMIAL DEFINING THE LENS DISTORTION MODEL (INPUT) */
   int Na,   /* DEGREE OF POLYNOMIAL MODEL FOR LENS DISTORTION */
   double x0,double y0, /* COORDINATES OF THE IMAGE CENTER */
   int k1,  /* COEFICIENT 1 OF THE POLYNOMIAL TO BE UPDATED*/
   int k2,  /* COEFICIENT 2 OF THE POLYNOMIAL TO BE UPDATED*/
   double **pol) /* 4 DEGREE 2 VARIABLE POLYNOM TO MINIMIZE (INPUT-OUTPUT) */
{
  int i, j/*,k,l*/;
  double *A, *x2, *y2, *d1, *d2,/*x_m,y_m,*/x2_m, y2_m, s_xx, s_yy/*,s_max*/, xA_m;
  double xd1_m, yA_m, yd1_m, xd2_m, yd2_m;
  double paso, **pol1, **pol2, **p_xx, **p_xy, **p_yy;
  if(Np < 3)
    return(-1);
  /* WE ALLOCATE MEMORY */
  A=(double*)malloc( sizeof(double)*Np );
  x2=(double*)malloc( sizeof(double)*Np );
  y2=(double*)malloc( sizeof(double)*Np );
  d1=(double*)malloc( sizeof(double)*Np );
  d2=(double*)malloc( sizeof(double)*Np );
  ami_calloc2d(pol1,double,2,2);
  ami_calloc2d(pol2,double,2,2);
  ami_calloc2d(p_xx,double,3,3);
  ami_calloc2d(p_xy,double,3,3);
  ami_calloc2d(p_yy,double,3,3);

  /* WE COMPUTE THE DISTANCE TO THE IMAGE CENTER */
  for(i=0; i<Np; i++)
    d1[i]=sqrt( (x[i]-x0)*(x[i]-x0)+(y[i]-y0)*(y[i]-y0) );

  /* WE COMPUTE THE POINT TRANSFORMATION WITH THE CURRENT LENS DISTORTION MODEL */
  for(i=0; i<Np; i++)
  {
    A[i]=ami_polynomial_evaluation(a,Na,d1[i]);
    x2[i]=x0+(x[i]-x0)*A[i];
    y2[i]=y0+(y[i]-y0)*A[i];
  }

  /* WE COMPUTE THE DISTANCE POWER k1 AND k2 (THE COEFICIENT OF THE LENS
  DISTORTION MODEL TO BE UPDATED */
  for(i=0; i<Np; i++)
  {
    paso=d1[i];
    d1[i]=pow(paso,(double) k1);
    d2[i]=pow(paso,(double) k2);
 }

  /* WE COMPUTE THE VARIANCE OF THE TRANSFORMED POINTS */
  x2_m=0;
  for(i=0; i<Np; i++)
    x2_m+=x2[i]; 
  x2_m/=Np;
  
  s_xx=0;
  for(i=0; i<Np; i++)
    s_xx+=(x2[i]-x2_m)*(x2[i]-x2_m);
  s_xx/=Np;
  
  y2_m=0; for(i=0;i<Np;i++) y2_m+=y2[i];  y2_m/=Np;
  s_yy=0; for(i=0;i<Np;i++) s_yy+=(y2[i]-y2_m)*(y2[i]-y2_m); s_yy/=Np;
  //s_max=s_xx>s_yy?s_xx:s_yy;

  /* WE COMPUTE SOME AVERAGES WE NEED */
  xA_m=0; for(i=0;i<Np;i++) xA_m+=(x[i]-x0)*A[i]; xA_m/=Np;
  xd1_m=0; for(i=0;i<Np;i++) xd1_m+=(x[i]-x0)*d1[i]; xd1_m/=Np;
  xd2_m=0; for(i=0;i<Np;i++) xd2_m+=(x[i]-x0)*d2[i]; xd2_m/=Np;
  yA_m=0; for(i=0;i<Np;i++) yA_m+=(y[i]-y0)*A[i]; yA_m/=Np;
  yd1_m=0; for(i=0;i<Np;i++) yd1_m+=(y[i]-y0)*d1[i]; yd1_m/=Np;
  yd2_m=0; for(i=0;i<Np;i++) yd2_m+=(y[i]-y0)*d2[i]; yd2_m/=Np;

 /* WE COMPUTE THE POLYNOMS OF THE SECOND ORDER MOMENT OF THE POINT
     p_xx p_xy AND p_yy DISTRIBUTION */
  for(i=0; i<Np; i++)
  {
    pol1[0][0]=(x[i]-x0)*A[i]-xA_m;
    pol1[1][0]=(x[i]-x0)*d1[i]-xd1_m;
    pol1[0][1]=(x[i]-x0)*d2[i]-xd2_m;
    pol2[0][0]=(y[i]-y0)*A[i]-yA_m;
    pol2[1][0]=(y[i]-y0)*d1[i]-yd1_m;
    pol2[0][1]=(y[i]-y0)*d2[i]-yd2_m;
    ami_2v_polynom_multiplication(pol1,1,pol1,1,p_xx);
    ami_2v_polynom_multiplication(pol1,1,pol2,1,p_xy);
    ami_2v_polynom_multiplication(pol2,1,pol2,1,p_yy);
  }
  /* WE COMPUTE p_xx * p_yy - p_xy * p_xy  AND WE ADD THE VALUE TO THE
  GIVEN ENERGY ERROR POLYNOM */

  for(i=0; i<=2; i++)
    for(j=0; j<=2; j++)
      p_xx[i][j]/=1. /*s_max*/ ;
  ami_2v_polynom_multiplication(p_xx,2,p_yy,2,pol);
  for(i=0; i<=2; i++)
    for(j=0; j<=2; j++)
      p_xx[i][j]=-p_xy[i][j]/1. /*s_max*/;
  ami_2v_polynom_multiplication(p_xy,2,p_xx,2,pol);

  /* WE FREE THE MEMORY */
  free(A); free(x2); free(y2); free(d1); free(d2);
  ami_free2d(p_xx); ami_free2d(p_xy);  ami_free2d(p_yy);
  ami_free2d(pol1); ami_free2d(pol2);

  return(0);
}




/**
 * \fn void ami_2v_polynom_derivatives(double **p,int N,double **p_x,double **p_y)
 *  \brief Function to compute the partial derivatives of a 2 variable polynom the degree of the derivative polynoms is assumed to be the same that the original one .
 *  \pre Any parameter can be null.
 *  \param[in] p Original polynom
 *  \param[in] N Degree of the original polybom
 *  \param[out] p_x Derivative of the polynom with respect to the first variable
 *  \param[out] p_y Derivative of the polynom with respect to the second variable
 */
AMI_DLL_CPP void ami_2v_polynom_derivatives(
  double **p,/* ORIGINAL POLYNOM (INPUT)*/
  int N, /* DEGREE OF THE ORIGINAL POLYBOM (INPUT) */
  double **p_x, /* DERIVATIVE OF THE POLYNOM WITH RESPECT TO THE FIRST VARIABLE (OUTPUT) */
  double **p_y) /* DERIVATIVE OF THE POLYNOM WITH RESPECT TO THE SECOND VARIABLE(OUTPUT) */
{
  int i, j;
  for(i=0; i<=N; i++)
    for(j=0; j<=N; j++)
      p_x[i][j]=p_y[i][j]=0;

  for(i=1; i<=N; i++)
    for(j=0; j<=N; j++)
      p_x[i-1][j]=i*p[i][j];

  for(i=0; i<=N; i++)
    for(j=1; j<=N; j++)
      p_y[i][j-1]=j*p[i][j];

}

/**
 * \fn double ami_determinante(double **A,int N)
 *  \brief Function to evaluate the determinant of a matrix
 *  \pre Any parameter can be null.
 *  \param[in] A Matrix
 *  \param[in] N Degree of A matrix
 *  \return Determinant of A
 */
AMI_DLL_CPP double ami_determinante(double **A,int N)
{
  int i, k, l, cont;
  double **B,paso;
  // printf("N=%d ",N);
  if(N == 1)
    return(A[0][0]);
  ami_calloc2d(B,double,N-1,N-1);
  paso=0;
  cont=-1;
  for(i=0; i<N; i++)
  {
    cont*=-1;
    for(k=0; k<N-1; k++)
    {
      for(l=0; l<N-1; l++)
      {
        B[k][l]=A[k+1][l>=i?l+1:l];
      }
    }
    paso+=cont*A[0][i]*ami_determinante(B,N-1);
  }
  ami_free2d(B);
  return(paso);
}


/**
 * \fn double ami_determinante(A[3][3])
 *  \brief Function to evaluate the determinant of a matrix
 *  \pre Any parameter can be null.
 *  \param[in] A Matrix
 *  \param[in] N Degree of A matrix
 *  \return Determinant of A
 */
AMI_DLL_CPP double ami_determinante(double A[3][3])
{
  int i, k, l, cont;
  double **B, paso;
  int N=3;
  ami_calloc2d(B,double,N-1,N-1);
  paso=0;
  cont=-1;
  for(i=0; i<N; i++)
  {
    cont*=-1;
    for(k=0; k<N-1; k++)
    {
      for(l=0; l<N-1; l++)
      {
        B[k][l]=A[k+1][l>=i?l+1:l];
      }
    }
    paso+=cont*A[0][i]*ami_determinante(B,N-1);
  }
  ami_free2d(B);
  return(paso);
}


/**
 * \fn void ami_polynom_determinant(double p[6][6][19],int Np,int Nd,double *q)
 *  \brief  Function to compute the determinant of a polynom matrix
 *  \pre Any parameter can be null.
 *  \param[in] p
 *  \param[in] Np
 *  \param[in] Nd
 *  \param[out] q
 */
AMI_DLL_CPP void ami_polynom_determinant(double p[6][6][19],int Np,int Nd,double *q)
{
  int i,j,k,l,m,cont;
  double *q2;
  double p2[6][6][19];
  if(Nd==1)
  { 
    for(i=0; i<=18; i++)
      q[i]=p[0][0][i]; 
    return;
  }
  
  for(i=0; i<6; i++)
    for(j=0; j<6; j++)
      for(m=0; m<19; m++)
        p2[i][j][m]=0.;
    
  q2=(double*)malloc(sizeof(double)* (Np+1));
  for(i=0; i<=Np; i++)
    q[i]=0;
  cont=-1;
  
  for(i=0; i<Nd; i++)
  {
    for(k=0; k<=Np; k++) 
      q2[k]=0;
    cont*=-1;
    for(k=0; k<(Nd-1); k++)
    {
      for(l=0; l<(Nd-1); l++)
      {
        for(m=0; m<=Np; m++)
        { 
          p2[k][l][m]= p[k+1][l >= i ? l+1 : l][m];
        }
      }
    }
    ami_polynom_determinant(p2,Np,Nd-1,q2);
    
    if(cont < 0) 
      for(m=0; m<=Np; m++)
        q2[m]=-q2[m];
    
    q=ami_1v_polynom_multiplication(p[0][i],Np,q2,Np,q);
  }
  free(q2);
}


/**
 * \fn double ami_2v_polynom_evaluation(double **pt1,int N1,double x,double y)
 *  \brief  Function to evaluate a 2 variable polynom in one point
 *  \pre Any parameter can be null.
 *  \param[in] pt1 2 variable polynom
 *  \param[in] N1  Degree of polynom 1
 *  \param[in] x,y Point coordinate where the polynom will be evaluated
 *  \return Value o evaluation
 */
AMI_DLL_CPP double ami_2v_polynom_evaluation(
   double **pt1, /* 2 VARIABLE POLYNOM (INPUT)*/
   int N1,   /* DEGREE OF POLYNOM 1 (INPUT)*/
   double x,double y)  /* POINT COORDINATE WHERE THE POLYNOM WILL BE EVALUATED (INPUT) */
{
  int i, j/*,k*/;
  double *p,*q,paso;
  p=(double*)malloc(sizeof(double)*(N1+1));
  q=(double*)malloc(sizeof(double)*(N1+1));
  for(i=0;i<=N1;i++)
  {
    for(j=0; j<=N1; j++) 
      p[j]=pt1[i][j];
    q[i]=ami_polynomial_evaluation(p,N1,y);
  }
  paso=ami_polynomial_evaluation(q,N1,x);
  free(p); free(q);
  return(paso);
}


/**
 * \fn void ami_2v_polynom_to_1v_polynom(double **pt1,int N1,double *p3,double z,int flat)
 *  \brief  Function to evaluate a 2 variable polynom in one of the variable value.
 *  \pre Any parameter can be null.
 *  \post The output is a 1 degree polynom
 *  \param[in] pt1 2 variable polynom
 *  \param[in] N1  Degree of polynom 1
 *  \param[out] p3 Output 1 variable polynom
 *  \param[in] z Point where the 2 variable polynom is going to be evaluated
 *  \param[in] flat Variable where the polynom is going to be evaluated
 */
AMI_DLL_CPP void ami_2v_polynom_to_1v_polynom(
   double **pt1, /* 2 VARIABLE POLYNOM (INPUT)*/
   int N1,   /* DEGREE OF POLYNOM 1 (INPUT)*/
   double *p3, /* OUTPUT 1 VARIABLE POLYNOM (OUTPUT)*/
   double z,  /* POINT WHERE THE 2 VARIABLE POLYNOM IS GOING TO BE EVALUATED */
   int flat)  /* VARIABLE WHERE THE POLYNOM IS GOING TO BE EVALUATED */
{
  int i, j/*,k*/;
  double *p;
  p=(double*)malloc(sizeof(double)*(N1+1));
  if(flat==1){
    for(i=0; i<=N1; i++)
    {
      for(j=0; j<=N1 ;j++)
        p[j]=pt1[i][j];
      p3[i]=ami_polynomial_evaluation(p,N1,z);
    }
  }
  else
  {
    for(i=0; i<=N1; i++)
    {
      for(j=0; j<=N1; j++) 
        p[j]=pt1[j][i];
      p3[i]=ami_polynomial_evaluation(p,N1,z);
    }
  }
  free(p);
}

/**
 * \fn double* ami_1v_polynom_multiplication(double *pt1,int N1,double *pt2,int N2,double *pt3)
 *  \brief  Function to multiply polinoms of 1 variable.
 *  \pre Any parameter can be null.
 *  \post The result is added to the output polynom coeficients.
 *  \param[in] pt1 POLYNOM 1
 *  \param[in] N1  Degree of polynom 1
 *  \param[in] pt2 POLYNOM 2
 *  \param[in] N2  Degree of polynom 2
 *  \param[out] pt3 Output polynom
 *  \return p3
 */
AMI_DLL_CPP double* ami_1v_polynom_multiplication(
  double *pt1, /* POLYNOM 1 (INPUT) */
  int N1, /* DEGREE OF POLYNOM 1 (INPUT)*/
  double *pt2,/* POLYNOM 2 (INPUT) */
  int N2, /* DEGREE OF POLYNOM 2 (INPUT) */
  double *pt3) /* OUTPUT POLYNOM (INPUT-OUTPUT)*/
{
  int i,j;
  /* WE MULTIPLY THE POLYNOMS */
  for(i=0; i<=N1; i++)
  {
    if(pt1[i] != 0)
    {
      for(j=0; j<=N2; j++)
        if(pt2[j] != 0)
          pt3[i+j] += pt1[i]*pt2[j];
    }
  }
  return(pt3);
 }


/**
 * \fn void ami_2v_polynom_multiplication(double **pt1,int N1,double **pt2,int N2,double **pt3)
 *  \brief  Function to multiply polinoms of 2 variables.
 *  \pre Any parameter can be null.
 *  \post The result is added to the output polynom coeficients.
 *  \param[in] pt1 POLYNOM 1
 *  \param[in] N1  Degree of polynom 1
 *  \param[in] pt2 POLYNOM 2
 *  \param[in] N2  Degree of polynom 2
 *  \param pt3 Output polynom
 */
AMI_DLL_CPP void ami_2v_polynom_multiplication(
  double **pt1, /* POLYNOM 1 (INPUT) */
  int N1, /* DEGREE OF POLYNOM 1 (INPUT)*/
  double **pt2,/* POLYNOM 2 (INPUT) */
  int N2, /* DEGREE OF POLYNOM 2 (INPUT) */
  double **pt3) /* OUTPUT POLYNOM (INPUT - OUTPUT)*/
{
  int i, j, k, l;
  for(i=0; i<=N1; i++)
  {
    for(j=0; j<=N1; j++)
    {
      if(pt1[i][j] != 0)
      {
        for(k=0; k<=N2; k++)
          for(l=0; l<=N2; l++)
            if(pt2[k][l] != 0 )
              pt3[i+k][j+l] += pt1[i][j]*pt2[k][l];
      }
    }
  }
}

/**
 * \fn int ami_RootCubicPolynomial(double *a,int N,double *x)
 *  \brief  Function to compute the real roots of a cubic polynomial.
 *  \pre Any parameter can be null.
 *  \pre N has to be 3.
 *  \post It returns the number of roots found sorted by magnitud.
 *  \param[in] a POLINOMIAL COEFICIENTS a[0]+a[1]x+a[2]x^2 +...
 *  \param[in] N  Degree of polinomial (it has to be 3)
 *  \param[out] x Polinomial roots
 * \author Luis Alvarez
 */
AMI_DLL_CPP int ami_RootCubicPolynomial(
    double *a, /* POLINOMIAL COEFICIENTS a[0]+a[1]x+a[2]x^2 +... */
    int N, /* DEGREE OF POLINOMIAL (IT HAS TO BE 3) */
    double *x) /* POLINOMIAL ROOTS */
{
  double a1,a2,a3,Q,R,S,T,D,A;
  if(N!=3 || a[3]==0)
    return(-100000);
  a1=a[2]/a[3];
  a2=a[1]/a[3];
  a3=a[0]/a[3];
  Q=(3*a2-a1*a1)/9.;
  R=(9*a1*a2-27*a3-2*a1*a1*a1)/54.;
  D=Q*Q*Q+R*R;
  if(D > 0)
  {
    S=R+sqrt(D);
    T=R-sqrt(D);
    if(S > 0)
      S=pow(S,(double)1./3.);
    else
      S=-pow(-S,(double)1./3.);
    if(T > 0)
      T=pow(T,(double)1./3.);
    else
      T=-pow(-T,(double)1./3.);
    x[0]=S+T-a1/3.;
    return(1);
  }
  else
  {
    double PI2=acos(-1.);
    if(Q != 0)
      A=acos(R/sqrt(-Q*Q*Q));
    else
      A = 0;
    //printf("Q=%lf,R=%lf,S=%lf,T=%lf,D=%lf,A=%lf,PI=%lf\n",Q,R,S,T,D,A,PI2);
    Q=2.*sqrt(-Q);
    x[0]=Q*cos(A/3.)-a1/3.;
    x[1]=Q*cos(A/3.+2.*PI2/3.)-a1/3.;
    x[2]=Q*cos(A/3+4.*PI2/3.)-a1/3.;

    if(fabs(x[0])>fabs(x[1])){ Q=x[1]; x[1]=x[0]; x[0]=Q; }
    if(fabs(x[0])>fabs(x[2])){ Q=x[2]; x[2]=x[0]; x[0]=Q; }
    if(fabs(x[1])>fabs(x[2])){ Q=x[2]; x[2]=x[1]; x[1]=Q; }

    return(3);
  }
}

/**
 * \fn double ami_polynomial_evaluation(double *a,int Na,double x)
 *  \brief  Evaluation of a polynom using horner algorithm.
 *  \pre Any parameter can be null.
 *  \pre Na has to be positive.
 *  \post It returns the number of roots found sorted by magnitud.
 *  \param[in] a Polinomial coeficients a[0]+a[1]x+a[2]x^2 +...
 *  \param[in] Na Polynom degree
 *  \param[in] x Point where the polynom is evaluated
 *  \return Evaluation of x in polynomial a
 * \author Luis Alvarez
 */
AMI_DLL_CPP double ami_polynomial_evaluation(
  const double *a, /* POLYNOM COEFICIENT */
  int Na, /* POLYNOM DEGREE */
  double x) /* POINT WHERE THE POLYNOM IS EVALUATED */
{
  double sol = a[Na];
  for(int i = Na - 1; i > -1; i--)
    sol = sol * x + a[i];
  return sol;
}

/**
 * \fn int ami_lens_distortion_polynomial_update(double *x, double *y,int Np,double *a,int Na,double x0,double y0,int k,double *pol)
 *  \brief  Function to add the information of a line point sequence to the 4 degree polynomial to compute the lens distortion model.
 *  \pre Any parameter can be null.
 *  \pre Np and Na have to be positive.
 *  \param[in] x,y Distorted line coordinates
 *  \param[in] Np Number of points
 *  \param[in] a Polynomial defining the lens distortion model
 *  \param[in] Na Degree of polynomial model for lens distortion
 *  \param[in] x0,y0 Coordinates of the image center
 *  \param[in] k Coeficient of the polynomial to be updated
 *  \param pol 4 degree polynom to minimize
 *  \return 0
 */
AMI_DLL_CPP int ami_lens_distortion_polynomial_update(
   double *x, double *y,  /* DISTORTED LINE COORDINATES (INPUT)*/
   int Np, /* NUMBER OF POINTS (INPUT)*/
   double *a, /* POLYNOMIAL DEFINING THE LENS DISTORTION MODEL (INPUT) */
   int Na,   /* DEGREE OF POLYNOMIAL MODEL FOR LENS DISTORTION */
   double x0,double y0, /* COORDINATES OF THE IMAGE CENTER */
   int k,  /* COEFICIENT OF THE POLYNOMIAL TO BE UPDATED*/
   double *pol) /* 4 DEGREE POLYNOM TO MINIMIZE (INPUT-OUTPUT) */
{
  int i, j;
  double *A, *x2, *y2, *d, x2_m, y2_m, s_xx, s_yy/*,x_m,y_m,s_max*/, xA_m, xd_m, yA_m, yd_m;
  double pol1[5],pol2[5],pol3[5];
   if(Np < 3)
     return(-1);
  /* WE ALLOCATE MEMORY */
  A=(double*)malloc( sizeof(double)*Np );
  x2=(double*)malloc( sizeof(double)*Np );
  y2=(double*)malloc( sizeof(double)*Np );
  d=(double*)malloc( sizeof(double)*Np );

  /* WE COMPUTE THE DISTANCE TO THE IMAGE CENTER */
  for(i=0; i<Np; i++)
    d[i]=sqrt( (x[i]-x0)*(x[i]-x0)+(y[i]-y0)*(y[i]-y0) );

  /* WE COMPUTE THE POINT TRANSFORMATION WITH THE CURRENT LENS DISTORTION MODEL */
  for(i=0; i<Np; i++){
    A[i]=ami_polynomial_evaluation(a,Na,d[i]);
    x2[i]=x0+(x[i]-x0)*A[i];
    y2[i]=y0+(y[i]-y0)*A[i];
  }

  /* WE COMPUTE THE DISTANCE POWER k (THE COEFICIENT OF THE LENS DISTORTION MODEL
    TO BE UPDATED */
  for(i=0; i<Np; i++)
    d[i]=pow(d[i],(double) k);

  /* WE COMPUTE THE VARIANCE OF THE TRANSFORMED POINTS */
  x2_m=0; 
  for(i=0; i<Np; i++) 
    x2_m += x2[i];
  
  x2_m/=Np;
  s_xx=0; 
  for(i=0;i<Np;i++) 
    s_xx += (x2[i]-x2_m)*(x2[i]-x2_m); 
  
  s_xx/=Np;
  y2_m=0; 
  for(i=0; i<Np; i++)
    y2_m+=y2[i];
  
  y2_m/=Np;
  s_yy=0; 
  for(i=0; i<Np; i++) 
    s_yy += (y2[i]-y2_m)*(y2[i]-y2_m); 
  
  s_yy/=Np;
//  s_max=s_xx>s_yy?s_xx:s_yy;

  /* WE COMPUTE SOME AVERAGES WE NEED */
  xA_m=0; for(i=0;i<Np;i++) xA_m+=(x[i]-x0)*A[i]; xA_m/=Np;
  xd_m=0; for(i=0;i<Np;i++) xd_m+=(x[i]-x0)*d[i]; xd_m/=Np;
  yA_m=0; for(i=0;i<Np;i++) yA_m+=(y[i]-y0)*A[i]; yA_m/=Np;
  yd_m=0; for(i=0;i<Np;i++) yd_m+=(y[i]-y0)*d[i]; yd_m/=Np;

  /* WE COMPUTE THE POLYNOMIAL TO MINIMIZE */
  for(i=0; i<5; i++)
    pol1[i]=pol2[i]=pol3[i]=0;
  for(i=0; i<Np; i++)
  {
    pol1[0]+=((x[i]-x0)*A[i]-xA_m)*((x[i]-x0)*A[i]-xA_m);
    pol1[1]+=2.*((x[i]-x0)*A[i]-xA_m)*((x[i]-x0)*d[i]-xd_m);
    pol1[2]+=((x[i]-x0)*d[i]-xd_m)*((x[i]-x0)*d[i]-xd_m);
    pol2[0]+=((y[i]-y0)*A[i]-yA_m)*((y[i]-y0)*A[i]-yA_m);
    pol2[1]+=2.*((y[i]-y0)*A[i]-yA_m)*((y[i]-y0)*d[i]-yd_m);
    pol2[2]+=((y[i]-y0)*d[i]-yd_m)*((y[i]-y0)*d[i]-yd_m);
    pol3[0]+=((y[i]-y0)*A[i]-yA_m)*((x[i]-x0)*A[i]-xA_m);
    pol3[1]+=((y[i]-y0)*A[i]-yA_m)*((x[i]-x0)*d[i]-xd_m)+
             ((y[i]-y0)*d[i]-yd_m)*((x[i]-x0)*A[i]-xA_m);
    pol3[2]+=((y[i]-y0)*d[i]-yd_m)*((x[i]-x0)*d[i]-xd_m);
  }
  for(i=0; i<3; i++)
  {
    for(j=0; j<3; j++)
    {
      pol[i+j] += (pol1[i]*pol2[j]-pol3[i]*pol3[j])/1. /*s_max*/;
    }
  }
  /* WE FREE MEMORY */
  free(A);
  free(x2);
  free(y2);
  free(d);

  return(0);
}

/**
 * \fn int ami_lens_distortion_model_update(double *a,int Na,int k,double *pol)
 *  \brief  Function to update the lens distortion model by minimizing a 4 degree polynom.
 *  \pre Any parameter can be null.
 *  \pre Na has to be positive.
 *  \param a Polynomial defining the lens distortion model
 *  \param[in] Na Degree of polynomial model for lens distortion
 *  \param[in] k Coeficient of the polynomial to be updated
 *  \param pol 4 degree polynom to minimize
 *  \return 0
 */
AMI_DLL_CPP int ami_lens_distortion_model_update(
   double *a, /* POLYNOMIAL DEFINING THE LENS DISTORTION MODEL (INPUT-OUTPUT) */
   int k,  /* COEFICIENT OF THE POLYNOMIAL TO BE UPDATED (INPUT)*/
   double *pol) /* 4 DEGREE POLYNOM TO MINIMIZE (INPUT) */
{
  int j, i, M;
  double *x, *b, p[3];
  x = (double*)malloc( sizeof(double)*3 );
  b = (double*)malloc( sizeof(double)*4 );
  b[0]=pol[1];
  b[1]=2*pol[2];
  b[2]=3.*pol[3];
  b[3]=4.*pol[4];
  M = ami_RootCubicPolynomial(b,3,x);
  for(i=0; i<M; i++)
    p[i]=ami_polynomial_evaluation(pol,4,x[i]);
  
  j = 0;
  if(M == 3)
  {
    if(p[j]>p[1] /*&& fabs(x[1])<1. */) j=1;
    if(p[j]>p[2] /*&& fabs(x[2])<1.*/ ) j=2;
  }

  //if(fabs(x[j])<1.) a[k]+=x[j];

  a[k] += x[j];

  //if(k==1 && a[1]<0) a[1]=0;
  //if(k==Na && a[k]<0) a[k]=0;
  free(x);
  free(b);
  return(0);
}

/**
 * \fn double ami_LensDistortionVarianceError(double *x,double *y,int Np,double x0,double y0,double *a,int Na)
 *  \brief  Function to compute the lens distortion energy error (the residual variance of the point distribution.
 *  \pre Any parameter can be null.
 *  \pre Np and Na have to be positive.
 *  \param[in] x,y Original point distribution
 *  \param[in] Np Number of points
 *  \param[in] x0,y0 Center of the image
 *  \param[in] a Lens Distortion Polynomial model
 *  \param[in] Na Degree of polynomial model
 *  \return The lens distortion energy error
 * \author Luis Alvarez
 */
AMI_DLL_CPP double ami_LensDistortionVarianceError(
  double *x,double *y, /* ORIGINAL POINT DISTRIBUTION (INPUT)*/
  int Np, /* NUMBER OF POINTS (INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *a, /* Lens Distortion Polynomial model (INPUT)*/
  int Na) /* Degree of Polynomial model (INPUT)*/
{
int i/*,j*/;
  double A,*x2, *y2, d, x2_m, y2_m, s_xx, s_yy, s_max, s_xy/*,y_m,x_m,xA_m,xd_m,yA_m,yd_m*/;
   if(Np < 3)
     return((double) 0.);
  
  /* WE ALLOCATE MEMORY */
  x2=(double*)malloc( sizeof(double)*Np );
  y2=(double*)malloc( sizeof(double)*Np );

  /* WE COMPUTE THE POINT TRANSFORMATION USING THE LENS DISTORTION MODEL*/
  for(i=0; i<Np; i++)
  {
    d=sqrt( (x[i]-x0)*(x[i]-x0)+(y[i]-y0)*(y[i]-y0) );
    A=ami_polynomial_evaluation(a,Na,d);
    x2[i]=x0+(x[i]-x0)*A;
    y2[i]=y0+(y[i]-y0)*A;
  }
  /* WE COMPUTE THE VARIANCE OF THE TRANSFORMED POINTS */
  x2_m=0; 
  for(i=0; i<Np; i++)
    x2_m += x2[i]; 
  x2_m/=Np;
  
  s_xx=0; 
  for(i=0; i<Np; i++) 
    s_xx += (x2[i]-x2_m)*(x2[i]-x2_m); 
  s_xx/=Np;
  y2_m=0; 
  for(i=0; i<Np; i++) 
    y2_m+=y2[i];  
  y2_m/=Np;
  
  s_yy=0;
  for(i=0; i<Np; i++)
    s_yy+=(y2[i]-y2_m)*(y2[i]-y2_m);
  s_yy/=Np;
  s_xy=0; 
  for(i=0;i<Np;i++) 
    s_xy+=(y2[i]-y2_m)*(x2[i]-x2_m); 
  s_xy/=Np;
  s_max=s_xx>s_yy ? s_xx : s_yy;

  /* WE FREE MEMORY */
  free(x2); free(y2);
//printf("Distortion Error =%e\n",(s_xx*s_yy-s_xy*s_xy)/s_max);
  return((s_xx*s_yy-s_xy*s_xy)/s_max);

}


/**
 * \fn double ami_LensDistortionEnergyError(double *x,double *y,int Np,double x0,double y0,double *a,int Na)
 *  \brief  Function to compute the lens distortion energy error.
 *  \pre Any parameter can be null.
 *  \pre Np and Na have to be positive.
 *  \param[in] x,y Original point distribution
 *  \param[in] Np Number of points
 *  \param[in] x0,y0 Center of the image
 *  \param[in] a Lens Distortion Polynomial model
 *  \param[in] Na Degree of polynomial model
 *  \return The lens distortion energy error
 * \author Luis Alvarez
 */
AMI_DLL_CPP double ami_LensDistortionEnergyError(
  double *x,double *y, /* ORIGINAL POINT DISTRIBUTION (INPUT)*/
  int Np, /* NUMBER OF POINTS (INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *a, /* Lens Distortion Polynomial model (INPUT)*/
  int Na) /* Degree of Polynomial model (INPUT)*/
{
int i/*,j*/;
  double A,*x2,*y2,d,x2_m,y2_m,s_xx,s_yy,s_xy/*,y_m,s_max,x_m,xA_m,xd_m,yA_m,yd_m*/;
  /* WE ALLOCATE MEMORY */
   if(Np < 3)
     return((double) 0.);
  x2=(double*)malloc( sizeof(double)*Np );
  y2=(double*)malloc( sizeof(double)*Np );

  /* WE COMPUTE THE POINT TRANSFORMATION USING THE LENS DISTORTION MODEL*/
  for(i=0; i<Np; i++)
  {
    d=sqrt( (x[i]-x0)*(x[i]-x0)+(y[i]-y0)*(y[i]-y0) );
    A=ami_polynomial_evaluation(a,Na,d);
    x2[i]=x0+(x[i]-x0)*A;
    y2[i]=y0+(y[i]-y0)*A;
  }
  /* WE COMPUTE THE ENERGY OF THE TRANSFORMED POINTS */
  x2_m=0;
  for(i=0; i<Np; i++)
    x2_m+=x2[i];
  x2_m/=Np;
  
  s_xx=0;
  for(i=0; i<Np; i++) 
    s_xx+=(x2[i]-x2_m)*(x2[i]-x2_m); 
  s_xx/=Np;
  
  y2_m=0;
  for(i=0; i<Np; i++)
    y2_m += y2[i];
  y2_m/=Np;
  
  s_yy=0; 
  for(i=0; i<Np; i++)
    s_yy += (y2[i]-y2_m)*(y2[i]-y2_m);
  s_yy/=Np;
  
  s_xy=0; 
  for(i=0; i<Np; i++)
    s_xy+=(y2[i]-y2_m)*(x2[i]-x2_m);
  s_xy/=Np;

  /* WE FREE MEMORY */
  free(x2); 
  free(y2);
//printf("Distortion Error =%e\n",(s_xx*s_yy-s_xy*s_xy));
  return((s_xx*s_yy-s_xy*s_xy));

}

/**
 * \fn double ami_points_to_line_equation(double *a, int Na, double xc,double yc, double *x,double *y,int N, double line[3])
 *  \brief  COMPUTATION OF THE LINE EQUATION INCLUDING THE POLINOMIAL DISTORTION MODEL
 *  \pre Any parameter can be null.
 *  \param[in]  a INPUT POLINOMIAL DISTORTION MODEL
 *  \param[in]  Na INPUT POLINOMIAL DISTORTION MODEL DEGREE
 *  \param[in]  xc, yc INPUT CENTER OF THE DISTORTION MODEL
 *  \param[in]  x, y INPUT POINTS COORDINATES
 *  \param[in]  N INPUT NUMBER OF POINTS
 *  \param[out] line OUTPUT LINE EQUATION
 *  return AVERAGE ERROR
 * \author Luis Alvarez
 */
AMI_DLL_CPP double ami_points_to_line_equation(
  double *a, // INPUT POLINOMIAL DISTORTION MODEL
  int Na, // INPUT POLINOMIAL DISTORTION MODEL DEGREE
  double xc,double yc, // INPUT CENTER OF THE DISTORTION MODEL
  double *x,double *y, //INPUT POINTS COORDINATES
  int N,  // INPUT NUMBER OF POINTS
  double line[3]) // OUTPUT LINE EQUATION
{
  int i, j, k;
  double suu, suv, svv/*,h0,h1*/, um, vm, h, r[4][3], min, paso, norma;
  double cero=1e-20;
  double *x2, *y2;
  double error=0;
  if(N < 2)
  {
    printf("Numero de puntos para el Calculo de la recta 2D menor que 2\n");
    return((double) -1.);
  }
  if(Na > 0)
  {
    x2=(double*)malloc(sizeof(double)*N);
    y2=(double*)malloc(sizeof(double)*N);
  }
  else
  {
    x2=x;
    y2=y;
  }
  suu=0; 
  suv=0; 
  svv=0; 
  um=0; 
  vm=0;
  
  for(i=0; i<N; i++)
  {
  ami_lens_distortion_model_evaluation(a,Na,xc,yc,x[i],y[i],&(x2[i]),&(y2[i]));
    um+=x2[i];
    vm+=y2[i];
  }
  um/=N;
  vm/=N;
  
  for(i=0; i<N; i++)
  {
    suu+=(x2[i]-um)*(x2[i]-um);
    svv+=(y2[i]-vm)*(y2[i]-vm);
    suv+=(x2[i]-um)*(y2[i]-vm);
  }
  suu/=N; 
  svv/=N; 
  suv/=N;
  if(fabs(suv)<= cero)
  {
    if(suu<svv && svv>cero)
    {
      line[0]=1; 
      line[1]=0; 
      line[2]=-um;
      for(i=0; i<N; i++)
        error+=fabs(line[0]*x2[i]+line[1]*y2[i]+line[2]);
      if(Na > 0)
      {
        free(x2);
        free(y2);
      }
      return(error/N);
    }
    if(svv<suu && suu>cero)
    {
      line[0]=0;
      line[1]=1;
      line[2]=-vm;
      for(i=0; i<N; i++)
        error += fabs(line[0]*x2[i]+line[1]*y2[i]+line[2]);
      if(Na > 0)
      {
        free(x2);
        free(y2);
      }
      return(error/N);
    }
    printf("No se pudo calcular la recta 2D\n");
    if(Na > 0)
    {
      free(x2);
      free(y2);
    }
    return((double)-1.);
  }

  r[2][1]=r[3][1]=r[0][0]=r[1][0]=1.;
  h=0.5*(suu-svv)/suv;
  if(h > 0)
  {
    r[0][1]=-h-sqrt(1.+h*h);
    r[0][2]=-(um+r[0][1]*vm);
    r[1][1]=-1./r[0][1];
    r[1][2]=-(um+r[1][1]*vm);
    r[2][0]=h+sqrt(1.+h*h);
    r[2][2]=-(r[2][0]*um+vm);
    r[3][0]=-1./r[2][0];
    r[3][2]=-(r[3][0]*um+vm);
  }
  else
  {
    r[0][1]=-h+sqrt(1+h*h);
    r[0][2]=-(um+r[0][1]*vm);
    r[1][1]=-1./r[0][1];
    r[1][2]=-(um+r[1][1]*vm);
    r[2][0]=h-sqrt(1+h*h);
    r[2][2]=-(r[2][0]*um+vm);
    r[3][0]=-1./r[2][0];
    r[3][2]=-(r[3][0]*um+vm);
  }

  for(j=0; j<4; j++)
  {
    norma=sqrt(r[j][0]*r[j][0]+r[j][1]*r[j][1]);
    for(i=0; i<3; i++)
      r[j][i]/=norma;
  }

  min=0.; k=0;
  for(i=0; i<N; i++)
  {
   paso=r[0][0]*x2[i]+r[0][1]*y2[i]+r[0][2];
   min+=paso*paso;
  }
  for(j=1; j<4; j++)
  {
    h=0;
    for(i=0; i<N; i++)
    {
      paso=r[j][0]*x2[i]+r[j][1]*y2[i]+r[j][2];
      h+=paso*paso;
    }
    if(h < min)
    {
      k=j;
      min=h;
    }
  }

  line[0]=r[k][0];
  line[1]=r[k][1];
  line[2]=r[k][2];
  for(i=0; i<N; i++)
    error += fabs(line[0]*x2[i]+line[1]*y2[i]+line[2]);
  if(Na > 0)
  {
    free(x2);
    free(y2);
  }
  return(error/N);

}

/**
 * \fn double ami_points_to_line_equation_outlier_elimination(double *a, int Na, double xc,double yc, double *x,double *y,int *Np, double line[3],double outlier_elimination_factor)
 *  \brief  CALCULA LA RECTA QUE MEJOR APROXIMA UN CONJUNTO DE PUNTOS 2D  Y ELIMINA PUNTOS ALEJADOS
  TENIENDO EN CUENTA EL MODELO DE DISTORSION
 *  \pre Any parameter can be null.
 *  \param [in] a INPUT POLINOMIAL DISTORTION MODEL
 *  \param [in] Na INPUT POLINOMIAL DISTORTION MODEL DEGREE
 *  \param [in]  xc, yc INPUT CENTER OF THE DISTORTION MODEL
 *  \param [in]  x, y INPUT POINTS COORDINATES
 *  \param Np  INPUT/OUTPUT NUMBER OF POINTS OF THE LINE
 *  \param [out]  line OUTPUT LINE EQUATION
 *  \param [in] outlier_elimination_factor INPUT WE ELIMINATE POINTS SUCH THAT THE DISTANCE TO THE LINE BE BIGGER THAN THE AVERAGE DISTANCE + outlier_elimination_factor*STANDARD DEVIATION
 *  return AVERAGE ERROR
 * \author Luis Alvarez
 */
AMI_DLL_CPP double ami_points_to_line_equation_outlier_elimination(
  double *a, // INPUT POLINOMIAL DISTORTION MODEL
  int Na, // INPUT POLINOMIAL DISTORTION MODEL DEGREE
  double xc,double yc, // INPUT CENTER OF THE DISTORTION MODEL
  double *x,double *y, //INPUT POINTS COORDINATES
  int *Np,  // INPUT/OUTPUT NUMBER OF POINTS
  double line[3], // OUTPUT LINE EQUATION
  double outlier_elimination_factor)  // INPUT WE ELIMINATE POINTS SUCH THAT THE DISTANCE TO THE
                                    //LINE BE BIGGER THAN THE AVERAGE DISTANCE +
                                    // + outlier_elimination_factor*STANDARD DEVIATION
{
  int i, j, k;
  double suu, suv, svv/*,h0,h1*/, um, vm, h, r[4][3], min, paso, norma;
  double cero=1e-20;
  double *x2, *y2;
  double error=0;
  int N = *Np;
  if(N < 2)
  {
    printf("Numero de puntos para el Calculo de la recta 2D menor que 2\n");
    return((double) -1.);
  }

   if(Na > 0)
   {
    x2=(double*)malloc(sizeof(double)*N);
    y2=(double*)malloc(sizeof(double)*N);
  }
  else
  {
    x2=x;
    y2=y;
  }
  suu=0;
  suv=0;
  svv=0;
  um=0;
  vm=0;
  
  for(i=0; i<N; i++)
  {
  ami_lens_distortion_model_evaluation(a,Na,xc,yc,x[i],y[i],&(x2[i]),&(y2[i]));
    um+=x2[i];
    vm+=y2[i];
  }
  um/=N;
  vm/=N;
  for(i=0; i<N; i++)
  {
    suu+=(x2[i]-um)*(x2[i]-um);
    svv+=(y2[i]-vm)*(y2[i]-vm);
    suv+=(x2[i]-um)*(y2[i]-vm);
  }
  suu/=N;
  svv/=N;
  suv/=N;
  if(fabs(suv)<= cero)
  {
    if(suu<svv && svv>cero)
    {
      line[0]=1; 
      line[1]=0; 
      line[2]=-um;
      for(i=0; i<N; i++)
        error += fabs(line[0]*x2[i]+line[1]*y2[i]+line[2]);
      if(Na > 0)
      {
        free(x2);
        free(y2);
      }
      return(error/N);
    }
    if(svv<suu && suu>cero)
    {
      line[0]=0;
      line[1]=1;
      line[2]=-vm;
      for(i=0; i<N; i++)
        error+=fabs(line[0]*x2[i]+line[1]*y2[i]+line[2]);
      if(Na > 0)
      {
        free(x2);
      free(y2);
      }
      return(error/N);
    }
    printf("No se pudo calcular la recta 2D\n");
    if(Na > 0)
    {
      free(x2);
      free(y2);
    }
    return((double)-1.);
  }

  r[2][1]=r[3][1]=r[0][0]=r[1][0]=1.;
  h=0.5*(suu-svv)/suv;
  if(h > 0)
  {
    r[0][1]=-h-sqrt(1.+h*h);
    r[0][2]=-(um+r[0][1]*vm);
    r[1][1]=-1./r[0][1];
    r[1][2]=-(um+r[1][1]*vm);
    r[2][0]=h+sqrt(1.+h*h);
    r[2][2]=-(r[2][0]*um+vm);
    r[3][0]=-1./r[2][0];
    r[3][2]=-(r[3][0]*um+vm);
  }
  else
  {
    r[0][1]=-h+sqrt(1+h*h);
    r[0][2]=-(um+r[0][1]*vm);
    r[1][1]=-1./r[0][1];
    r[1][2]=-(um+r[1][1]*vm);
    r[2][0]=h-sqrt(1+h*h);
    r[2][2]=-(r[2][0]*um+vm);
    r[3][0]=-1./r[2][0];
    r[3][2]=-(r[3][0]*um+vm);
  }

  for(j=0; j<4; j++)
  {
    norma=sqrt(r[j][0]*r[j][0]+r[j][1]*r[j][1]);
    for(i=0;i<3;i++)
      r[j][i]/=norma;
  }

  min=0.; k=0;
  for(i=0; i<N; i++)
  {
   paso=r[0][0]*x2[i]+r[0][1]*y2[i]+r[0][2];
   min+=paso*paso;
  }
  for(j=1; j<4; j++)
  {
    h=0;
    for(i=0; i<N; i++)
    {
      paso=r[j][0]*x2[i]+r[j][1]*y2[i]+r[j][2];
      h+=paso*paso;
    }
    if(h < min)
    {
      k=j;
      min=h;
    }
  }

  line[0]=r[k][0]; 
  line[1]=r[k][1]; 
  line[2]=r[k][2];
  // OUTLIER ELIMINATION
  for(i=0; i<N; i++)
    error += fabs(line[0]*x2[i]+line[1]*y2[i]+line[2]);
  double average=error/N;
  double standard_deviation=0;
  for(i=0; i<N; i++)
  {
    paso=fabs(line[0]*x2[i]+line[1]*y2[i]+line[2])-average;
    standard_deviation += paso*paso;
  }
  standard_deviation=sqrt(standard_deviation/N);
  for(i=0;i<N;i++)
  {
    paso=fabs(line[0]*x2[i]+line[1]*y2[i]+line[2]);
    if(paso > (average+standard_deviation*outlier_elimination_factor))
    {
      x[i]=x[N-1];
      y[i]=y[N-1];
      x2[i]=x2[N-1];
      y2[i]=y2[N-1];
      i--;
      N--;
    }
  }
  *Np=N;
  if(Na > 0)
  {
    free(x2);
  free(y2);
  }
  return(ami_points_to_line_equation(a,Na,xc,yc,x,y,N,line));

}






/**
 * \fn double ami_distortion_model_estimation_2p(double xc2,double yc2,double **x,double **y,int Nl,int *Np,double **a, int *Na)
 *  \brief  DISTORTION MODEL ESTIMATION FROM A SET OF LINES
 *  \pre Any parameter can be null.
 *  \param [in] xc2, yc2  DISTORTION CENTER (THE CENTER OF THE IMAGE)
 *  \param [in]  x, y LINE POINT INFORMATION
 *  \param [in]  Nl NUMBER OF LINES
 *  \param [in]  Np INPUT NUMBER OF POINTS IN EACH LINE
 *  \param [out]  a OUTPUT POLINOMIAL DISTORTION MODEL
 *  \param [out]  Na OUTPUT DEGREE OF POLINOMIAL DISTORTION MODEL (Na=4)
 *  return AVERAGE ERROR
 * \author Luis Alvarez
 */

AMI_DLL_CPP double ami_distortion_model_estimation_2p(
  double xc2,double yc2, //INPUT DISTORTION CENTER (THE CENTER OF THE IMAGE)
  double **x,double **y, // INPUT LINE POINT INFORMATION
  int Nl,                // INPUT NUMBER OF LINES
  int *Np,               // INPUT NUMBER OF POINTS IN EACH LINE
  double **a, //OUTPUT POLINOMIAL DISTORTION MODEL
  int *Na) // OUTPUT DEGREE OF POLINOMIAL DISTORTION MODEL (Na=4)
{
  double **xx,**yy;
  double factor_n;
  int i,m,cont;
  double xc=xc2;
  double yc=yc2;
  xx=(double**) malloc(sizeof(double*)*Nl);
  yy=(double**) malloc(sizeof(double*)*Nl);
  *Na=4;
  *a=(double*) malloc(sizeof(double)*(*Na+1));
  (*a)[0]=1.; for(i=1;i<=(*Na);i++) (*a)[i]=0.;

  // WE COMPUTE THE DISTORSION ERROR IN THE ORIGINAL POINTS
  double Emin=0;
  int cont2=0;
  for(m=0; m<Nl; m++)
  {
    //Emin+=sqrt(ami_LensDistortionVarianceError(x[m],y[m],Np[m],xc,yc,(*a),*Na));
    Emin += ami_LensDistortionEnergyError(x[m],y[m],Np[m],xc,yc,(*a),*Na);
    cont2 += Np[m];
  }
  /* NORMALIZAMOS LAS COORDENADAS DE LOS PUNTOS */
  factor_n=0.; cont=0;
  for(m=0; m<Nl; m++)
  {
    xx[m]=(double*)malloc(sizeof(double)*Np[m]);
    yy[m]=(double*)malloc(sizeof(double)*Np[m]);
    for(i=0;i<Np[m];i++){
      xx[m][i]=x[m][i]-xc;
      yy[m][i]=y[m][i]-yc;
      cont++;
      factor_n += xx[m][i] * xx[m][i]+yy[m][i] * yy[m][i];
     }
  }

  factor_n=sqrt(factor_n/cont);
  for(m=0; m<Nl; m++)
    for(i=0; i<Np[m]; i++)
    { 
      xx[m][i]/=factor_n;
      yy[m][i]/=factor_n; 
    }

  // WE COMPUTE THE DISTORTION ERROR FOR THE NORMALIZED POINT COORDINATES
  Emin=0;
  for(m=0; m<Nl; m++)
  {
    //Emin+=sqrt(ami_LensDistortionVarianceError(x[m],y[m],Np[m],xc,yc,(*a),*Na));
    Emin += ami_LensDistortionEnergyError(xx[m],yy[m],Np[m],0.,0.,(*a),*Na);
  }

  /* WE RESET THE POLYNOMIAL LENS DISTORTION MODEL */
  (*a)[0]=1; for(i=1;i<=*Na;i++) (*a)[i]=0;
  /* WE COMPUTE THE DISTORSION MODEL */
   /*double residual_variance=*/ami_lens_distortion_estimation_2v(xx,yy,Nl,Np,(double) 0.,(double) 0.,*a,*Na,2,4,(double) 0.);
  //printf("Residual Variance =%e with Na=%d\n",residual_variance,*Na); system("PAUSE");
  /* LLAMO A LA NORMALIZACION DEL ZOOM */
  // ami_lens_distortion_zoom_normalization(xx,yy,Nl,Np,(double) 0.,(double) 0.,a,Na);

  // WE COMPUTE THE DISTORTION ERROR FOR THE NORMALIZED POINT COORDINATES
  Emin=0;
  for(m=0; m<Nl; m++)
  {
    //Emin+=sqrt(ami_LensDistortionVarianceError(x[m],y[m],Np[m],xc,yc,(*a),*Na));
    Emin += ami_LensDistortionEnergyError(xx[m],yy[m],Np[m],0.,0.,(*a),*Na);
  }
  /* WE REMOVE POINT NORMALIZATION FROM THE DISTORSION PARAMETERS */
  double paso=1.;
  for(i=0; i<=*Na; i++)
  {
    (*a)[i]=(*a)[i]*paso;
    paso/=factor_n;
  }
  // ami_printf1d("a",(*a),*Na+1);
  // WE COMPUTE THE DISTORSION ERROR IN THE ORIGINAL POINTS
  Emin=0;
  for(m=0; m<Nl; m++)
  {
    //Emin+=sqrt(ami_LensDistortionVarianceError(x[m],y[m],Np[m],xc,yc,(*a),*Na));
    Emin += ami_LensDistortionEnergyError(x[m],y[m],Np[m],xc,yc,(*a),*Na);
  }
  // WE FREE THE MEMORY
  for(m=0; m<Nl; m++)
  {
    free(xx[m]);
    free(yy[m]);
  }
  free(xx);
  free(yy);
  return(Emin/Nl);

}


/**
 * \fn void ami_lens_distortion_model_evaluation(double *a,int Na, double xc,double yc,double x_input,double y_input,double *x_output,double *y_output)
 *  \brief  COMPUTE THE LENS DISTORTION MODEL IN A POINT
 *  \pre Any parameter can be null.
 *  \pre Na has to be positive.
 *  \param [in] a INPUT POLINOMIAL DISTORTION MODEL
 *  \param [in] Na INPUT DEGREE OF POLINOMIAL DISTORTION MODEL
 *  \param [in] xc,yc INPUT CENTER OF DISTORTION
 *  \param [in] x_input,y_input INPUT POINT
 *  \param [out] x_output,y_output  OUTPUT UNDISTORTED POINT
 * \author Luis Alvarez
 */
AMI_DLL_CPP void ami_lens_distortion_model_evaluation(
  double *a, // INPUT POLINOMIAL DISTORTION MODEL
  int Na, // INPUT DEGREE OF POLINOMIAL DISTORTION MODEL
  double xc,double yc,  // INPUT CENTER OF DISTORTION
  double x_input,double y_input,  // INPUT POINT
  double *x_output,double *y_output  // OUTPUT UNDISTORTED POINT
)
{
  double norm = sqrt((x_input-xc)*(x_input-xc)+(y_input-yc)*(y_input-yc));
  double A = ami_polynomial_evaluation(a,Na,norm);
  *x_output=xc+(x_input-xc)*A;
  *y_output=yc+(y_input-yc)*A;
}



/**************************************************************************************
FUNCION AUXILIAR ami_inverse_lens_distortion()
**************************************************************************************/
AMI_DLL_CPP double ami_inverse_lens_distortion_newton_raphson(
  double x,double y, /* POINT TO INVERSE (INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *xt,double *yt, /* INVERVE POINT TRANSFORMED (OUTPUT) */
  double *a, /* LENS DISTORTION MODEL POLYNOM */
  int Na) /* DEGREE OF THE LENS DISTORTION MODEL POLYNOM */
{
  int i;
  double paso, d, *b, *b2, root;
  if(a[Na] == 0.) 
    return(-1);
  /* WE ALLOCATE MEMORY */
  b=(double*)malloc( sizeof(double)*(Na+2) );
  b2=(double*)malloc( sizeof(double)*(Na+2) );
  d=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
  b[0]=-1;
  b[1]=a[0];//1.;
  paso=d;
  for(i=2; i<(Na+2); i++)
  {
    b[i]=a[i-1]*paso;
    paso*=d;
  }
  for(i=0; i<(Na+2); i++)
    b2[i]=b[Na+1-i];

  // WE IMPROVE SOLUTION USING NEWTON RAPHSON
  double norm2;
  //double xn=*xt,yn=*yt,xp,yp,xn2,yn2,xn3,yn3,error,tol=1e-8;
  double xn=x, yn=y, xp, yp, xn2=0., yn2=0., xn3, yn3, error, tol=1e-8;
  ami_lens_distortion_model_evaluation(a, Na, x0, y0, xn, yn, &xp, &yp);
  norm2=(xp-x)*(xp-x)+(yp-y)*(yp-y);
  error=tol+1;
  int iter=0;
  root=1;
  double lambda=1;
  while(error>tol && ++iter<100)
  {
    xn2=x0+(xn-x0)*(root+1e-6);
    yn2=y0+(yn-y0)*(root+1e-6);
    xn3=x0+(xn-x0)*(root-1e-6);
    yn3=y0+(yn-y0)*(root-1e-6);
    ami_lens_distortion_model_evaluation(a,Na,x0,y0,xn2,yn2,&xp,&yp);
    double norm3=(xp-x)*(xp-x)+(yp-y)*(yp-y);
    ami_lens_distortion_model_evaluation(a,Na,x0,y0,xn3,yn3,&xp,&yp);
    double norm4=(xp-x)*(xp-x)+(yp-y)*(yp-y);
    double derivative=(norm3-norm4)/2e-6;
    double derivative2=(norm3+norm4-2.*norm2)/1e-12;
    if(derivative2 == 0)
      break;
    double root2=root-lambda*derivative/derivative2;
    error=fabs(root2-root);
    xn2=x0+(xn-x0)*root2;
    yn2=y0+(yn-y0)*root2;
    ami_lens_distortion_model_evaluation(a,Na,x0,y0,xn2,yn2,&xp,&yp);
    norm3=(xp-x)*(xp-x)+(yp-y)*(yp-y);
    if(norm3 < norm2)
    {
      root=root2;
      norm2=norm3;
      if(lambda < 1.)
        lambda*=2.;
    }
    else 
    {
      lambda/=2.;
    }
  }
  *xt=xn2;
  *yt=yn2;

  free(b); 
  free(b2);
  //free(rx); free(ry);  OJO comentado Pedro, a veces da segmentation fault
  return(norm2);
}

/**
 * \fn int ami_inverse_lens_distortion(double x,double y,double x0,double y0,double *xt,double *yt,double *a,int Na)
 *  \brief  Function to inverse the lens distortion transformation.
 *  \pre Any parameter can be null.
 *  \pre Na has to be positive.
 *  \param [in] x,y Point to inverse
 *  \param [in] x0,y0 Center of the image
 *  \param [out] xt,yt Inverve point transformed
 *  \param [in] a Lens distortion model polynom
 *  \param [in] Na Degree of the lens distortion model polynom
 *  \return 0
 * \author Luis Alvarez
 */
AMI_DLL_CPP int ami_inverse_lens_distortion(
  double x,double y, /* POINT TO INVERSE (INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *xt,double *yt, /* INVERVE POINT TRANSFORMED (OUTPUT) */
  double *a, /* LENS DISTORTION MODEL POLYNOM */
  int Na) /* DEGREE OF THE LENS DISTORTION MODEL POLYNOM */
{
  int i, Nr;
  double paso, d, *b, *b2, *rx, *ry, root;
  if(a[Na]==0.) return(-1);
  /* WE ALLOCATE MEMORY */
  b=(double*)malloc( sizeof(double)*(Na+2) );
  b2=(double*)malloc( sizeof(double)*(Na+2) );
  rx=(double*)malloc( sizeof(double)*(Na+2) );
  ry=(double*)malloc( sizeof(double)*(Na+2) );
  /* WE DEFINE THE POLYNOM WE NEED TO COMPUTE ROOTS */
  d=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
  b[0]=-1;
  b[1]=a[0];
  paso=d;
  for(i=2; i<(Na+2); i++)
  {
    b[i]=a[i-1]*paso;
    paso*=d;
  }

  if(Na == 2)
  {
    Nr=ami_RootCubicPolynomial(b,3,rx);
    for(i=0; i<Na; i++)
      ry[i]=0;
  }
  else
  {
    for(i=0; i<(Na+2); i++)
      b2[i]=b[Na+1-i];
    //for(i=0;i<(Na+2);i++) printf("b2[%d]=%1.20e\n",i,b2[i]); system("pause");
    Nr=ami_polynomial_root(b2,Na+1,rx,ry);
  }
  /* WE SELECT THE REAL ROOT NEAR TO 1 */
  root=10e5;
  for(i=0; i<Nr; i++)
  {
    //printf("rx[%d]=%e ry[%d]=%e\n",i,rx[i],i,ry[i]);
    if(fabs(ry[i])<0.00000000001 && fabs(root-1)>fabs(rx[i]-1))
      root=rx[i];
  }
  if(Nr == 0)
  {
    root=1.;
  }

  /* WE TRANSFORM THE POINT COORDINATES */
  *xt=x0+(x-x0)*root;
  *yt=y0+(y-y0)*root;

  double x2, y2;
  double error=ami_inverse_lens_distortion_newton_raphson(x,y,x0,y0,&x2,&y2,a,Na);

  if(error < 2.)
  {
    if( ((x-x2)*(x-x2)+(y-y2)*(y-y2)) < ((x-(*xt))*(x-(*xt))+(y-(*yt))*(y-(*yt))) )
    {
      *xt=x2; *yt=y2;
    }
  }

  free(b); 
  free(rx); 
  free(ry); 
  free(b2);
  
  return(0);
}

/**
 * \fn int ami_inverse_lens_distortion(double x,double y,double x0,double y0,double *xt,double *yt,double *a,int Na, double dl1r)
 *  \brief  Function to inverse the lens distortion transformation.
 *  \pre Any parameter can be null.
 *  \pre Na has to be positive.
 *  \param [in] x,y Point to inverse
 *  \param [in] x0,y0 Center of the image
 *  \param [out] xt,yt Inverve point transformed
 *  \param [in] a Lens distortion model polynom
 *  \param [in] Na Degree of the lens distortion model polynom
 *  \param [in] dl1r coeficient interpolated from vector of max distorsion distances
 *  \return 0
 * \author Luis Alvarez, Pedro Henriquez
 */
AMI_DLL_CPP int ami_inverse_lens_distortion_fast(
  double x,double y, /* POINT TO INVERSE (INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *xt,double *yt, /* INVERVE POINT TRANSFORMED (OUTPUT) */
  double *a, /* LENS DISTORTION MODEL POLYNOM */
  int Na,/* DEGREE OF THE LENS DISTORTION MODEL POLYNOM */
  double dl1r /* COEFICIENT INTERPOLATED FROM VECTOR OF MAX DISTORSION DISTANCES*/)
{
  /* WE TRANSFORM THE POINT COORDINATES */
  *xt=x0+(x-x0)*dl1r;
  *yt=y0+(y-y0)*dl1r;

  double x2,y2;
  double error=ami_inverse_lens_distortion_newton_raphson(x,y,x0,y0,&x2,&y2,a,Na);

  if(error < 2.)
  {
    if( ((x-x2)*(x-x2)+(y-y2)*(y-y2)) < ((x-(*xt))*(x-(*xt))+(y-(*yt))*(y-(*yt))) )
    {
      *xt=x2; *yt=y2;
    }
  }

  return(0);
}









AMI_DLL_CPP int ami_inverse_lens_distortion_old(
  double x,double y, /* POINT TO INVERSE (INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *xt,double *yt, /* INVERVE POINT TRANSFORMED (OUTPUT) */
  double *a, /* LENS DISTORTION MODEL POLYNOM */
  int Na) /* DEGREE OF THE LENS DISTORTION MODEL POLYNOM */
{
  int i,Nr;
  double paso, d, *b, *b2, *rx, *ry, root;
  if(a[Na]==0.)
    return(-1);
  /* WE ALLOCATE MEMORY */
  b=(double*)malloc( sizeof(double)*(Na+2) );
  b2=(double*)malloc( sizeof(double)*(Na+2) );
  rx=(double*)malloc( sizeof(double)*(Na+2) );
  ry=(double*)malloc( sizeof(double)*(Na+2) );
  /* WE DEFINE THE POLYNOM WE NEED TO COMPUTE ROOTS */
  d=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
  b[0]=-1; b[1]=a[0];//1.;
  paso=d;
  for(i=2; i<(Na+2); i++)
  {
    b[i]=a[i-1]*paso;
    paso*=d;
  }

  if(Na == 2)
  {
    Nr=ami_RootCubicPolynomial(b,3,rx);
    for(i=0; i<Na; i++)
      ry[i]=0;
  }
  else
  {
    for(i=0; i<(Na+2); i++)
      b2[i]=b[Na+1-i];
    //for(i=0;i<(Na+2);i++) printf("b2[%d]=%1.20e\n",i,b2[i]); system("pause");
    Nr=ami_polynomial_root(b2,Na+1,rx,ry);
  }
  /* WE SELECT THE REAL ROOT NEAR TO 1 */
  root=10e5;
  for(i=0; i<Nr; i++)
  {
    //printf("rx[%d]=%e ry[%d]=%e\n",i,rx[i],i,ry[i]);
    if(fabs(ry[i])<0.00000000001 && fabs(root-1)>fabs(rx[i]-1))
      root=rx[i];
  }
  if(Nr == 0)
    root=1.;

  /* WE TRANSFORM THE POINT COORDINATES */
  *xt=x0+(x-x0)*root;
  *yt=y0+(y-y0)*root;
  //printf("*xt=%lf,*yt=%lf\n",*xt,*yt);

  free(b); 
  free(rx); 
  free(ry); 
  free(b2);
  
  return(0);
}

/**
 * \fn double ami_lens_distortion_estimation(double **x,double **y,int Nl,int *Np,double x0,double y0,double *a,int Na,int k,double alfa)
 *  \brief  Function to compute the lens distortion model
 *  \pre Any parameter can be null.
 *  \pre Na has to be positive.
 *  \param [out] x,y Original coleccion of lines distribution
 *  \param [in] Nl Number of lines
 *  \param [in] Np Number of points for each line
 *  \param [in] x0,y0 Center of the image
 *  \param [out] a Lens Distortion Polynomial model
 *  \param [in] Na Degree of Polynomial model
 *  \param [in] k Coeficient of the lens distortion polynom model to be updated
 *  \param [in] alfa Weight for minimizing the square of the distance bewteen distorted and undistorted points
 *  \return Error
 * \author Luis Alvarez
 */
AMI_DLL_CPP double ami_lens_distortion_estimation(
  double **x,double **y, /* ORIGINAL COLECCION OF LINES DISTRIBUTION (INPUT)*/
  int Nl, /* NUMBER OF LINES */
  int *Np, /* NUMBER OF POINTS FOR EACH LINE(INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *a, /* Lens Distortion Polynomial model (INPUT-OUTPUT)*/
  int Na, /* Degree of Polynomial model (INPUT)*/
  int k, /* COEFICIENT OF THE LENS DISTORTION POLYNOM MODEL TO BE UPDATED*/
  double alfa) /* WEIGHT FOR MINIMIZING THE SQUARE OF THE DISTANCE BEWTEEN
              DISTORTED AND UNDISTORTED POINTS */
{
  double *pol, suma_dd, suma_Ad, d, A, Error=0;
  int m, i, j;

  /* WE ALLOCATE MEMORY */
  pol=(double*)malloc( sizeof(double)*5 );
  for(i=0; i<=4; i++)
    pol[i]=0.;

  /* WE ADAPT a[0] TO MINIMIZE THE SQUARE OF THE DISTANCE BEWTEEN
              DISTORTED AND UNDISTORDED POINTS */
  if(alfa > 0)
  {
    suma_dd=suma_Ad=0;
    for(m=0; m<Nl; m++)
    {
      for(i=0; i<Np[m]; i++)
      {
        d=sqrt( (x[m][i]-x0)*(x[m][i]-x0)+(y[m][i]-y0)*(y[m][i]-y0) );
        A=0;
        for(j=1; j<=Na; j++)
          A+=a[j]*pow(d,(double) j+1);
        suma_dd+=d*d;
        suma_Ad+=A*d;
      }
    }
    a[0]=1-suma_Ad/suma_dd;
  }


  /* WE COMPUTE THE LENS DISTORTION MODEL */
  for(i=0; i<Nl; i++)
  {
    ami_lens_distortion_polynomial_update(x[i],y[i],Np[i],a,Na,x0,y0,k,pol);
  }
  ami_lens_distortion_model_update(a,k,pol);
  //for(i=0;i<=Na;i++) printf("k=%d a[%d]=%e\n",k,i,a[i]);
  /* WE FREE THE MEMORY */
  free(pol);

  for(i=0; i<Nl; i++)
    Error += ami_LensDistortionEnergyError(x[i],y[i],Np[i],x0,y0,a,Na);
  return(Error/Nl);
}

/**
 * \fn void ami_lens_distortion_zoom_normalization(double **x,double **y,int Nl,int *Np,double x0,double y0,double *a,int Na)
 *  \brief  Not described
 *  \pre Any parameter can be null.
 *  \pre .
 * \author Luis Alvarez
 */
AMI_DLL_CPP void ami_lens_distortion_zoom_normalization(double **x/** Not described */,double **y/** Not described */,int Nl/** Not described */,int *Np/** Not described */,
   double x0/** Not described */,double y0/** Not described */,double *a/** Not described */,int Na/** Not described */)
{
  int i, k, m, N=0;
  double Z, d, suma_Ad, A;

  /* WE UPDATE a BY ESTIMATING A ZOOM FACTOR Z */
    suma_Ad=0;
    for(m=0; m<Nl; m++)
    {
      for(i=0; i<Np[m]; i++)
      {
        N++;
        d=sqrt( (x[m][i]-x0)*(x[m][i]-x0)+(y[m][i]-y0)*(y[m][i]-y0) );
        A=a[0];
        for(k=1; k<=Na; k++)
          A+=a[k]*pow(d,(double) k);
        suma_Ad += A*A;
      }
    }
    Z=sqrt((double) N/suma_Ad);
    for(k=0; k<=Na; k++)
      a[k]*=Z;
}
