/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


 #ifndef AMI_DLL_H
  #define AMI_DLL_H
#endif

/**
 * \file lens_distortion.h
 * \brief Functions for lens distortion model basic operations
 * \author Luis Alvarez \n \n
*/
#ifndef LENS_DISTORTION_H
#define LENS_DISTORTION_H

#ifdef __AMIDEBUG__
  #include "wxAmiDebugLog.h"
#endif

/*****************************************************************************************/
/*CALCULA LA RECTA QUE MEJOR APROXIMA UN CONJUNTO DE PUNTOS 2D  Y ELIMINA PUNTOS ALEJADOS
  TENIENDO EN CUENTA EL MODELO DE DISTORSION                                             */
/*****************************************************************************************/
AMI_DLL_H double ami_points_to_line_equation_outlier_elimination(
                                        double *a, int Na, double xc,double yc,
                                        double *x,double *y,int *N,
                                        double line[3],
                                        double outlier_elimination_factor);

/****************************************************************/
/*CALCULA LA RECTA QUE MEJOR APROXIMA UN CONJUNTO DE PUNTOS 2D
  INCLUDING THE POLINOMIAL DISTORTION MODEL */
/****************************************************************/
AMI_DLL_H double ami_points_to_line_equation(double *a, int Na, double xc,
                                             double yc, double *x,double *y,
                                             int N, double line[3]);


//=====================================================================
// DISTORTION MODEL ESTIMATION FROM A SET OF LINES
//=====================================================================
AMI_DLL_H double ami_distortion_model_estimation_2p(double  xc,double yc,
                                                    double **x,double **y,
                                                    int Nl,int *Np,double **a,
                                                    int *Na);

// ====================================================================
// COMPUTE THE LENS DISTORTION MODEL IN A POINT
// ====================================================================
AMI_DLL_H void ami_lens_distortion_model_evaluation(double *a,int Na, double xc,
                                                    double yc,double x_input,
                                                    double y_input,
                                                    double *x_output,
                                                    double *y_output);


/***************************************************************************
  FUNCTION TO ADD THE INFORMATION OF A LINE POINT SEQUENCE TO THE 4 DEGREE
  POLYNOMIAL TO COMPUTE THE LENS DISTORTION MODEL
***************************************************************************/
AMI_DLL_H int ami_lens_distortion_polynomial_update_distance_2v(
                                          double *x, double *y,int Np,double *a,
                                          int Na,double x0,double y0,int k1,
                                          int k2,double **pol,double alfa);

/*****************************************************************************
UPDATE OF THE LENS DISTORTION POLYNOMIAL MODEL FOR 2 VARIABLES. IF alfa>0
WE ADAPT a[0] TO MINIMIZE THE SQUARE DISTANCE BEEWTEN DISTORTED AND
UNDISTORTED POINTS AND WE ADD A TERM TO THE POLYNOMIAL ALSO MINIMIZING SUCH
DISTANCE WITH WEIGHT alfa
****************************************************************************/
AMI_DLL_H double ami_lens_distortion_estimation_2v(double **x,double **y,int Nl,
                                                   int *Np,double x0,double y0,
                                                   double *a,int Na,int k1,
                                                   int k2,double alfa);

/***************************************************************************
  FUNCTION TO UPDATE THE LENS DISTORTION MODEL BY MINIMIZING A 4 DEGREE
  2 VARIABLE POLYNOM
***************************************************************************/
AMI_DLL_H int ami_lens_distortion_model_update_2v(double *a,int k1,int k2,
                                                  double **pol);

/***************************************************************************
  FUNCTION TO ADD THE INFORMATION OF A LINE POINT SEQUENCE TO THE 4 DEGREE
  POLYNOMIAL TO COMPUTE THE LENS DISTORTION MODEL
***************************************************************************/
AMI_DLL_H int ami_lens_distortion_polynomial_update_2v(double *x, double *y,
                                                       int Np,double *a,int Na,
                                                       double x0,double y0,
                                                       int k1,int k2,
                                                       double **pol);

/****************************************************************************
 FUNCTION TO COMPUTE THE PARTIAL DERIVATIVES OF A 2 VARIABLE POLYNOM
THE DEGREE OF THE DERIVATIVE POLYNOMs IS ASSUMED TO BE THE SAME THAT
THE ORIGINAL ONE
****************************************************************************/
AMI_DLL_H void ami_2v_polynom_derivatives(double **p,int N,double **p_x,
                                          double **p_y);

/**********************************************************************
 FUNCTION TO EVALUATE THE DETERMINANT OF A MATRIX
 *********************************************************************/
AMI_DLL_H double ami_determinante(double **A,int N);
AMI_DLL_H double ami_determinante(double A[3][3]);

/************************************************************
 FUNCTION TO COMPUTE THE DETERMINANT OF A POLYNOM MATRIX
 ************************************************************/
AMI_DLL_H void ami_polynom_determinant(double p[6][6][19],int Np,int Nd,
                                       double *q);

/************************************************************************
 FUNCTION TO EVALUATE A 2 VARIABLE POLYNOM IN ONE POINT
 ***********************************************************************/
AMI_DLL_H double ami_2v_polynom_evaluation(double **pt1,int N1,double x,
                                           double y);

/**********************************************************************
 FUNCTION TO EVALUATE A 2 VARIABLE POLYNOM IN ONE OF THE VARIABLE VALUE. THE OUTPUT
IS A 1 DEGREE POLYNOM
*****************************************************************************/
AMI_DLL_H void ami_2v_polynom_to_1v_polynom(double **pt1,int N1,double *p3,
                                            double z,int flat);

/*****************************************************************************
 FUNCTION TO MULTIPLY POLINOMS OF 1 VARIABLE. THE RESULT IS ADDED TO THE OUTPUT
POLYNOM COEFICIENTS
*****************************************************************************/
AMI_DLL_H double* ami_1v_polynom_multiplication(double *pt1,int N1,double *pt2,
                                                int N2,double *pt3);

/*****************************************************************************
 FUNCTION TO MULTIPLY POLINOMS OF 2 VARIABLES.
*****************************************************************************/
AMI_DLL_H void ami_2v_polynom_multiplication(double **pt1,int N1,double **p2,
                                             int N2,double **p3);

/***********************************************************************
  FUNCTION TO COMPUTE THE REAL ROOTS OF A CUBIC POLYNOMIAL. IT RETURNS
  THE NUMBER OF ROOTS FOUND SORTED BY MAGNITUD
************************************************************************/
AMI_DLL_H int ami_RootCubicPolynomial(double *a,int N,double *x);

/**************************************************************
 EVALUATION OF A POLYNOM USING HORNER ALGORITHM
 ************************************************************/
AMI_DLL_H double ami_polynomial_evaluation(double *a,int Na,double x);

/***************************************************************************
  FUNCTION TO ADD THE INFORMATION OF A LINE POINT SEQUENCE TO THE 4 DEGREE
  POLYNOMIAL TO COMPUTE THE LENS DISTORTION MODEL
***************************************************************************/
AMI_DLL_H int ami_lens_distortion_polynomial_update(double *x, double *y,int Np,
                                                    double *a,int Na,double x0,
                                                    double y0,int k,double *pol);

/***************************************************************************
  FUNCTION TO UPDATE THE LENS DISTORTION MODEL BY MINIMIZING A 4 DEGREE
  POLYNOM
***************************************************************************/
AMI_DLL_H int ami_lens_distortion_model_update(double *a,int k,double *pol);

/****************************************************************************
    FUNCTION TO COMPUTE THE LENS DISTORTION ENERGY ERROR (THE RESIDUAL VARIANCE
    OF THE POINT DISTRIBUTION
*****************************************************************************/
AMI_DLL_H double ami_LensDistortionVarianceError(double *x,double *y,int Np,
                                                 double x0,double y0,double *a,
                                                 int Na);

/****************************************************************************
    FUNCTION TO COMPUTE THE LENS DISTORTION ENERGY ERROR
*****************************************************************************/
AMI_DLL_H double ami_LensDistortionEnergyError(double *x,double *y,int Np,
                                               double x0,double y0,double *a,
                                               int Na);

/*************************************************************************
       FUNCTION TO INVERSE THE LENS DISTORTION TRANSFORMATION
***********************************************************************/
AMI_DLL_H int ami_inverse_lens_distortion(double x,double y,double x0,double y0,
                                          double *xt,double *yt,double *a,
                                          int Na);

/*************************************************************************
  FUNCTION TO COMPUTE THE LENS DISTORTION MODEL
*************************************************************************/
AMI_DLL_H double ami_lens_distortion_estimation(double **x,double **y,int Nl,
                                                int *Np,double x0,double y0,
                                                double *a,int Na,int k,
                                                double alfa);

AMI_DLL_H void ami_lens_distortion_zoom_normalization(double **x,double **y,
                                                      int Nl,int *Np,double x0,
                                                      double y0,double *a,int Na);

int ami_inverse_lens_distortion_fast(double x,double y,double x0,double y0,
                                     double *xt,double *yt, double *a, int Na,
                                     double dl1r);

double ami_inverse_lens_distortion_newton_raphson(
  double x,double y, /* POINT TO INVERSE (INPUT)*/
  double x0,double y0, /* CENTER OF THE IMAGE (INPUT)*/
  double *xt,double *yt, /* INVERVE POINT TRANSFORMED (OUTPUT) */
  double *a, /* LENS DISTORTION MODEL POLYNOM */
  int Na); /* DEGREE OF THE LENS DISTORTION MODEL POLYNOM */
#endif
