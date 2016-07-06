/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#include "lens_distortion_procedures.h"
#include "lens_distortion.h"
#include "../ami_utilities/utilities.h"
#include <iostream>
using namespace std;

//------------------------------------------------------------------------------

/**
 * \fn double distortion_points_to_line_equation(lens_distortion_model &d,
																								 line_points &l)
 * \brief LINE ESTIMATION FROM POINT AFTER APPLYING A LENS DISTORTION MODEL TO 
					THE POINTS. RETURN -1 IF IT DOES NOT WORK. OTHERWISE RETURN THE 
					AVERAGE OF THE SQUARED DISTANCE OF THE POINTS TO THE LINE
					AFTER APPLYING A LENS DISTORTION MODEL
 * \param d INPUT DISTORTION CENTER
 * \param l OUTPUT LINE
 * \author Luis Alvarez
 */
double distortion_points_to_line_equation(
    lens_distortion_model &d,line_points &l)
{
  int i, j, k;
  long double suu, suv, svv/*,h0,h1*/, um, vm, h, r[4][3], min, paso, norma;
  long double cero=10e-100;
  int N=l.get_points().size();
  double a,b,c;

  if(N < 2)
  {
    printf("Numero de puntos para el Calculo de la recta 2D menor que 2\n");
    return(-1.);
  }

  // WE CREATE AN AUXILIAR VECTOR OF POINTS
  vector< point2d<double> > points(l.get_points().size());
  for(i=0; i<((int)points.size()); i++)
    points[i]=d.evaluation(l.get_points()[i]);

  suu=0;
  suv=0;
  svv=0;
  um=0;
  vm=0;
  
  for(i=0; i<N; i++)
  {
    um+=points[i].x;
    vm+=points[i].y;
  }
  um/=N; vm/=N;
  for(i=0; i<N; i++)
  {
    suu+=(points[i].x-um)*(points[i].x-um);
    svv+=(points[i].y-vm)*(points[i].y-vm);
    suv+=(points[i].x-um)*(points[i].y-vm);
  }
  suu/=N; 
  svv/=N;
  suv/=N;
  
  if(fabs(suv)<= cero)
  {
    if(suu<svv && svv>cero)
    {
      //a=1; b=0; c=-um;
      a=1.;
      b=0.;
      c=-um;
      l.set_abc(a,b,c);
      return(0.);
    }
    if(svv<suu && suu>cero)
    {
      //a=0; b=1; c=-vm;
      a=0.;
      b=1.;
      c=-vm;
      l.set_abc(a,b,c);
      return(0.);
    }
    printf("No se pudo calcular la recta 2D\n");
    return(-1);
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

  min=0.;
  k=0;
  for(i=0; i<N; i++)
  {
    paso=r[0][0]*points[i].x+r[0][1]*points[i].y+r[0][2];
    min += paso*paso;
  }
  
  for(j=1; j<4; j++)
  {
    h=0;
    for(i=0; i<N; i++)
    {
      paso=r[j][0]*points[i].x+r[j][1]*points[i].y+r[j][2];
      h += paso*paso;
    }
    if(h < min)
    {
      k=j;
      min=h;
    }
  }

  l.set_abc((double) r[k][0],(double)r[k][1],(double)r[k][2]);

  return min;
}

//------------------------------------------------------------------------------

/**
 * \fn distortion_points_to_line_equation_quotient(lens_distortion_model &d,
                                                    line_points &l)
 * \brief LINE ESTIMATION FROM POINT AFTER APPLYING A LENS DISTORTION MODEL TO
          THE POINTS. RETURN -1 IF IT DOES NOT WORK. OTHERWISE RETURN THE
          AVERAGE OF THE SQUARED DISTANCE OF THE POINTS TO THE LINE
          AFTER APPLYING A LENS DISTORTION MODEL
 * \param d INPUT DISTORTION CENTER
 * \param l OUTPUT LINE
 * \author Luis Alvarez
 */
double distortion_points_to_line_equation_quotient(
    lens_distortion_model &d,line_points &l)
{
  int i, j, k;
  long double suu, suv, svv/*,h0,h1*/, um, vm, h, r[4][3], min, paso, norma;
  long double cero=10e-100;
  int N=l.get_points().size();
  double a,b,c;

  if(N < 2)
  {
    printf("Numero de puntos para el Calculo de la recta 2D menor que 2\n");
    return(-1.);
  }

  // WE CREATE AN AUXILIAR VECTOR OF POINTS
  vector< point2d<double> > points(l.get_points().size());
  for(i=0;i<((int)points.size());i++)
    points[i]=d.evaluation_quotient(l.get_points()[i]);

  suu=0; 
  suv=0; 
  svv=0; 
  um=0; 
  vm=0;
  for(i=0; i<N; i++)
  {
    um+=points[i].x;
    vm+=points[i].y;
  }
  um/=N;
  vm/=N;
  for(i=0; i<N; i++)
  {
    suu+=(points[i].x-um)*(points[i].x-um);
    svv+=(points[i].y-vm)*(points[i].y-vm);
    suv+=(points[i].x-um)*(points[i].y-vm);
  }
  suu/=N; 
  svv/=N; 
  suv/=N;

  if(fabs(suv) <= cero)
  {
    if(suu<svv && svv>cero)
    {
      //a=1; b=0; c=-um;
      a=1.;
      b=0.;
      c=-um;
      l.set_abc(a,b,c);
      return(0.);
    }
    if(svv<suu && suu>cero)
    {
      //a=0; b=1; c=-vm;
      a=1.;
      b=0.;
      c=-vm;
      l.set_abc(a,b,c);
      return(0.);
    }
    printf("No se pudo calcular la recta 2D\n");
    return(-1);
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
    paso=r[0][0]*points[i].x+r[0][1]*points[i].y+r[0][2];
    min+=paso*paso;
  }
  for(j=1; j<4; j++)
  {
    h=0;
    for(i=0; i<N; i++)
    {
      paso=r[j][0]*points[i].x+r[j][1]*points[i].y+r[j][2];
      h+=paso*paso;
    }
    if(h < min)
    {
      k=j;
      min=h;
    }
  }

  l.set_abc((double) r[k][0],(double)r[k][1],(double)r[k][2]);

  return min;
}

//------------------------------------------------------------------------------

//DISTORTION PARAMETER AND CENTER OPTIMIZATION WITH POLYNOMIAL MODEL (2 PARAMETERS)
double model_center_estimation_2p_polynomial(
   point2d<double>  distortion_center /** INPUT DISTORTION CENTER */,
   std::vector<line_points> &lines /** INPUT/OUTPUT VECTOR OF LINES TO COMPUTE */,
   lens_distortion_model &d /** OUTPUT DISTORTION MODEL */,
   int w, int h,
   const std::vector<bool>& v)
{
  lens_distortion_model ldm;
  bool debug = false;
  bool filedebug = false;
  double h1 = 1e-4, h2 = 1e-2;
  double E_p1_p2_xc_yc, E_p1plush_p2_xc_yc, E_p1_p2plush_xc_yc, E_p1_p2_xcplush_yc,
         E_p1_p2_xc_ycplush, E_p1minush_p2_xc_yc, E_p1_p2minush_xc_yc,
         E_p1_p2_xcminush_yc, E_p1_p2_xc_ycminush, E_p1plush_p2plush_xc_yc,
         E_p1plush_p2_xcplush_yc, E_p1plush_p2_xc_ycplush, E_p1_p2plush_xcplush_yc,
         E_p1_p2plush_xc_ycplush, E_p1_p2_xcplush_yc_plush, E1_p1_p2_xc_yc;
	E_p1_p2_xc_yc=E_p1plush_p2_xc_yc=E_p1_p2plush_xc_yc=E_p1_p2_xcplush_yc=
  E_p1_p2_xc_ycplush=E_p1minush_p2_xc_yc=E_p1_p2minush_xc_yc=
  E_p1_p2_xcminush_yc=E_p1_p2_xc_ycminush=E_p1plush_p2plush_xc_yc=
  E_p1plush_p2_xcplush_yc=E_p1plush_p2_xc_ycplush=E_p1_p2plush_xcplush_yc=
  E_p1_p2plush_xc_ycplush=E_p1_p2_xcplush_yc_plush=E1_p1_p2_xc_yc=0.;

  double gamma = 10.0;

  double prevk1 = d.get_d()[1], prevk2 = d.get_d()[2];
  double nextk1 = 1.0,          nextk2 = 1.0;
  double dmi = update_rsqmax(distortion_center,w,h);
  //We compute p values using the relation r1=2r2
  double r2 = (sqrt(dmi))/2;
  double r2_2 = r2*r2;
  double r2_4 = r2*r2*r2*r2;
  //p1 = k14r_2^2 + k216r_2^4
  double prevp1 = prevk1*4*r2_2 + prevk2*16*r2_4;
  //p2 = k1r_2^2 + k2r_2^4
  double prevp2 = prevk1*r2_2 + prevk2*r2_4;
  double nextp1 = 0.0, nextp2 = 0.0;
  double d0[4] = {prevp1, prevp2, distortion_center.x, distortion_center.y};
  if(debug)
  {
    cout << "Input: p1 = " << d0[0] <<
                 " p2 = " << d0[1] <<
                 " xc = " << d0[2] <<
                 " yc = " << d0[3] << endl;
  }
  double d1[4];
  point2d<double> dc0 = distortion_center;
  point2d<double> dc1 = distortion_center;
  dc1.x += h2;
  dc1.y += h2;

  int convergence_it = 0;
  ldm.get_d().resize(3);
  ldm.get_d()[0] = 1.;
  double *gradE;
  ami2_malloc1d(gradE,double,4);
  double **hessianE;
  ami2_malloc2d(hessianE,double,4,4);
  //Debug file
  ofstream res;//("name.txt", ios::app);

	double TOL = 1e-4;
	
  while((((fabs(prevk1-nextk1))>((fabs(nextk1)+1e-30)*TOL))
        ||((fabs(prevk2-nextk2))>((fabs(nextk2)+1e-30)*TOL))
        ||((fabs(dc0.x-dc1.x))>((fabs(dc1.x)+1e-30)*TOL))
        ||((fabs(dc0.y-dc1.y))>((fabs(dc1.y)+1e-30)*TOL)))
        && convergence_it<=100)
  {
    //----1ST WE COMPUTE THE Es
    //INITIALIZATION OF THE ACCUMULATORS
    E_p1_p2_xc_yc = E_p1plush_p2_xc_yc = E_p1_p2plush_xc_yc =
    E_p1_p2_xcplush_yc = E_p1_p2_xc_ycplush = E_p1minush_p2_xc_yc =
    E_p1_p2minush_xc_yc = E_p1_p2_xcminush_yc = E_p1_p2_xc_ycminush =
    E_p1plush_p2plush_xc_yc = E_p1plush_p2_xcplush_yc =
    E_p1plush_p2_xc_ycplush = E_p1_p2plush_xcplush_yc =
    E_p1_p2plush_xc_ycplush = E_p1_p2_xcplush_yc_plush = 0.0;
    dc0 = point2d<double>(d0[2],d0[3]);
    point2d<double> dc;
    dmi = update_rsqmax(dc0,w,h);
    r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
    prevk1 = (prevp1-16*prevp2)/(4*r2_2-16*r2_2);
    prevk2 = (4*prevp2-prevp1)/(4*r2_4-16*r2_4);
    for(int i=0; i<(int)lines.size(); i++)
    {
      dc = dc0;
      //E(p1,p2,xc,yc)
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = (prevp1-16*prevp2)/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*prevp2-prevp1)/(4*r2_4-16*r2_4);
      E_p1_p2_xc_yc += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1+h1,p2,xc,yc)
      ldm.get_d()[1] = ((prevp1+h1)-16*prevp2)/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*prevp2-(prevp1+h1))/(4*r2_4-16*r2_4);
      E_p1plush_p2_xc_yc += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1,p2+h1,xc,yc)
      ldm.get_d()[1] = (prevp1-16*(prevp2+h1))/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*(prevp2+h1)-prevp1)/(4*r2_4-16*r2_4);
      E_p1_p2plush_xc_yc += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1,p2,xc+h2,yc)
      dc.x = dc0.x+h2;
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = (prevp1-16*prevp2)/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*prevp2-prevp1)/(4*r2_4-16*r2_4);
      E_p1_p2_xcplush_yc += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1,p2,xc,yc+h2)
      dc.x = dc0.x;
      dc.y = dc0.y+h2;
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = (prevp1-16*prevp2)/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*prevp2-prevp1)/(4*r2_4-16*r2_4);
      E_p1_p2_xc_ycplush += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1-h1,p2,xc,yc)
      dc.y = dc0.y;
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = ((prevp1-h1)-16*prevp2)/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*prevp2-(prevp1-h1))/(4*r2_4-16*r2_4);
      E_p1minush_p2_xc_yc += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1,p2-h1,xc,yc)
      ldm.get_d()[1] = (prevp1-16*(prevp2-h1))/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*(prevp2-h1)-prevp1)/(4*r2_4-16*r2_4);
      E_p1_p2minush_xc_yc += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1,p2,xc-h2,yc)
      dc.x = dc0.x-h2;
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = (prevp1-16*prevp2)/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*prevp2-prevp1)/(4*r2_4-16*r2_4);
      E_p1_p2_xcminush_yc += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1,p2,xc,yc-h2)
      dc.x = dc0.x;
      dc.y = dc0.y-h2;
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = (prevp1-16*prevp2)/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*prevp2-prevp1)/(4*r2_4-16*r2_4);
      E_p1_p2_xc_ycminush += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1+h1,p2+h1,xc,yc)
      dc.y = dc0.y;
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = ((prevp1+h1)-16*(prevp2+h1))/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*(prevp2+h1)-(prevp1+h1))/(4*r2_4-16*r2_4);
      E_p1plush_p2plush_xc_yc += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1+h1,p2,xc+h2,yc)
      dc.x = dc0.x+h2;
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = ((prevp1+h1)-16*prevp2)/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*prevp2-(prevp1+h1))/(4*r2_4-16*r2_4);
      E_p1plush_p2_xcplush_yc += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1+h1,p2,xc,yc+h2)
      dc.x = dc0.x;
      dc.y = dc0.y+h2;
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = ((prevp1+h1)-16*prevp2)/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*prevp2-(prevp1+h1))/(4*r2_4-16*r2_4);
      E_p1plush_p2_xc_ycplush += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1,p2+h1,xc+h2,yc)
      dc.x = dc0.x+h2;
      dc.y = dc0.y;
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = (prevp1-16*(prevp2+h1))/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*(prevp2+h1)-prevp1)/(4*r2_4-16*r2_4);
      E_p1_p2plush_xcplush_yc += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1,p2+h1,xc,yc+h2)
      dc.x = dc0.x;
      dc.y = dc0.y+h2;
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = (prevp1-16*(prevp2+h1))/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*(prevp2+h1)-prevp1)/(4*r2_4-16*r2_4);
      E_p1_p2plush_xc_ycplush += distortion_points_to_line_equation(ldm,lines[i]);

      //E(p1,p2,xc+h2,yc+h2)
      dc.x = dc0.x+h2;
      ldm.set_distortion_center(dc);
      dmi = update_rsqmax(dc,w,h);
      r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
      ldm.get_d()[1] = (prevp1-16*prevp2)/(4*r2_2-16*r2_2);
      ldm.get_d()[2] = (4*prevp2-prevp1)/(4*r2_4-16*r2_4);
      E_p1_p2_xcplush_yc_plush += distortion_points_to_line_equation(ldm,lines[i]);
    }
    //GRADIENT VECTOR
    /*-(E(p1+h1,p2,xc,yc) - E(p1,p2,xc,yc))/h1
      -(E(p1,p2+h1,xc,yc) - E(p1,p2,xc,yc))/h1
      -(E(p1,p2,xc+h2,yc) - E(p1,p2,xc,yc))/h2
      -(E(p1,p2,xc,yc+h2) - E(p1,p2,xc,yc))/h2*/
    if(v[0])
      gradE[0] = -((E_p1plush_p2_xc_yc - E_p1_p2_xc_yc)/h1);
    else
      gradE[0] = 0.;
    if(v[1])
      gradE[1] = -((E_p1_p2plush_xc_yc - E_p1_p2_xc_yc)/h1);
    else
      gradE[1] = 0.;
    if(v[2])
      gradE[2] = -((E_p1_p2_xcplush_yc - E_p1_p2_xc_yc)/h2);
    else
      gradE[2] = 0.;
    if(v[3])
      gradE[3] = -((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2);
    else
      gradE[3] = 0.;
    if(filedebug)
    {
      res << "----------------------" << endl;
      res << "v[0] = " << v[0] << " v[1] = " << v[1]
          << " v[2] = " << v[2] << " v[3] = " << v[3] << endl;
      res << "Iter. = " << (int)convergence_it << endl;
      res << "grad[0] = " << gradE[0] << endl;
      res << "grad[1] = " << gradE[1] << endl;
      res << "grad[2] = " << gradE[2] << endl;
      res << "grad[3] = " << gradE[3] << endl;
    }

    //HESSIAN MATRIX
    /*HE00: ((E(p1+h1,p2,xc,yc)+E(p1-h1,p2,xc,yc)-2E(p1,p2,xc,yc))
             /h1^2)+gamma*/
    hessianE[0][0] = ((E_p1plush_p2_xc_yc + E_p1minush_p2_xc_yc -
                       2*E_p1_p2_xc_yc)/(h1*h1))+gamma;

    /*HE01=HE10: (1/h1)*(((E(p1+h1,p2+h1,xc,yc)-E(p+h1,p2,xc,yc))/h1)-
                         ((E(p1,p2+h1,xc,yc)-E(p1,p2,xc,yc))/h1))*/
    hessianE[0][1] = hessianE[1][0] = (1/h1)*
                    (((E_p1plush_p2plush_xc_yc - E_p1plush_p2_xc_yc)/h1)-
                     ((E_p1_p2plush_xc_yc - E_p1_p2_xc_yc)/h1));

    /*HE02=HE20: (1/h1)*(((E(p1+h1,p2,xc+h2,yc)-E(p1+h1,p2,xc,yc))/h2)-
                         ((E(p1,p2,xc+h2,yc)-E(p1,p2,xc,yc))/h2))*/
    hessianE[0][2] = hessianE[2][0] = (1/h1)*
                    (((E_p1plush_p2_xcplush_yc - E_p1plush_p2_xc_yc)/h2)-
                     ((E_p1_p2_xcplush_yc - E_p1_p2_xc_yc)/h2));

    /*HE03=HE30: (1/h1)*(((E(p1+h1,p2,xc,yc+h2)-E(p1+h1,p2,xc,yc))/h2)-
                         ((E(p1,p2,xc,yc+h2)-E(p1,p2,xc,yc))/h2))*/
    hessianE[0][3] = hessianE[3][0] = (1/h1)*
                    (((E_p1plush_p2_xc_ycplush - E_p1plush_p2_xc_yc)/h2)-
                     ((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2));

    /*HE11:((E(p1,p2+h1,xc,yc)+E(p1,p2-h1,xc,yc)-2E(p1,p2,xc,yc))
             /h1^2)+gamma*/
    hessianE[1][1] = ((E_p1_p2plush_xc_yc + E_p1_p2minush_xc_yc -
                       2*E_p1_p2_xc_yc)/(h1*h1))+gamma;

    /*HE12=HE21: (1/h1)*(((E(p1,p2+h1,xc+h2,yc)-E(p1,p2+h1,xc,yc))/h2)-
                         ((E(p1,p2,xc+h2,yc)-E(p1,p2,xc,yc))/h2))*/
    hessianE[1][2] = hessianE[2][1] = (1/h1)*
                    (((E_p1_p2plush_xcplush_yc - E_p1_p2plush_xc_yc)/h2)-
                     ((E_p1_p2_xcplush_yc - E_p1_p2_xc_yc)/h2));

    /*HE13=HE31: (1/h1)*(((E(p1,p2+h1,xc,yc+h2)-E(p1,p2+h1,xc,yc))/h2)-
                         ((E(p1,p2,xc,yc+h2)-E(p1,p2,xc,yc))/h2))*/
    hessianE[1][3] = hessianE[3][1] = (1/h1)*
                    (((E_p1_p2plush_xc_ycplush - E_p1_p2plush_xc_yc)/h2)-
                     ((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2));

    /*HE22:((E(p1,p2,xc+h2,yc)+E(p1,p2,xc-h2,yc)-2E(p1,p2,xc,yc))
             /h2^2)+gamma*/
    hessianE[2][2] = ((E_p1_p2_xcplush_yc + E_p1_p2_xcminush_yc -
                       2*E_p1_p2_xc_yc)/(h2*h2))+gamma;

    /*HE23=HE32: (1/h2)*(((E(p1,p2,xc+h2,yc+h2)-E(p1,p2,xc+h2,yc))/h2)-
                         ((E(p1,p2,xc,yc+h2)-E(p1,p2,xc,yc))/h2))*/
    hessianE[2][3] = hessianE[3][2] = (1/h2)*
                    (((E_p1_p2_xcplush_yc_plush - E_p1_p2_xcplush_yc)/h2)-
                     ((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2));

    /*HE33:((E(p1,p2,xc,yc+h2)+E(p1,p2,xc,yc-h2)-2E(p1,p2,xc,yc))
             /h2^2)+gamma*/
    hessianE[3][3] = ((E_p1_p2_xc_ycplush + E_p1_p2_xc_ycminush -
                       2*E_p1_p2_xc_yc)/(h2*h2))+gamma;
    if(v[0]==false)
    {
      hessianE[0][0]=1.+gamma;
      hessianE[0][1]=hessianE[0][2]=hessianE[0][3]=
      hessianE[1][0]=hessianE[2][0]=hessianE[3][0]=0.;
    }
    if(v[1]==false)
    {
      hessianE[1][1]=1.+gamma;
      hessianE[1][0]=hessianE[1][2]=hessianE[1][3]=
      hessianE[0][1]=hessianE[2][1]=hessianE[3][1]=0.;
    }
    if(v[2]==false)
    {
      hessianE[2][2]=1.+gamma;
      hessianE[2][0]=hessianE[2][1]=hessianE[2][3]=
      hessianE[0][2]=hessianE[1][2]=hessianE[3][2]=0.;
    }
    if(v[3]==false)
    {
      hessianE[3][3]=1.+gamma;
      hessianE[3][0]=hessianE[3][1]=hessianE[3][2]=
      hessianE[0][3]=hessianE[1][3]=hessianE[2][3]=0.;
    }
    if(filedebug)
    {
      res << "hessian[0][0] = " << hessianE[0][0] << endl;
      res << "hessian[0][1] = hessian[1][0] = " << hessianE[0][1] << endl;
      res << "hessian[0][2] = hessian[2][0] = " << hessianE[0][2] << endl;
      res << "hessian[0][3] = hessian[3][0] = " << hessianE[0][3] << endl;
      res << "hessian[1][1] = " << hessianE[1][1] << endl;
      res << "hessian[1][2] = hessian[2][1] = " << hessianE[1][2] << endl;
      res << "hessian[1][3] = hessian[3][1] = " << hessianE[1][3] << endl;
      res << "hessian[2][2] = " << hessianE[2][2] << endl;
      res << "hessian[2][3] = hessian[3][2] = " << hessianE[2][3] << endl;
      res << "hessian[3][3] = " << hessianE[3][3] << endl;
      res << "With gamma = " << gamma << endl;
      res << "----------------------" << endl;
    }

    //WE SOLVE THE SYSTEM
    if(ami2_gauss(hessianE,gradE,4)!=0)
        cout << "The system hasn't solution" << endl;

    //WE COMPUTE THE ERROR WITH THE RESULT OF THE SYSTEM
    //d1 = d0 + z
    d1[0] = d0[0] + gradE[0];
    d1[1] = d0[1] + gradE[1];
    d1[2] = d0[2] + gradE[2];
    d1[3] = d0[3] + gradE[3];
    dc1.x = d1[2];
    dc1.y = d1[3];
    ldm.set_distortion_center(dc1);
    dmi = update_rsqmax(dc1,w,h);
    r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;

    nextp1 = d1[0];
    nextp2 = d1[1];
    nextk1 = ldm.get_d()[1] = (nextp1-16*nextp2)/(4*r2*r2-16*r2*r2);
    nextk2 = ldm.get_d()[2] = (4*nextp2-nextp1)/(4*r2*r2*r2*r2-16*r2*r2*r2*r2);
    E1_p1_p2_xc_yc = 0.0;
    for(int i=0; i<(int)lines.size(); i++)
      E1_p1_p2_xc_yc += distortion_points_to_line_equation(ldm,lines[i]);

    if(filedebug)
    {
      res << "E0(p1,p2,xc0,yc0) = " << E_p1_p2_xc_yc  << " p1 = "  << d0[0]
                                                      << " p2 = " << d0[1]
                                                      << " xc0 = " << d0[2]
                                                      << " yc0 = " << d0[3]
                                                      << endl;
      res << "E1(p1,p2,xc1,yc1) = " << E1_p1_p2_xc_yc << " p1 = "  << d1[0]
                                                      << " p2 = " << d1[1]
                                                      << " xc1 = " << d1[2]
                                                      << " yc1 = " << d1[3]
                                                      << endl;
    }

		//TEMPORAL MODEL FOR CHECKING THE INVERTIBILITY
		dmi = update_rsqmax(point2d<double>(d1[2],d1[3]),w,h);
		r2 = (sqrt(dmi))/2;//, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
		lens_distortion_model tmp_ldm;
    tmp_ldm.get_d().resize(3);
    tmp_ldm.get_d()[0] = 0.;
		tmp_ldm.get_d()[1] = (d1[0]-16*d1[1])/(4*r2*r2-16*r2*r2);
		tmp_ldm.get_d()[2] = (4*d1[1]-d1[0])/(4*r2*r2*r2*r2-16*r2*r2*r2*r2);
		tmp_ldm.set_distortion_center(point2d<double>(d1[2],d1[3]));
    tmp_ldm.set_type(POLYNOMIAL);
    
		if(E1_p1_p2_xc_yc < E_p1_p2_xc_yc && check_invertibility(tmp_ldm,w,h)==true)
    {
      prevp1 = nextp1;
      prevp2 = nextp2;
      d0[0]  = d1[0];
      d0[1]  = d1[1];
      d0[2]  = d1[2];
      d0[3]  = d1[3];
      gamma /= 10;
    }
    else
    {
      gamma *= 10;
    }
    if(filedebug)
    {
      res << "gamma update = " << gamma << endl;
      res << "---------------------------------------------------------------" << endl;
    }
    convergence_it++;
  }
  free(gradE);
  ami2_free2d(hessianE);
  if(filedebug)
  {
    res << "-------------------------------------------------------------------" << endl;
  }
  res.close();
  if(debug)
  {
    cout << "No. of iterations = " << convergence_it << endl;
    cout << "gamma = " << gamma << endl;
    cout << "Output: p1 = " << d0[0] <<
                  " p2 = " << d0[1] <<
                  " xc = " << d0[2] <<
                  " yc = " << d0[3] << endl;
  }

  //WE UPDATE THE FINAL MODEL
  dmi = update_rsqmax(point2d<double>(d0[2],d0[3]),w,h);
  r2 = (sqrt(dmi))/2, r2_2 = r2*r2, r2_4 = r2*r2*r2*r2;
  d.get_d()[1] = (d0[0]-16*d0[1])/(4*r2*r2-16*r2*r2);
  d.get_d()[2] = (4*d0[1]-d0[0])/(4*r2*r2*r2*r2-16*r2*r2*r2*r2);
  d.set_distortion_center(point2d<double>(d0[2],d0[3]));

  //AND RETURN THE ERROR VALUE
  if(E1_p1_p2_xc_yc<E_p1_p2_xc_yc)
      return E1_p1_p2_xc_yc;
  else
      return E_p1_p2_xc_yc;
}

//------------------------------------------------------------------------------

//DISTORTION PARAMETER AND CENTER OPTIMIZATION WITH DIVISION MODEL (2 PARAMETERS)
double model_center_estimation_2p_quotient(
   point2d<double>  distortion_center /** INPUT DISTORTION CENTER */,
   std::vector<line_points> &lines /** INPUT/OUTPUT VECTOR OF LINES TO COMPUTE */,
   lens_distortion_model &d /** OUTPUT DISTORTION MODEL */,
   int w /**IMAGE WIDTH*/,
   int h /**IMAGE HEIGHT*/,
   const std::vector<bool>& v)
{
	lens_distortion_model ldm;
	bool debug = false;
	bool filedebug = false;
	double h1 = 1e-4, h2 = 1e-2;
	double E_p1_p2_xc_yc, E_p1plush_p2_xc_yc, E_p1_p2plush_xc_yc, E_p1_p2_xcplush_yc,
				 E_p1_p2_xc_ycplush, E_p1minush_p2_xc_yc, E_p1_p2minush_xc_yc,
				 E_p1_p2_xcminush_yc, E_p1_p2_xc_ycminush, E_p1plush_p2plush_xc_yc,
				 E_p1plush_p2_xcplush_yc, E_p1plush_p2_xc_ycplush, E_p1_p2plush_xcplush_yc,
				 E_p1_p2plush_xc_ycplush, E_p1_p2_xcplush_yc_plush, E1_p1_p2_xc_yc;
	E_p1_p2_xc_yc=E_p1plush_p2_xc_yc=E_p1_p2plush_xc_yc=E_p1_p2_xcplush_yc=
	E_p1_p2_xc_ycplush=E_p1minush_p2_xc_yc=E_p1_p2minush_xc_yc=
	E_p1_p2_xcminush_yc=E_p1_p2_xc_ycminush=E_p1plush_p2plush_xc_yc=
	E_p1plush_p2_xcplush_yc=E_p1plush_p2_xc_ycplush=E_p1_p2plush_xcplush_yc=
	E_p1_p2plush_xc_ycplush=E_p1_p2_xcplush_yc_plush=E1_p1_p2_xc_yc=0.;

	double gamma = 10.0;

	double prevk1 = d.get_d()[1], prevk2 = d.get_d()[2];
	double nextk1 = 1.0,          nextk2 = 1.0;
	double dmi = update_rsqmax(distortion_center,w,h);
	//We compute r1 and r2
	double r1 = sqrt(dmi);
	double r2 = r1/2;
	//p1 = k14r_2^2 + k216r_2^4
	double prevp1 = (1/(1+prevk1*r1*r1+prevk2*r1*r1*r1*r1))-1;
	//p2 = k1r_2^2 + k2r_2^4
	double prevp2 = (1/(1+prevk1*r2*r2+prevk2*r2*r2*r2*r2))-1;
	double nextp1 = 0.0, nextp2 = 0.0;
	double d0[4] = {prevp1, prevp2, distortion_center.x, distortion_center.y};
	if(debug)
	{
		cout << "Input: p1 = " << d0[0] <<
								 " p2 = " << d0[1] <<
								 " xc = " << d0[2] <<
								 " yc = " << d0[3] << endl;
	}
	double d1[4];
	point2d<double> dc0 = distortion_center;
	point2d<double> dc1 = distortion_center;
	dc1.x += h2;
	dc1.y += h2;

	int convergence_it = 0;
	ldm.get_d().resize(3);
	ldm.get_d()[0] = 1.;
	double *gradE;
	ami2_malloc1d(gradE,double,4);
	double **hessianE;
	ami2_malloc2d(hessianE,double,4,4);
	//Debug file
	ofstream res;//("name.txt", ios::app);

	double TOL = 1e-4;
	
	while((((fabs(prevk1-nextk1))>((fabs(nextk1)+1e-30)*TOL))
				||((fabs(prevk2-nextk2))>((fabs(nextk2)+1e-30)*TOL))
				||((fabs(dc0.x-dc1.x))>((fabs(dc1.x)+1e-30)*TOL))
				||((fabs(dc0.y-dc1.y))>((fabs(dc1.y)+1e-30)*TOL)))
        && convergence_it<=100)
	{
		//----1ST WE COMPUTE THE Es
		//INITIALIZATION OF THE ACCUMULATORS
		E_p1_p2_xc_yc = E_p1plush_p2_xc_yc = E_p1_p2plush_xc_yc =
		E_p1_p2_xcplush_yc = E_p1_p2_xc_ycplush = E_p1minush_p2_xc_yc =
		E_p1_p2minush_xc_yc = E_p1_p2_xcminush_yc = E_p1_p2_xc_ycminush =
		E_p1plush_p2plush_xc_yc = E_p1plush_p2_xcplush_yc =
		E_p1plush_p2_xc_ycplush = E_p1_p2plush_xcplush_yc =
		E_p1_p2plush_xc_ycplush = E_p1_p2_xcplush_yc_plush = 0.0;
		dc0 = point2d<double>(d0[2],d0[3]);
		point2d<double> dc;
		dmi = update_rsqmax(dc0,w,h);
		r2 = sqrt(dmi)/2;
		prevk1 = (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2)))/(-12*r2*r2);
		prevk2 = (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1)))/(-12*r2*r2*r2*r2);
		for(int i=0; i<(int)lines.size(); i++)
		{
			dc = dc0;
			//E(p1,p2,xc,yc)
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1)))/(-12*r2*r2*r2*r2);
			E_p1_p2_xc_yc += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1+h1,p2,xc,yc)
			ldm.get_d()[1] = (((-(prevp1+h1))/(1+prevp1+h1))+((16*prevp2)/(1+prevp2)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*prevp2)/(1+prevp2))+((prevp1+h1)/(1+prevp1+h1)))/(-12*r2*r2*r2*r2);
			E_p1plush_p2_xc_yc += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1,p2+h1,xc,yc)
			ldm.get_d()[1] = (((-prevp1)/(1+prevp1))+((16*(prevp2+h1))/(1+prevp2+h1)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*(prevp2+h1))/(1+prevp2+h1))+(prevp1/(1+prevp1)))/(-12*r2*r2*r2*r2);
			E_p1_p2plush_xc_yc += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1,p2,xc+h2,yc)
			dc.x = dc0.x+h2;
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1)))/(-12*r2*r2*r2*r2);
			E_p1_p2_xcplush_yc += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1,p2,xc,yc+h2)
			dc.x = dc0.x;
			dc.y = dc0.y+h2;
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1)))/(-12*r2*r2*r2*r2);
			E_p1_p2_xc_ycplush += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1-h1,p2,xc,yc)
			dc.y = dc0.y;
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-(prevp1-h1))/(1+prevp1-h1))+((16*prevp2)/(1+prevp2)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*prevp2)/(1+prevp2))+((prevp1-h1)/(1+prevp1-h1)))/(-12*r2*r2*r2*r2);
			E_p1minush_p2_xc_yc += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1,p2-h1,xc,yc)
			ldm.get_d()[1] = (((-prevp1)/(1+prevp1))+((16*(prevp2-h1))/(1+prevp2-h1)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*(prevp2-h1))/(1+prevp2-h1))+(prevp1/(1+prevp1)))/(-12*r2*r2*r2*r2);
			E_p1_p2minush_xc_yc += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1,p2,xc-h2,yc)
			dc.x = dc0.x-h2;
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1)))/(-12*r2*r2*r2*r2);
			E_p1_p2_xcminush_yc += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1,p2,xc,yc-h2)
			dc.x = dc0.x;
			dc.y = dc0.y-h2;
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1)))/(-12*r2*r2*r2*r2);
			E_p1_p2_xc_ycminush += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1+h1,p2+h1,xc,yc)
			dc.y = dc0.y;
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-(prevp1+h1))/(1+prevp1+h1))+((16*(prevp2+h1))/(1+prevp2+h1)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*(prevp2+h1))/(1+prevp2+h1))+((prevp1+h1)/(1+prevp1+h1)))/(-12*r2*r2*r2*r2);
			E_p1plush_p2plush_xc_yc += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1+h1,p2,xc+h2,yc)
			dc.x = dc0.x+h2;
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-(prevp1+h1))/(1+prevp1+h1))+((16*prevp2)/(1+prevp2)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*prevp2)/(1+prevp2))+((prevp1+h1)/(1+prevp1+h1)))/(-12*r2*r2*r2*r2);
			E_p1plush_p2_xcplush_yc += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1+h1,p2,xc,yc+h2)
			dc.x = dc0.x;
			dc.y = dc0.y+h2;
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-(prevp1+h1))/(1+prevp1+h1))+((16*prevp2)/(1+prevp2)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*prevp2)/(1+prevp2))+((prevp1+h1)/(1+prevp1+h1)))/(-12*r2*r2*r2*r2);
			E_p1plush_p2_xc_ycplush += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1,p2+h1,xc+h2,yc)
			dc.x = dc0.x+h2;
			dc.y = dc0.y;
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-prevp1)/(1+prevp1))+((16*(prevp2+h1))/(1+prevp2+h1)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*(prevp2+h1))/(1+prevp2+h1))+(prevp1/(1+prevp1)))/(-12*r2*r2*r2*r2);
			E_p1_p2plush_xcplush_yc += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1,p2+h1,xc,yc+h2)
			dc.x = dc0.x;
			dc.y = dc0.y+h2;
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-prevp1)/(1+prevp1))+((16*(prevp2+h1))/(1+prevp2+h1)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*(prevp2+h1))/(1+prevp2+h1))+(prevp1/(1+prevp1)))/(-12*r2*r2*r2*r2);
			E_p1_p2plush_xc_ycplush += distortion_points_to_line_equation_quotient(ldm,lines[i]);

			//E(p1,p2,xc+h2,yc+h2)
			dc.x = dc0.x+h2;
			ldm.set_distortion_center(dc);
			dmi = update_rsqmax(dc,w,h);
			r2 = sqrt(dmi)/2;
			ldm.get_d()[1] = (((-prevp1)/(1+prevp1))+((16*prevp2)/(1+prevp2)))/(-12*r2*r2);
			ldm.get_d()[2] = (((-4*prevp2)/(1+prevp2))+(prevp1/(1+prevp1)))/(-12*r2*r2*r2*r2);
			E_p1_p2_xcplush_yc_plush += distortion_points_to_line_equation_quotient(ldm,lines[i]);
		}
		//GRADIENT VECTOR
		/*-(E(p1+h1,p2,xc,yc) - E(p1,p2,xc,yc))/h1
			-(E(p1,p2+h1,xc,yc) - E(p1,p2,xc,yc))/h1
			-(E(p1,p2,xc+h2,yc) - E(p1,p2,xc,yc))/h2
			-(E(p1,p2,xc,yc+h2) - E(p1,p2,xc,yc))/h2*/
		if(v[0])
			gradE[0] = -((E_p1plush_p2_xc_yc - E_p1_p2_xc_yc)/h1);
		else
			gradE[0] = 0.;
		if(v[1])
			gradE[1] = -((E_p1_p2plush_xc_yc - E_p1_p2_xc_yc)/h1);
		else
			gradE[1] = 0.;
		if(v[2])
			gradE[2] = -((E_p1_p2_xcplush_yc - E_p1_p2_xc_yc)/h2);
		else
			gradE[2] = 0.;
		if(v[3])
			gradE[3] = -((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2);
		else
			gradE[3] = 0.;
		if(filedebug)
		{
			res << "----------------------" << endl;
			res << "v[0] = " << v[0] << " v[1] = " << v[1]
				 << " v[2] = " << v[2] << " v[3] = " << v[3] << endl;
			res << "Iter. = " << convergence_it << endl;
			res << "grad[0] = " << gradE[0] << endl;
			res << "grad[1] = " << gradE[1] << endl;
			res << "grad[2] = " << gradE[2] << endl;
			res << "grad[3] = " << gradE[3] << endl;
		}

		//HESSIAN MATRIX
		/*HE00: ((E(p1+h1,p2,xc,yc)+E(p1-h1,p2,xc,yc)-2E(p1,p2,xc,yc))
					 /h1^2)+gamma*/
		hessianE[0][0] = ((E_p1plush_p2_xc_yc + E_p1minush_p2_xc_yc -
										 2*E_p1_p2_xc_yc)/(h1*h1))+gamma;

		/*HE01=HE10: (1/h1)*(((E(p1+h1,p2+h1,xc,yc)-E(p+h1,p2,xc,yc))/h1)-
												 ((E(p1,p2+h1,xc,yc)-E(p1,p2,xc,yc))/h1))*/
		hessianE[0][1] = hessianE[1][0] = (1/h1)*
										(((E_p1plush_p2plush_xc_yc - E_p1plush_p2_xc_yc)/h1)-
										 ((E_p1_p2plush_xc_yc - E_p1_p2_xc_yc)/h1));

		/*HE02=HE20: (1/h1)*(((E(p1+h1,p2,xc+h2,yc)-E(p1+h1,p2,xc,yc))/h2)-
												 ((E(p1,p2,xc+h2,yc)-E(p1,p2,xc,yc))/h2))*/
		hessianE[0][2] = hessianE[2][0] = (1/h1)*
										(((E_p1plush_p2_xcplush_yc - E_p1plush_p2_xc_yc)/h2)-
										 ((E_p1_p2_xcplush_yc - E_p1_p2_xc_yc)/h2));

		/*HE03=HE30: (1/h1)*(((E(p1+h1,p2,xc,yc+h2)-E(p1+h1,p2,xc,yc))/h2)-
												 ((E(p1,p2,xc,yc+h2)-E(p1,p2,xc,yc))/h2))*/
		hessianE[0][3] = hessianE[3][0] = (1/h1)*
										(((E_p1plush_p2_xc_ycplush - E_p1plush_p2_xc_yc)/h2)-
										 ((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2));
		/*HE11:((E(p1,p2+h1,xc,yc)+E(p1,p2-h1,xc,yc)-2E(p1,p2,xc,yc))
					 /h1^2)+gamma*/
		hessianE[1][1] = ((E_p1_p2plush_xc_yc + E_p1_p2minush_xc_yc -
											 2*E_p1_p2_xc_yc)/(h1*h1))+gamma;

		/*HE12=HE21: (1/h1)*(((E(p1,p2+h1,xc+h2,yc)-E(p1,p2+h1,xc,yc))/h2)-
												 ((E(p1,p2,xc+h2,yc)-E(p1,p2,xc,yc))/h2))*/
		hessianE[1][2] = hessianE[2][1] = (1/h1)*
										(((E_p1_p2plush_xcplush_yc - E_p1_p2plush_xc_yc)/h2)-
										 ((E_p1_p2_xcplush_yc - E_p1_p2_xc_yc)/h2));

		/*HE13=HE31: (1/h1)*(((E(p1,p2+h1,xc,yc+h2)-E(p1,p2+h1,xc,yc))/h2)-
												 ((E(p1,p2,xc,yc+h2)-E(p1,p2,xc,yc))/h2))*/
		hessianE[1][3] = hessianE[3][1] = (1/h1)*
										(((E_p1_p2plush_xc_ycplush - E_p1_p2plush_xc_yc)/h2)-
										 ((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2));
		/*HE22:((E(p1,p2,xc+h2,yc)+E(p1,p2,xc-h2,yc)-2E(p1,p2,xc,yc))
					 /h2^2)+gamma*/
		hessianE[2][2] = ((E_p1_p2_xcplush_yc + E_p1_p2_xcminush_yc -
											 2*E_p1_p2_xc_yc)/(h2*h2))+gamma;

		/*HE23=HE32: (1/h2)*(((E(p1,p2,xc+h2,yc+h2)-E(p1,p2,xc+h2,yc))/h2)-
												 ((E(p1,p2,xc,yc+h2)-E(p1,p2,xc,yc))/h2))*/
		hessianE[2][3] = hessianE[3][2] = (1/h2)*
										(((E_p1_p2_xcplush_yc_plush - E_p1_p2_xcplush_yc)/h2)-
										 ((E_p1_p2_xc_ycplush - E_p1_p2_xc_yc)/h2));
		/*HE33:((E(p1,p2,xc,yc+h2)+E(p1,p2,xc,yc-h2)-2E(p1,p2,xc,yc))
						 /h2^2)+gamma*/
		hessianE[3][3] = ((E_p1_p2_xc_ycplush + E_p1_p2_xc_ycminush -
											 2*E_p1_p2_xc_yc)/(h2*h2))+gamma;

		if(v[0]==false)
		{
			hessianE[0][0]=1.+gamma;
			hessianE[0][1]=hessianE[0][2]=hessianE[0][3]=
			hessianE[1][0]=hessianE[2][0]=hessianE[3][0]=0.;
		}
		if(v[1]==false)
		{
			hessianE[1][1]=1.+gamma;
			hessianE[1][0]=hessianE[1][2]=hessianE[1][3]=
			hessianE[0][1]=hessianE[2][1]=hessianE[3][1]=0.;
		}
		if(v[2]==false)
		{
			hessianE[2][2]=1.+gamma;
			hessianE[2][0]=hessianE[2][1]=hessianE[2][3]=
			hessianE[0][2]=hessianE[1][2]=hessianE[3][2]=0.;
		}
		if(v[3]==false)
		{
			hessianE[3][3]=1.+gamma;
			hessianE[3][0]=hessianE[3][1]=hessianE[3][2]=
			hessianE[0][3]=hessianE[1][3]=hessianE[2][3]=0.;
		}

		if(filedebug)
		{
			res << "hessian[0][0] = " << hessianE[0][0] << endl;
			res << "hessian[0][1] = hessian[1][0] = " << hessianE[0][1] << endl;
			res << "hessian[0][2] = hessian[2][0] = " << hessianE[0][2] << endl;
			res << "hessian[0][3] = hessian[3][0] = " << hessianE[0][3] << endl;
			res << "hessian[1][1] = " << hessianE[1][1] << endl;
			res << "hessian[1][2] = hessian[2][1] = " << hessianE[1][2] << endl;
			res << "hessian[1][3] = hessian[3][1] = " << hessianE[1][3] << endl;
			res << "hessian[2][2] = " << hessianE[2][2] << endl;
			res << "hessian[2][3] = hessian[3][2] = " << hessianE[2][3] << endl;
			res << "hessian[3][3] = " << hessianE[3][3] << endl;
			res << "----------------------" << endl;
		}

		//WE SOLVE THE SYSTEM
		if(ami2_gauss(hessianE,gradE,4)!=0)
				cout << "The system hasn't solution" << endl;

		//WE COMPUTE THE ERROR WITH THE RESULT OF THE SYSTEM
		//d1 = d0 + z
		d1[0] = d0[0] + gradE[0];
		d1[1] = d0[1] + gradE[1];
		d1[2] = d0[2] + gradE[2];
		d1[3] = d0[3] + gradE[3];
		dc1.x = d1[2];
		dc1.y = d1[3];
		ldm.set_distortion_center(dc1);
		dmi = update_rsqmax(dc1,w,h);
		r2 = sqrt(dmi)/2;

		nextp1 = d1[0];
		nextp2 = d1[1];
		nextk1 = ldm.get_d()[1] = (((-nextp1)/(1+nextp1))+((16*nextp2)/(1+nextp2)))/(-12*r2*r2);
		nextk2 = ldm.get_d()[2] = (((-4*nextp2)/(1+nextp2))+(nextp1/(1+nextp1)))/(-12*r2*r2*r2*r2);
		E1_p1_p2_xc_yc = 0.0;
		for(int i=0; i<(int)lines.size(); i++)
			E1_p1_p2_xc_yc += distortion_points_to_line_equation_quotient(ldm,lines[i]);

		if(filedebug)
		{
			res << "E0(p1,p2,xc0,yc0) = " << E_p1_p2_xc_yc  << " p1 = "  << d0[0]
																											<< " p2 = "  << d0[1]
																											<< " xc0 = " << d0[2]
																											<< " yc0 = " << d0[3]
																											<< endl;
			res << "E1(p1,p2,xc1,yc1) = " << E1_p1_p2_xc_yc << " p1 = "  << d1[0]
																											<< " p2 = "  << d1[1]
																											<< " xc1 = " << d1[2]
																											<< " yc1 = " << d1[3]
																											<< endl;
			res << "---------------------------------------------------------------" << endl;
		}

		//TEMPORAL MODEL FOR CHECKING THE INVERTIBILITY
		dmi = update_rsqmax(point2d<double>(d1[2],d1[3]),w,h);
		r2 = sqrt(dmi)/2;
		lens_distortion_model tmp_ldm;
    tmp_ldm.get_d().resize(3);
    tmp_ldm.get_d()[0] = 0.;
		tmp_ldm.get_d()[1] = (((-d1[0])/(1+d1[0]))+((16*d1[1])/(1+d1[1])))/(-12*r2*r2);
		tmp_ldm.get_d()[2] = (((-4*d1[1])/(1+d1[1]))+(d1[0]/(1+d1[0])))/(-12*r2*r2*r2*r2);
		tmp_ldm.set_distortion_center(point2d<double>(d1[2],d1[3]));
    tmp_ldm.set_type(DIVISION);
	
		if(E1_p1_p2_xc_yc < E_p1_p2_xc_yc && check_invertibility(tmp_ldm,w,h)==true)
		{
			prevp1 = nextp1;
			prevp2 = nextp2;
			d0[0]  = d1[0];
			d0[1]  = d1[1];
			d0[2]  = d1[2];
			d0[3]  = d1[3];
			gamma /= 10;
		}
		else
		{
			gamma *= 10;
		}
		convergence_it++;
	}
	free(gradE);
	ami2_free2d(hessianE);
	if(filedebug)
	{
		res << "-------------------------------------------------------------------" << endl;
	}
	res.close();
	if(debug)
	{
		cout << "No. of iterations = " << convergence_it << endl;
		cout << "gamma = " << gamma << endl;
		cout << "Output: p1 = " << d0[0] <<
									" p2 = " << d0[1] <<
									" xc = " << d0[2] <<
									" yc = " << d0[3] << endl;
	}

	//WE UPDATE THE FINAL MODEL
	dmi = update_rsqmax(point2d<double>(d0[2],d0[3]),w,h);
	r2 = sqrt(dmi)/2;
	d.get_d()[1] = (((-d0[0])/(1+d0[0]))+((16*d0[1])/(1+d0[1])))/(-12*r2*r2);
	d.get_d()[2] = (((-4*d0[1])/(1+d0[1]))+(d0[0]/(1+d0[0])))/(-12*r2*r2*r2*r2);
	d.set_distortion_center(point2d<double>(d0[2],d0[3]));

	//AND RETURN THE ERROR VALUE
	if(E1_p1_p2_xc_yc<E_p1_p2_xc_yc)
		return E1_p1_p2_xc_yc;
	else
		return E_p1_p2_xc_yc;
}

//------------------------------------------------------------------------------

/**
 * \fn int build_l1r_vector(vector<double> &l1r,point2d<double> &dc, 
														double max_distance_corner,int Na, double *a)
 * \brief Build an intermediate vector with values of L-1(r) for d = (0 to max_distance_corner)
 * \param [in] [out] l1r vector to store the roots
 * \param [in] dc distortion center
 * \param [in] max_distance_corner Maximum distance from distortion center to a corner
 * \param [in] a Lens distortion model polynom
 * \param [in] Na Degree of the lens distortion model polynom
 * \author Luis Alvarez, Pedro Henriquez
 */
int build_l1r_vector(std::vector<double> &l1r, 
										 double max_distance_corner,int Na, double *a)
{
	//BUILD INTERMEDIATE VECTOR WITH L-1(r) FROM CENTER TO FURTHEST CORNER
	if(a[Na]==0.) return(-1);
	l1r.resize((int)(max_distance_corner+1.5));

	// AUXILIARY VARIABLES
	double *b,*b2,root2,root=1.;

	// WE ALLOCATE MEMORY
	b=(double*)malloc( sizeof(double)*(Na+2) );
	b2=(double*)malloc( sizeof(double)*(Na+2) );

	// WE DEFINE THE POLYNOM WE NEED TO COMPUTE ROOTS AND THE DERIVATIVE */
	for(int i=1; i<(Na+2); i++){ b[i]=a[i-1];}
	for(int i=0; i<(Na+1); i++){ b2[i]=a[i]*(i+1);}

	// WE BUILD THE VECTOR OF INVERSE LENS DISTORTION FUNCTION
	for(int dist=1; dist<(int)(max_distance_corner+1.5); dist++)
	{
		// WE DEFINE THE POLYNOM WE NEED TO COMPUTE ROOTS AND THE DERIVATIVE */
		b[0]=-dist;

		//NEWTON-RAPHSON TO COMPUTE THE POLYNOMIAL ROOT
		for(int k=0;k<10000;k++){
			double pol_eval=ami_polynomial_evaluation(b,Na+1,root);
			double pol_der=ami_polynomial_evaluation(b2,Na,root);
			if(pol_der==0.) break;
			root2=root-pol_eval/pol_der;
			if(fabs(root-root2)<(fabs(root)*1e-12)){
				root=root2;
				//printf("dist=%d   l1r=%1.8lf pol_eval=%1.2e\n",dist,root/dist,pol_eval);
				break;
			}
			root=root2;
		}

		//PUSH RESULT
		l1r[dist]=root/dist;
	}
	free(b);
	free(b2);
	l1r[0]=l1r[1];

	// WE NORMALIZE THE VECTOR IN ORDER THE CORRECTED IMAGE INCLUDE ALL ORIGINAL IMAGE

	// double paso = 1/l1r[l1r.size()-1];//(l1r.size()-1)/l1r[l1r.size()-1];

	// for(int k=0;k<l1r.size();k++) l1r[k]*=paso;


	return 0;
}

//------------------------------------------------------------------------------

/**
 * \fn int build_l1r_quotient_vector(vector<double> &l1r,point2d<double> &dc,
                                     double max_distance_corner,int Na,
                                     double *a)
 * \brief Build an intermediate vector with values of L-1(r) for
          d = (0 to max_distance_corner)
 * \param [in] [out] l1r vector to store the roots
 * \param [in] dc distortion center
 * \param [in] max_distance_corner Maximum distance from distortion center to a
                corner
 * \param [in] a Lens distortion model polynom
 * \param [in] Na Degree of the lens distortion model polynom
 * \author Luis Alvarez, Pedro Henriquez
 */
int build_l1r_quotient_vector(std::vector<double> &l1r,
                              double max_distance_corner, double *a,int Na)
{
    //BUILD INTERMEDIATE VECTOR WITH L-1(r) FROM CENTER TO FURTHEST CORNER
    if(a[Na]==0. || Na<2) return(-1);
    l1r.resize((int)(max_distance_corner+1.5));

    // WE BUILD THE VECTOR OF INVERSE LENS DISTORTION FUNCTION
    double root,disc,root1,root2;
    int dist=(int)(max_distance_corner+0.5);
    disc= 1.-4*dist*dist*a[2];
    if(disc<0) return(-2);
    disc=sqrt(disc);
    root1=(1-disc)/(2*dist*a[2]);
    root2=(1+disc)/(2*dist*a[2]);
    if(root1>0 ) { l1r[dist]=root1/dist; root=root1;}
    else {l1r[dist]=root2/dist; root=root2; }

    while( (dist--)>0)
    {
        disc= 1.-4*dist*dist*a[2];
        if(disc<0) return(-2);
        disc=sqrt(disc);
        root1=(1-disc)/(2*dist*a[2]);
        root2=(1+disc)/(2*dist*a[2]);
        if(fabs(root-root1)<fabs(root-root2) ) { l1r[dist]=root1/dist; root=root1;}
        else {l1r[dist]=root2/dist; root=root2; }
    }

    l1r[0]=l1r[1];

    return 0;
}
