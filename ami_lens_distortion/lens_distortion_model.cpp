/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */

 #ifndef AMI_DLL_CPP
  #define AMI_DLL_CPP
#endif

/**
 * \file lens_distortion_model.cpp
 * \brief Lens distortion model basic procedures
 * \author Luis Alvarez \n \n
 */

#include "lens_distortion_model.h"
#include "../ami_primitives/line_points.h"
#include "../ami_pol/ami_pol.h"
#include "../ami_primitives/point2d.h"
#include "lens_distortion.h"
#include <vector>
#include <iostream>
#include <cstdio>

using namespace ami;

/**
 * \fn point2d<double> lens_distortion_model::evaluation_quotient(const point2d<double> &p) const
 * \brief COMPUTE THE MODEL IN A POINT (WE CORRECT THE DISTORTION)
 * \param p point to be evaluated
 * \author Luis Alvarez
 */
 // COMPUTE THE MODEL IN A POINT
AMI_DLL_CPP point2d<double> lens_distortion_model::evaluation_quotient(const point2d<double> &p) const
 {
  if(d.size()==0) return(p);
  else{
    double norm=(p.x-c.x)*(p.x-c.x)+(p.y-c.y)*(p.y-c.y);
    double A=ami_polynomial_evaluation(( double*)&d[0],(int) d.size()-1,norm);
    if(A==0) return(p);
    return( point2d<double>( c.x+(p.x-c.x)/A,c.y+(p.y-c.y)/A));
  }
 }

 /**
 * \fn point2d<double> lens_distortion_model::evaluation_quotient(point2d<double> &p)
 * \brief COMPUTE THE MODEL IN A POINT (WE CORRECT THE DISTORTION)
 * \param p point to be evaluated
 * \author Luis Alvarez
 */
 // COMPUTE THE MODEL IN A POINT
AMI_DLL_CPP point2d<double> lens_distortion_model::evaluation_quotient(point2d<double> &p)
 {
  if(d.size()==0) return(p);
  else{
    double norm=(p.x-c.x)*(p.x-c.x)+(p.y-c.y)*(p.y-c.y);
    double A=ami_polynomial_evaluation(( double*)&d[0],(int) d.size()-1,norm);
    if(A==0) return(p);
    return( point2d<double>( c.x+(p.x-c.x)/A,c.y+(p.y-c.y)/A));
  }
 }


/**
 * \fn point2d<double> lens_distortion_model::evaluation(const point2d<double> &p) const
 * \brief COMPUTE THE MODEL IN A POINT (WE CORRECT THE DISTORTION)
 * \param p point to be evaluated
 * \author Luis Alvarez
 */
 // COMPUTE THE MODEL IN A POINT
 AMI_DLL_CPP point2d<double> lens_distortion_model::evaluation(const point2d<double> &p) const
 {

  if(d.size()==0) return(p);
  if(type==DIVISION){
    double norm=(p.x-c.x)*(p.x-c.x)+(p.y-c.y)*(p.y-c.y);
    double A=ami_polynomial_evaluation(( double*)&d[0],(int) d.size()-1,norm);
    if(A==0) return(p);
    return( point2d<double>( c.x+(p.x-c.x)/A,c.y+(p.y-c.y)/A));
  }
  double norm=(p.x-c.x)*(p.x-c.x)+(p.y-c.y)*(p.y-c.y);
  double A=ami_polynomial_evaluation(( double*)&d[0],(int) d.size()-1,norm);
  return( point2d<double>((c.x+(p.x-c.x)*A),(c.y+(p.y-c.y)*A)));
 }
 

/**
 * \fn point2d<double> lens_distortion_model::inverse_evaluation(point2d<double> &p)
 * \brief COMPUTE THE INVERSE MODEL IN A POINT
 * \param p point to be evaluated
 * \author Luis Alvarez
 */
// COMPUTE THE INVERSE MODEL IN A POINT

AMI_DLL_CPP point2d<double> lens_distortion_model::inverse_evaluation(point2d<double> &p){
 if(type==DIVISION){
   return((*this).inverse_evaluation_quotient(p));
 }
 double x,y;
 double *a;
 unsigned int i;
 int Na=2*(d.size()-1);
 if(d.size()==0) return(p);
 a=(double*)malloc(sizeof(double)*(Na+1));
 a[0]=d[0];
 for(i=1;i<d.size();i++){a[2*i-1]=0.; a[2*i]=d[i]; }
 int cont2=Na;
 while(a[cont2]==0){ cont2--; if (cont2 == 0)break;}
 Na=cont2;
 if(Na==0) return(p);

 ami_inverse_lens_distortion(p.x,p.y,c.x,c.y,&x,&y,a,Na);
 free(a);
 return( point2d<double>(x,y));
}

/**
 * \fn point2d<double> lens_distortion_model::inverse_evaluation_quotient(point2d<double> &p)
 * \brief COMPUTE THE INVERSE RATIONAL MODEL IN A POINT
 * \param p point to be evaluated
 * \author Luis Alvarez, Daniel Santana-Cedrés
 */
// COMPUTE THE INVERSE RATIONAL MODEL IN A POINT

AMI_DLL_CPP point2d<double> lens_distortion_model::inverse_evaluation_quotient(point2d<double> &p){
    if(d.size()==2){
      //Evaluated coordinates
      double x,y;
      //Distortion parameter
      double k  = d[1];
      double x_ = p.x-c.x;
      double y_ = p.y-c.y;
      //Radius from the input point to the distortion center
      //double r_ = sqrt(((x_)*(x_))+((y_)*(y_)));
      double r_ = ((x_)*(x_))+((y_)*(y_));
      //Discriminant of the quadratic equation
      double discriminant = 1-4*k*r_;//*r_;
      double denominator = 2*k*r_;//*r_;

      //If the radius is 0, we return the input point
      if(r_==0) return p;

      //We show a message if the discriminant is negative and return the input point
      if(discriminant<0)
      {
          cout << "It's no possible to compute the inverse_evaluation_quotient: discriminant<0" << endl;
          return p;
      }

      //We take only the negative solution of the quadratic equation
      double new_r = (1-sqrt(discriminant))/denominator;

      //Values of the coordinates
      x = c.x + new_r * x_;
      y = c.y + new_r * y_;
      return point2d<double>(x,y);
    }
    else{
      //Evaluated coordinates
      double x,y;
      //Distortion parameter
      double x_ = p.x-c.x;
      double y_ = p.y-c.y;
      //Radius from the input point to the distortion center
      //double r_ = sqrt(((x_)*(x_))+((y_)*(y_)));
      double r_ = ((x_)*(x_))+((y_)*(y_));
      if(r_<=0.) return(p);
      double norm=sqrt(r_);

      int i,Nr,degree=2*d.size()-2;
      double *b,*rx,*ry,root;
      /* WE ALLOCATE MEMORY */
      b=(double*)malloc( sizeof(double)*(degree+1) );
      rx=(double*)malloc( sizeof(double)*(degree+1) );
      ry=(double*)malloc( sizeof(double)*(degree+1) );
      /* WE DEFINE THE POLYNOMIAL WE NEED TO COMPUTE ROOTS */
      b[degree]=norm*d[0]; b[degree-1]=-1;
      int m=1;
      for(i=degree-2;i>=0;i--){
        if(i%2==1) b[i]=0;
        else {
         //printf("i=%d, d[%d]=%e\n",i,m,d[m]);
         b[i]=norm*d[m++];
        }

      }
      //for(i=0;i<=degree;i++) printf("b[%d]=%e\n",i,b[i]);

      Nr=ami_polynomial_root(b,degree,rx,ry);
      //printf("degree=%d Nr=%d\n",degree,Nr);
      //system("pause");
      /* WE SELECT THE smaller REAL ROOT */
      root=1e30;
      for(i=0;i<Nr;i++){
        //printf("rx[%d]=%e ry[%d]=%e\n",i,rx[i],i,ry[i]);
        //if(fabs(ry[i])<0.00000000001 && fabs(root-norm)>fabs(rx[i]-norm)) root=rx[i];
        if(fabs(ry[i])<0.00000000001 && fabs(root)>fabs(rx[i])) root=rx[i];
      }
      if(Nr==0){
        //We show a message if the discriminant is negative and return the input point
        cout << "It's no possible to compute the inverse_evaluation_quotient: Number of polynomial roots == 0" << endl;
        //system("pause");
        //return point2d<double>(0,0);
        return p;
      }

      //printf("root=%lf,norm=%lf,norm2=%lf\n",root,norm,norm2);

      /* WE TRANSFORM THE POINT COORDINATES */
      x = c.x + root * x_ / norm;
      y = c.y + root * y_ / norm;
      free(b); free(rx); free(ry);
      return point2d<double>(x,y);
    }
}

/**
 * \fn point2d<double> lens_distortion_model::inverse_evaluation_fast(point2d<double> &p,double dl1r, double *a, int Na)
 * \brief COMPUTE THE INVERSE MODEL IN A POINT
 * \param p point to be evaluated
 *  \param [in] dl1r coeficient interpolated from vector of max distorsion distances
 *  \param [in] a Lens distortion model polynom
 *  \param [in] Na Degree of the lens distortion model polynom
 * \author Luis Alvarez, Pedro Henriquez
 */
// COMPUTE THE INVERSE MODEL IN A POINT

point2d<double> lens_distortion_model::inverse_evaluation_fast(point2d<double> &p,double dl1r, double *a, int Na){
 double x,y;
 ami_inverse_lens_distortion_fast(p.x,p.y,c.x,c.y,&x,&y,a,Na,dl1r);
 return( point2d<double>(x,y));
}

/**
 * \fn point2d<double> lens_distortion_model::evaluation(point2d<double> &p)
 * \brief COMPUTE THE MODEL IN A POINT
 * \param p
 * \author Luis Alvarez
 */
 // COMPUTE THE MODEL IN A POINT
AMI_DLL_CPP point2d<double> lens_distortion_model::evaluation(point2d<double> &p)
 {

  if(d.size()==0) return(p);
  if(type==DIVISION){
    double norm=(p.x-c.x)*(p.x-c.x)+(p.y-c.y)*(p.y-c.y);
    double A=ami_polynomial_evaluation(( double*)&d[0],(int) d.size()-1,norm);
    if(A==0) return(p);
    return( point2d<double>( c.x+(p.x-c.x)/A,c.y+(p.y-c.y)/A));
  }
  double norm=(p.x-c.x)*(p.x-c.x)+(p.y-c.y)*(p.y-c.y);
  double A=ami_polynomial_evaluation(( double*)&d[0],(int) d.size()-1,norm);
  return( point2d<double>((c.x+(p.x-c.x)*A),(c.y+(p.y-c.y)*A)));
 }

/**
 * \fn point2d<double> lens_distortion_model::inverse_evaluation(const point2d<double> &p) const
 * \brief COMPUTE THE INVERSE MODEL IN A POINT
 * \param p
 * \author Luis Alvarez
 */
// COMPUTE THE INVERSE MODEL IN A POINT

AMI_DLL_CPP point2d<double> lens_distortion_model::inverse_evaluation(const point2d<double> &p) const{
 if(type==DIVISION){
   return((*this).inverse_evaluation_quotient(p));
 }
 double x,y;
 double *a;
 unsigned int i;
 int Na=2*(d.size()-1);
 if(d.size()==0) return(p);
 a=(double*)malloc(sizeof(double)*(Na+1));
 a[0]=d[0];
 for(i=1;i<d.size();i++){a[2*i-1]=0.; a[2*i]=d[i]; }
 ami_inverse_lens_distortion(p.x,p.y,c.x,c.y,&x,&y,a,Na);
 free(a);
 return( point2d<double>(x,y));
}

/**
 * \fn point2d<double> lens_distortion_model::inverse_evaluation_quotient(const point2d<double> &p) const
 * \brief COMPUTE THE INVERSE RATIONAL MODEL IN A POINT
 * \param p
 * \author Luis Alvarez, Daniel Santana-Cedrés
 */
// COMPUTE THE INVERSE MODEL IN A POINT

AMI_DLL_CPP point2d<double> lens_distortion_model::inverse_evaluation_quotient(const point2d<double> &p) const{
   if(d.size()==2){
      //Evaluated coordinates
      double x,y;
      //Distortion parameter
      double k  = d[1];
      double x_ = p.x-c.x;
      double y_ = p.y-c.y;
      //Radius from the input point to the distortion center
      //double r_ = sqrt(((x_)*(x_))+((y_)*(y_)));
      double r_ = ((x_)*(x_))+((y_)*(y_));
      //Discriminant of the quadratic equation
      double discriminant = 1-4*k*r_;//*r_;
      double denominator = 2*k*r_;//*r_;

      //If the radius is 0, we return the input point
      if(r_==0) return p;

      //We show a message if the discriminant is negative and return the input point
      if(discriminant<0)
      {
          cout << "It's no possible to compute the inverse_evaluation_quotient: discriminant<0" << endl;
          return p;
      }

      //We take only the negative solution of the quadratic equation
      double new_r = (1-sqrt(discriminant))/denominator;

      //Values of the coordinates
      x = c.x + new_r * x_;
      y = c.y + new_r * y_;
      return point2d<double>(x,y);
    }
    else{
      //Evaluated coordinates
      double x,y;
      //Distortion parameter
      double x_ = p.x-c.x;
      double y_ = p.y-c.y;
      //Radius from the input point to the distortion center
      //double r_ = sqrt(((x_)*(x_))+((y_)*(y_)));
      double r_ = ((x_)*(x_))+((y_)*(y_));
      if(r_<=0.) return(p);
      double norm=sqrt(r_);

      int i,Nr,degree=2*d.size()-2;
      double *b,*rx,*ry,root;
      /* WE ALLOCATE MEMORY */
      b=(double*)malloc( sizeof(double)*(degree+1) );
      rx=(double*)malloc( sizeof(double)*(degree+1) );
      ry=(double*)malloc( sizeof(double)*(degree+1) );
      /* WE DEFINE THE POLYNOMIAL WE NEED TO COMPUTE ROOTS */
      b[degree]=norm*d[0]; b[degree-1]=-1;
      int m=1;
      for(i=degree-2;i>=0;i--){
        if(i%2==1) b[i]=0;
        else {
         //printf("i=%d, d[%d]=%e\n",i,m,d[m]);
         b[i]=norm*d[m++];
        }

      }
      //for(i=0;i<=degree;i++) printf("b[%d]=%e\n",i,b[i]);

      Nr=ami_polynomial_root(b,degree,rx,ry);
      //printf("degree=%d Nr=%d\n",degree,Nr);
      /* WE SELECT THE smaller REAL ROOT */
      root=1e30;
      for(i=0;i<Nr;i++){
        //printf("rx[%d]=%e ry[%d]=%e\n",i,rx[i],i,ry[i]);
        //if(fabs(ry[i])<0.00000000001 && fabs(root-norm)>fabs(rx[i]-norm)) root=rx[i];
        if(fabs(ry[i])<0.00000000001 && fabs(root)>fabs(rx[i])) root=rx[i];
      }
      if(Nr==0){
        //We show a message if the discriminant is negative and return the input point
        cout << "It's no possible to compute the inverse_evaluation_quotient: Problem with polynomial roots<0" << endl;
        return p;
      }

      //printf("root=%lf,norm=%lf,norm2=%lf\n",root,norm,norm2);

      /* WE TRANSFORM THE POINT COORDINATES */
      x = c.x + root * x_ / norm;
      y = c.y + root * y_ / norm;
      free(b); free(rx); free(ry);
      return point2d<double>(x,y);
    }
}

/**
 * \fn int lens_distortion_model::read(char name[300])
 * \brief Function to read from disk lens_distortion model information
 * \author Luis Alvarez
 */
AMI_DLL_CPP int lens_distortion_model::read(
  char name[300] /**INPUT FILE NAME */)
{
  int i,N;
  FILE *f;
  char name_paso[300];
  double a;

  /* Open file */
  if ((f = fopen(name,"rb")) == NULL) {
    fprintf(stderr,"Unable to read file \"%s\"\n",name);
    return(-1);
  }

  // WE CLEAR DISTORTION MODEL
  d.clear();

	int val = fscanf(f,"%s\n",name_paso);
  if(val==EOF) printf("Read error: distortion_center\n");

  val = fscanf(f,"%lf %lf \n",&(c.x),&(c.y));
  if(val==EOF) printf("Error reading the distortion center value\n");


  val = fscanf(f,"%s\n",name_paso);
  if(val==EOF) printf("Read error: first_distortion_model\n");

  val = fscanf(f,"%d\n",&N);
  if(val==EOF) printf("Error reading number of components\n");

  for(i=0;i<N;i++){

    val = fscanf(f,"%lf\n",&a);
    if(val==EOF) printf("Error reading lens distortion parameter %d\n",i);
    d.push_back(a);

  }

  fclose(f);
  return 0;
}



/**
 * \fn int lens_distortion_model::write(char name[300])
 * \brief Function to write to disk lens_distortion model information
 * \author Luis Alvarez
 */
AMI_DLL_CPP int lens_distortion_model::write(
  char name[300] /**INPUT FILE NAME */)
{
         int i;
  FILE *f;

  /* Open file */
  if ((f = fopen(name,"wb")) == NULL) {
    fprintf(stderr,"Unable to read file \"%s\"\n",name);
    return(-1);
  }
  fprintf(f,"distortion_center\n");
  fprintf(f,"%lf %lf \n",c.x,c.y);

  fprintf(f,"first_distortion_model\n");
  fprintf(f,"%d\n",(int)d.size());
  for(i=0;i<(int)d.size();i++){
    fprintf(f,"%1.16e\n",d[i]);
  }

  fclose(f);
  return 0;
}


