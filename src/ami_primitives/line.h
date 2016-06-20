/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#ifndef AMI_DLL_H
  #define AMI_DLL_H
#endif

/**
 * \file line.h
 * \brief class  for straight line
 * \author Luis Alvarez \n \n
*/
#ifndef LINE_H
#define LINE_H
#include <iostream>
#include <math.h>
#include "point2d.h"

namespace ami
{
/**
 * \class  line
 * \brief class  for straight lines
 */
// class AMI_DLL_H FOR STRAIGHT LINES
class AMI_DLL_H  line{
 double a,b,c/** LINE EQUATIONS (ax+by+c=0)*/;
public:
 line():a(0.0),b(0.0),c(0.0){}; /** CONSTRUCTOR TO INITIALIZE THE LINE*/
 line(double xx, double yy,double zz):a(xx), b(yy), c(zz){}
 line(const ami::point2d<double> &u,const ami::point2d<double> &v){
   a=v.y-u.y; b=u.x-v.x; c=v.x*u.y-u.x*v.y;}

 /**
  * \fn const double &get_a()
  * \brief Get coefficient "a"
  * \author Luis Alvarez
  */
 const double &get_a() const {return a;}
 double &get_a() {return a;}
 /**
  * \fn const double &get_b()
  * \brief Get coefficient "b"
  * \author Luis Alvarez
  */
 const double &get_b() const {return b;}
 double &get_b() {return b;}
 /**
  * \fn const double &get_c()
  * \brief Get coefficient "c"
  * \author Luis Alvarez
  */
 const double &get_c()  const {return c;}
 double &get_c() {return c;}

 /**
  * \fn get_abc(double &a2,double &b2,double &c2)
  * \brief Get the three coefficients "a,b,c"
  * \author Luis Alvarez
  */
 void get_abc(double &a2,double &b2,double &c2) {a2=a; b2=b; c2=c;}
 void get_abc(double &a2,double &b2,double &c2) const {a2=a; b2=b; c2=c;}


  /**
  * \fn void set_a(double a2)
  * \brief Set coefficient "a"
  * \author Luis Alvarez
  */
 void set_a(double a2){a=a2;}
 /**
  * \fn void set_b(double b2)
  * \brief Set coefficient "b"
  * \author Luis Alvarez
  */
 void set_b(double b2){b=b2;}
  /**
  * \fn void set_c(double c2)
  * \brief Set coefficient "c"
  * \author Luis Alvarez
  */
 void set_c(double c2){c=c2;}
 /**
  * \fn void set_abc(double a2,double b2,double c2)
  * \brief Set the three coefficients "a,b,c"
  * \author Luis Alvarez
  */
 void set_abc(double a2,double b2,double c2){a=a2; b=b2; c=c2;}


 /**
  * \fn double distance(const point2d<double> &p) const
  * \brief Calculates the distance between a point and the line.
  * \author Luis Alvarez
  */
  //template <class  T>
 double distance(const point2d<double> &p) const {
        double d = sqrt(a*a+b*b);
        return (a*p.x+b*p.y+c)/d;
 };

 /**
  * \fn double evaluation(const point2d<double> &p) const
  * \brief Evaluates a point with the line equation.
  * \author Luis Alvarez
  */
 //template <class  T>
 double evaluation(const point2d<double> &p) const {
        return (a*p.x+b*p.y+c);
 };

 friend point2d <double> line_intersection(const line &l1,const line &l2);

};


 point2d <double> line_intersection(const line &l1,const line &l2);

}

 std::ostream &operator <<(std::ostream &s, const ami::line &p);
 std::istream &operator >>(std::istream &s, ami::line &p) ;

#endif
