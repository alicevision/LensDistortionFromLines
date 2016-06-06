/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#ifndef AMI_DLL_H
  #define AMI_DLL_H
#endif

/**
 * \file point2d.h
 * \brief point2d class  definition
 * \author Luis Alvarez \n \n
*/
#ifndef POINT2D_H
#define POINT2D_H
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

namespace ami
{

/**
 * \class  point2d
 * \brief class  to store 2D points and basic methods
 * \author Luis Alvarez
 */
template < class  T >
class AMI_DLL_H point2d {

  public :
  T x /** point x coordinate */;
  T y /** point y coordinate */;
	point2d():x( (T) 0), y( (T) 0){};
	~point2d(){};
	point2d(const T xx , const T yy){x=xx; y=yy;}
	point2d(const T scalar){x=y=scalar;}
	point2d & operator=(const point2d &p){ x=p.x; y=p.y; return *this;}
	point2d & operator=(const T scalar){ x=scalar; y=scalar; return *this;}
  point2d (const point2d<T> &p){x=p.x; y=p.y;}
	point2d operator+(const point2d &p)const {return point2d(x+p.x,y+p.y);}
	point2d operator-(const point2d &p)const {return point2d(x-p.x,y-p.y);}
	point2d operator*(const T &a)const {return point2d(a*x,a*y);}
	double operator*(const point2d &p)const {return ( (double) x*p.x+y*p.y);}
	inline T norm(){T paso=x*x+y*y; return(paso>0.?sqrtf((float) paso):0.);}
	inline T norm2(){ return(x*x+y*y);}
	void print(){ std::cout << "point2d : (" << x << "," << y << ")" <<
                std::endl;}
	int find_nearest_point(std::vector< point2d <T> > primitive);
	void find_nearest_point(std::vector< point2d <T> > primitive,int &id,
                         T &distance);
	bool operator!=(const T scalar) const {return(x!=scalar && y!=scalar);}
};

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class  T>  int point2d<T>::
                        find_nearest_point(std::vector< point2d <T> > primitive)
 * \brief Find the point in th primitive that is the neareast to point, and
          returns its index
 * \param [in] primitive Vector with points 2d, it can represente a line,
               a polygon,a contour,etc..
 * \author Pedro Henriquez
*/
template <class  T>
  int point2d<T>::find_nearest_point(std::vector< point2d <T> > primitive)
  {
    //FIND THE NEAREST POINT TO THE ACTUAL POINT
    T min_dist,dist;
    int min=-1;
    min_dist=99999;
    point2d<T> point_p;
    for(unsigned int i=0;i<primitive.size();i++)
    {
      point_p.x=primitive[i].x;
      point_p.y=primitive[i].y;
      dist = ((*this)-point_p).norm2(); //DISTANCE
      if(dist < min_dist)
      {
        min_dist=dist;
        min=i;
      }
    }
    return min;
  }


/**
 * \fn template <class  T> void point2d<T>::
  find_nearest_point(std::vector< point2d <T> > primitive,int &id, T &distance)
 * \brief Find the point in th primitive that is the neareast to point, and
          returns its index
 * \param [in] primitive Vector with points 2d, it can represente a line,
               a polygon,a contour,etc..
 * \author Pedro Henriquez
*/
template <class  T>
  void point2d<T>::find_nearest_point(std::vector< point2d <T> >primitive/**.*/,
                                      int &id/**.*/, T &distance/**.*/)
  {
    //FIND THE NEAREST POINT TO THE ACTUAL POINT
    T dist;
    id=-1;
    distance=99999;
    point2d<T> point_p;
    for(unsigned int i=0;i<primitive.size();i++)
    {
      point_p.x=primitive[i].x;
      point_p.y=primitive[i].y;
      dist = ((*this)-point_p).norm2(); //DISTANCE
      if(dist < distance)
      {
        distance=dist;
        id=i;
      }
    }
  }

}
#endif
