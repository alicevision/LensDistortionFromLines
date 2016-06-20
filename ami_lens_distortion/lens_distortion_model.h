/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


 #ifndef AMI_DLL_H
  #define AMI_DLL_H
#endif

/**
 * \file lens_distortion_model.h
 * \brief Lens distortion model class AMI_DLL_H definition
 * \author Luis Alvarez \n \n
*/
#ifndef LENS_DISTORTION_MODEL_H
#define LENS_DISTORTION_MODEL_H

#define POLYNOMIAL 0
#define DIVISION 1

#include "../ami_primitives/point2d.h"

#ifdef __AMIDEBUG__
  #include "wxAmiDebugLog.h"
#endif

using namespace ami;
using namespace std;

namespace ami
{
/**
 * \class  lens_distortion_model
 * \brief class  to store distortion model and basic methods
 * \author Luis Alvarez
 */

class AMI_DLL_H  lens_distortion_model{
 std::vector<double> d /** RADIAL DISTORSION MODEL POLYNOMIAL */;
 ami::point2d<double> c /** CENTER OF THE DISTORSION MODEL */;

 int type; /** MODEL TYPE (==0 MEANS POLYNOMIAL AND ==1 MEANS DIVISION) */

public:
 lens_distortion_model(){c.x=0; c.y=0; type=POLYNOMIAL;}; // CONSTRUCTOR WITHOUT TAKING MEMORY
  ~lens_distortion_model(){d.clear();} /** DESTRUCTOR TO FREE MEMORY */;

  /**
  * \fn std::vector<double> get_d()
	* \brief This function returns the vector with the lens distortion model parameters
  * \author Luis Alvarez
  */
 std::vector<double> &get_d(){return d;}


  /**
  * \fn const std::vector<double> &get_d()
  * \author Luis Alvarez
  */
 const std::vector<double> &get_d() const {return d;}

  /**
  * \fn void set_d(const std::vector<double> &d2)
	*	\brief Method for setting the lens distortion model parameters
  * \author Luis Alvarez
  */
 void set_d(const std::vector<double> &d_2){d=d_2;}
 void set_d(std::vector<double> &d_2){d=d_2;}

  /**
  * \fn void set_distortion_center(const point2d<double> &c2)
	* \brief Method for setting the center of the lens distortion model
  * \author Luis Alvarez
  */
  void set_distortion_center(const point2d<double> &c2)/**  SET DISTORTION MODEL CENTER*/{c=c2;}

 /**
  * \fn point2d get_distortion_center()
	* \brief This function returns the center of the lens distortion model
  * \author Luis Alvarez
  */

 point2d<double> &get_distortion_center(){return c;}

  /**
  * \fn const point2d &get_distortion_center()
  * \author Luis Alvarez
  */
 const point2d<double> &get_distortion_center() const {return c;}
 
  /**
  * \fn point2d<double> evaluation_quotient(point2d<double> &p)
	* \brief For a given point, this function computes the division model in that 
	*        point and returns its new position
  * \author Luis Alvarez
  */
 point2d<double> evaluation_quotient(const point2d<double> &p) const;
 
  /**
  * \fn point2d<double> evaluation(point2d<double> &p)
	* \brief For a given point, this function computes the polynomial model in that
	*        point and returns its new position
  * \author Luis Alvarez
  */
 point2d<double> evaluation(const point2d<double> &p) const;
 point2d<double> evaluation(point2d<double> &p);
 
 /**
  * \fn point2d<double> inverse_evaluation(point2d<double> &p)
	* \brief This method computes the inverse polynomial model in a point and returns
	*        its original position
  * \author Luis Alvarez
  */
 point2d<double> inverse_evaluation(const point2d<double> &p) const;
 point2d<double> inverse_evaluation(point2d<double> &p);
 
 /**
  * \fn point2d<double> inverse_evaluation_quotient(point2d<double> &p)
	* \brief This method computes the inverse division model in a point and returns
	*        its original position
  * \author Luis Alvarez
  */
 point2d<double> inverse_evaluation_quotient(const point2d<double> &p) const;
 point2d<double> inverse_evaluation_quotient(point2d<double> &p);

 
  /**
  * \fn point2d<double> inverse_evaluation_fast(point2d<double> &p,double dl1r, double *a, int Na)
	* \brief Accelerated version of the method that computes the inverse polynomial 
	*        model and returns the original position of the point
  * \author Luis Alvarez
  */
 point2d<double> inverse_evaluation_fast(point2d<double> &p,double dl1r, double *a, int Na);
 
 /**
  * \fn std::vector<point2d<double> > evaluation(const std::vector<point2d<double> > & ptl)
	* \brief For a given vector of points, this function computes the polynomial
	*        model of each one and returns their new position inside a vector
  * \author Luis Alvarez
  */
 std::vector<point2d<double> > evaluation(const std::vector<point2d<double> > & ptl) const{
 std::vector<point2d<double> > res;
 for(int i=0;i<((int)ptl.size());i++)
 {
         res.push_back(evaluation(ptl[i]));
 }
 return res;
 }

 void set_type(const int type2){type=type2;}
 int get_type(){return(type);}
 
  /**
  * \fn int read(char name[300])
	* \brief This function reads the lens distortion model from a file
  * \author Luis Alvarez
  */
 int read(char name[300]);
 
  /**
  * \fn write(char name[300])
	* \brief This function writes the lens distortion model to a file
  * \author Luis Alvarez
  */
 int write(char name[300]);

  /**
  * \fn void reset()
  * \brief Reset to initial state
  * \author Pedro Henriquez
  */

 void reset(){d.clear();  c.x=0; c.y=0; /*H.reset();*/}

};
}

#endif
