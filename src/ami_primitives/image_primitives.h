/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#ifndef AMI_DLL_H
  #define AMI_DLL_H
#endif

/**
 * \file  image_primitives.h
 * \brief class  to store image primitives (points,lines,ellipses and distortion
                                            model)
 * \author Luis Alvarez \n \n
*/
#ifndef IMAGE_PRIMITIVES_H
#define IMAGE_PRIMITIVES_H

#include "../ami_lens_distortion/lens_distortion_model.h"
#include "point2d.h"
#include "line_points.h"

using namespace std;

namespace ami {
/**
 * \class  image_primitives
 * \brief class  to store the image primitives (points, lines,
                                                ellipses and ditortion model)
 * \author Luis Alvarez
 */
// * \test \ref test_image_primitives_test.cpp "image_primitives_test.cpp"

class AMI_DLL_H  image_primitives
{
 //homography H /** HOMOGRAPHY WITH RESPECT TO THE REFERENCE VIEW*/;
 std::vector< point2d<double> > points /**  REFERENCE POINTS (CORNERS OF SOCCER
                                                            STADIUM LOCATION)*/;
 std::vector<line_points> lines /** LINES AND ASSOCIATED POINTS */;
// std::vector<ellipse_points> ellipses /** ELLIPSES AND ASSOCIATED POINTS */;
 lens_distortion_model distortion /** LENS DISTORTION MODEL */;
public :

  image_primitives(){};
  image_primitives(double *a,double *b, double *c,int size);

 std::vector< point2d<double> > projected_3dpoints /**  COORDINATES OF PROJECTED
                                   3D POINTS LOCATED OUTSIDE PRIMITIVE PLANE */;

 image_primitives & operator=(const image_primitives &primitiva)
 {
  points     = primitiva.points;
  lines      = primitiva.lines;
  distortion = primitiva.distortion;
  return *this;
 }

 ~image_primitives(){};

 void initialize(unsigned int npoints=0, unsigned int nlines=0,
                 unsigned int nellipses=0);

 bool IsProjected_3dpointsInvalid ();



/**
  * \fn homography get_homography()
  * \brief Return Homography
  * \author Luis Alvarez
  */
 //const homography &get_homography(void) const {return H;}
 //homography &get_homography(void)  {return H;}




/**
  * \fn const point2d<double> get_distorsion_center()const
  * \brief Return distortion center
  * \author Carlos Falcon
  */
  const point2d<double> get_distorsion_center()const {
    return distortion.get_distortion_center();}


/**
  * \fn std::vector<point2d> &get_points()
  * \brief Return points
  * \author Luis Alvarez
  */
 std::vector<point2d<double> > &get_points(){return points;}

 const std::vector<point2d<double> > &get_points() const {return points;}

/**
  * \fn std::vector<line_points> &get_lines()
  * \brief Return line_points collection
  * \author Luis Alvarez
  */
 std::vector<line_points> &get_lines(){return lines;}
/**
  * \fn const std::vector<line_points> &get_lines()
  * \brief Return lines_poinst collection
  */
 const std::vector<line_points> &get_lines() const {return lines;}


/**
  * \fn std::vector<ellipse> &get_ellipses()
  * \brief Return ellipse_points collection
  * \author Luis Alvarez
  */
 //std::vector<ellipse_points> &get_ellipses(){return ellipses;}
/**
  * \fn const std::vector<ellipse_points> &get_ellipses()
  * \brief Return ellipse_points collection.
  * \author Luis Alvarez
  */
 //const std::vector<ellipse_points> &get_ellipses() const {return ellipses;}


/**
  * \fn lens_distortion_model get_distortion()
  * \brief Return distortion model.
  * \author Luis Alvarez
  */
 const lens_distortion_model &get_distortion() const {return distortion;}
 lens_distortion_model &get_distortion(){return distortion;}

/**
  * \fn void set_homography(const homography &H2)
  * \brief Set homography.
  * \author Luis Alvarez
  */
 //void set_homography(const homography &H2/** Homography */){H=H2;}

/**
  * \fn void set_points(const std::vector<point2d<double> > &p2)
  * \brief Set points collection.
  * \author Luis Alvarez
  */
 void set_points(const std::vector<point2d<double> > &p2/** point2d */){
   if(points.size()>0) points.clear(); points=p2;}

/**
  * \fn void set_lines(const std::vector<line_points> &lines2)
  * \brief Set Lines collection.
  * \author Luis Alvarez
  */
 void set_lines(const std::vector<line_points> &lines2 /**collection of lines*/)
  {if(lines.size()>0) lines.clear();lines=lines2;}

/**
  * \fn void set_ellipses(const std::vector<ellipse_points> &ellipses2)
  * \brief Set ellipse collection.
  * \author Luis Alvarez
  */
// void set_ellipses(const std::vector<ellipse_points> &ellipses2/** collection
//                     of ellipse points */){
//      if(ellipses.size()>0) ellipses.clear(); ellipses=ellipses2;}

/**
  * \fn void set_distortion(const lens_distortion_model &distortion2)
  * \brief Set distortion model.
  * \author Luis Alvarez
  */
 void set_distortion(const lens_distortion_model &distortion2/** Distortion
                     model */) {distortion=distortion2;}

 int read (char *name);

 int write(const std::string& name) { return write(name.c_str()); }
 int write(const char *name/**Filename*/)
{
  FILE *f;
  double v1,v2,v3;

  //OPEN FILE
  f = fopen(name, "w+");
  if (!f)
     return -1;

   //SAVE DISTORTION MODEL
  fprintf(f,"DISTORTION MODEL\n");
  fprintf(f,"distortion center = (%.20g %.20g)\n",distortion.get_distortion_center().x,distortion.get_distortion_center().y);
  fprintf(f,"Number of distortion parameters = %d\n",(int)distortion.get_d().size());
  for(unsigned int i=0;i<distortion.get_d().size();i++)
  {
    fprintf(f, "value of distortion parameter %d  = %.20g\n",i,distortion.get_d()[i]);
  }
  //SAVE POINTS
  for(unsigned int i=0;i<points.size();i++)
  {
          fprintf(f, "%.20g %.20g\n",points[i].x,points[i].y);
  }

  //SAVE STRAIGHT
  fprintf(f,"\nNUMBER OF LINES = %d \n",(int)lines.size());
  for(unsigned int i = 0; i < lines.size(); i++)
  {
    fprintf(f,"\nline %d \n",(int) i);
    v1 = lines[i].get_a();
    v2 = lines[i].get_b();
    v3 = lines[i].get_c();
    fprintf(f,"line parameters = %.20g %.20g %.20g\n",v1,v2,v3);
    fprintf(f,"line points : \n");
    const std::vector<point2d<double> >& puntos = lines[i].get_points();
    fprintf(f,"%d\n",(int)puntos.size());
    for(unsigned int j=0;j<puntos.size();j++)
    {
        fprintf(f, "%.20g %.20g\n",puntos[j].x,puntos[j].y);
    }
  }

  /* WE STORE 3D POINTS PROJECTED COORDINATES*/
  //fprintf(f,"\n%d 3D_POINTS_PROJECTED_COORDINATES \n",(int)projected_3dpoints.size());
  /*for(unsigned int i=0;i<projected_3dpoints.size();i++)
  {
          fprintf(f, "%.15f %.15f\n",projected_3dpoints[i].x,projected_3dpoints[i].y);
  }*/

  //CLOSE FILE
  fclose(f);

  return 0;
}

 unsigned int Num_Puntos() const;

friend inline std::ostream &operator <<(std::ostream &s,
                                        const image_primitives &ip) {
  //double d = /*ip.H.H[2][2]*/0.;
  s << "H " << std::endl;
  for(int i=0;i<3;i++)
  {
          for(int j=0;j<3;j++)
          {
                  //s << ip.H.H[i][j] / ip.H.H[2][2] << " ";
          }
          s << std::endl;
  }
  s << "Points " << ip.points.size() << std::endl;
  for(unsigned int i=0;i<ip.points.size();i++)
  {
          s << "(" << ip.points[i].x << ", " << ip.points[i].y << ") ";
  }
  s << std::endl;
  s << "Lines " << ip.lines.size() << std::endl;
  for(unsigned int i=0; i<ip.lines.size();i++)
  {
          s << "(" << ip.lines[i].get_a() << ", " << ip.lines[i].get_b() <<
          ", " << ip.lines[i].get_c() << ", " << std::endl;
  }
  s << std::endl;
  s << "Distortion model " << std::endl;
  s << "( " << ip.distortion.get_distortion_center().x << ", " <<
  ip.distortion.get_distortion_center().y << ") ";
  for(unsigned int i=0;i<ip.distortion.get_d().size();i++)
  {
          s << ip.distortion.get_d()[i] << " ";
  }
  s << std::endl;
  return s;
}

bool is_empty();
void clear()
{
	points.clear();
  lines.clear();
  distortion.reset();
}

void find_nearest_point_in_primitives(point2d<double> &pto,
                                      point2d<double> &pto_out,double &distance,
                                      int &type_prim,int &prim,int &pos);

void find_nearest_point_in_others_primitives(point2d<double> &pto,
                                             point2d<double> &pto_out,
                                             double &distance, int &type_prim,
                                             int &prim,int &pos);
void find_nearest_line_equation(point2d<double> &pto,double &distance,
                                int &prim);
};

}

#endif
