/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#ifndef AMI_DLL_CPP
  #define AMI_DLL_CPP
#endif

/**
 * \file subpixel_image_contours.cpp
 * \brief Subpixel image contour basic methods
 * \author Luis Alvarez \n \n
 */

#include "subpixel_image_contours.h"
#include "point2d.h"
#include <math.h>

using namespace ami;

/**
 * \fn  subpixel_image_contours::subpixel_image_contours(int width_c,
                                                         int height_c)
 * \brief Constructor taking memory
 * \author Luis Alvarez
 */
AMI_DLL_CPP subpixel_image_contours::
  subpixel_image_contours(int width_c/** New subpixel width */,
                          int height_c/** New subpixel height */)
{
  width   = width_c;
  height  = height_c;
  c       = new bool[width * height];
  x       = new float[width * height];
  y       = new float[width * height];
  d       = new float[width * height];
  coseno  = new float[width * height];
  seno    = new float[width * height];
 }

/**
 * \fn  subpixel_image_contours::
          subpixel_image_contours(const subpixel_image_contours &subpixel)
 * \brief Copy constructor
 */
AMI_DLL_CPP subpixel_image_contours::
  subpixel_image_contours(const subpixel_image_contours &other)
{
  operator=(other);
}

AMI_DLL_CPP subpixel_image_contours& subpixel_image_contours::
  operator=(const subpixel_image_contours &other)
{
  //FREE MEMORY
  delete[] c;
  delete[] x;
  delete[] y;
  delete[] d;
  delete[] coseno;
  delete[] seno;
  //TAKE MEMORY AND COPY CURRENT VALUES
  N       = other.N;
  width   = other.width;
  height  = other.height;
  index   = other.index;
  int lim = width * height;
  c      = new bool [lim];
  x      = new float[lim];
  y      = new float[lim];
  d      = new float[lim];
  coseno = new float[lim];
  seno   = new float[lim];
  //COPY CURRENT VALUES
  for (int i = 0; i < lim; i++)
  {
    c[i]      = other.c[i];
    x[i]      = other.x[i];
    y[i]      = other.y[i];
    d[i]      = other.d[i];
    coseno[i] = other.coseno[i];
    seno[i]   = other.seno[i];
  }
  return *this;
}

/**
 * \fn subpixel_image_contours::~subpixel_image_contours()
 * \brief Destructor to free memory
 */
AMI_DLL_CPP subpixel_image_contours::~subpixel_image_contours()
{
  delete[] c;
  delete[] x;
  delete[] y;
  delete[] d;
  delete[] coseno;
  delete[] seno;
}

/**
 * \fn bool subpixel_image_contours::subpixel_empty()
 * \brief Determine if a subpixel_image_contours is empty
 * \author Pedro Henríquez, Carlos Falcon
 */
AMI_DLL_CPP bool subpixel_image_contours::subpixel_empty()
{
  if( get_width() == 0 && get_height() == 0 )
    return false;
  else
    return true;
}


/**
 * \fn point2d<double> subpixel_image_contours::
        find_nearest_subpixel(point2d<double>  point )
 * \brief Find the nearest point to point px,py in subpixel image contours.
 * \param [in] point coordinates of selected point
 * \author Pedro Henríquez, Carlos Falcon
 */
AMI_DLL_CPP point2d<double> subpixel_image_contours::
  find_nearest_subpixel(point2d<double>  point )
{
  float distance = 99999;
  float dist;
  point2d<double> coordenadas;
  point2d<double> salida;
  salida.x = -1;
  salida.y = -1;
  int size = width * height;
  for(int i = 0;i < size ; i++)
  {
    if(c[i]==1)
    {
      coordenadas.x = i % width;
      coordenadas.y = i / width;
      dist = (point-coordenadas).norm2();
      if(dist < distance)
      {
         distance = dist;
         salida.x = coordenadas.x;
         salida.y = coordenadas.y;
      }
    }
  }
  return salida;
}

/**
 * \fn void subpixel_image_contours::clean(const int neigborhood_radius,
																					 const int min_neigbor_points,
																					 const double min_orientation_value,
																					 const int min_distance_point);
 * \brief Remove outlier contour points
 * \author Luis Alvarez
 */
void subpixel_image_contours::clean(
const int neighborhood_radius, /** radius of neighborhood to take into account */
const int min_neighbor_points, /** min number of contour points in a neighborhood */
const double min_orientation_value, /** min average scalar product of neighborhood point orientation */
const int min_distance_point) /** minimum distance between contour points */
{
  // WE CHECK THAT THE CONTOURS ARE NOT EMPTY
  if(c==NULL || x==NULL || y==NULL || coseno==NULL || seno==NULL)
    return;

  // WE CHECK PARAMETER neighborhood_radius
  if(neighborhood_radius <= 0)
    return;

  int Nedges=index.size();
  if(Nedges == 0)
    return;

  // AUXILIARY VARIABLES
  vector<double> scalar_product(width * height, 0.);
  vector<int> number_neighborhood_points(width * height, 0);
  int i_min = neighborhood_radius;
  int i_max = width - neighborhood_radius;
  int j_min = neighborhood_radius;
  int j_max = height - neighborhood_radius;

  // WE COMPUTE THE NUMBER OF NEIGHBORHOOD POINTS,
  // THE NEIGBORHOOD SCALAR PRODUCT AND WE FILL index VECTOR
  int p=0;
  for(int j=j_min; j<j_max; j++)
  {
    int mj = j * width;
    for(int i=i_min; i<i_max; i++)
    {
      int m = mj+i;
      if(c[m] == 0) continue;
      index[p++] = m;
      for(int k=-neighborhood_radius; k<=neighborhood_radius; k++)
      {
        int k2 = m+k*width;
        for(int l=-neighborhood_radius; l<=neighborhood_radius; l++)
        {
          int n = k2+l;
          if(c[n]==0 || n==m) continue;
          scalar_product[m] += coseno[n]*coseno[m]+seno[n]*seno[m];
          number_neighborhood_points[m]++;
        }
      }
    }
  }
  // WE RESIZE index VECTOR
  index.resize(p);

  // WE REMOVE CONTOUR POINTS ACCORDING TO THE NUMBER OF NEIGHBORHOOD POINTS
  // AND STABILITY OF POINT ORIENTATION (CORNERS ARE REMOVED)
  for(int Np=index.size(), q=0; q<Np; q++){
    int m=index[q];
    if(number_neighborhood_points[m] < min_neighbor_points)
    {
      c[m]=0;
      continue;
    }
    if(scalar_product[m] < (number_neighborhood_points[m]*min_orientation_value)){
      c[m]=0;
      continue;
    }
  }

  // WE REMOVE ISOLATED POINTS USING AND ITERATIVE PROCEDURE
  for(int s=0; s<4; s++)
  {
    // WE COMPUTE THE NUMBER OF POINTS IN THE NEIGHBORHOOD
    for(int Np=index.size(), q=0; q<Np; q++)
    {
      int m=index[q];
      if(c[m] == 0) 
        continue;
      number_neighborhood_points[m]=0;
      for(int k=-neighborhood_radius; k<=neighborhood_radius; k++)
      {
        int k2 = m+k*width;
        for(int l=-neighborhood_radius;l<=neighborhood_radius;l++)
        {
          int n = k2+l;
          if(c[n]==0 || n==m) 
            continue;
          number_neighborhood_points[m]++;
        }
      }
    }
    // WE REMOVE ISOLATED CONTOUR POINTS
    for(int Np=index.size(), q=0; q<Np; q++)
    {
      int m=index[q];
      if(c[m] == 0) 
        continue;
      if(number_neighborhood_points[m] < neighborhood_radius)
      {
       c[m]=0;
      }
    }
  }

  // WE KEEP ONLY 1 CONTOUR POINT EN EACH WINDOW OF RADIUS min_distance_point
  if(min_distance_point > 0)
  {
    i_min = min_distance_point;
    i_max = width-min_distance_point;
    j_min = min_distance_point;
    j_max = height-min_distance_point;
    int window_size = 2*min_distance_point+1;
    for(int j=j_min; j<j_max; j+=window_size)
    {
      int mj=j*width;
      for(int i=i_min; i<i_max; i+=window_size)
      {
        int m=mj+i;
        double max=scalar_product[m];
        int m_max=m;
        for(int k=-min_distance_point; k<=min_distance_point; k++)
        {
          int k2=m+k*width;
          for(int l=-min_distance_point; l<=min_distance_point; l++)
          {
            int n=k2+l;
            if(c[n]==0 || m==n) 
              continue;
            if(scalar_product[n]>max)
            {
              max=scalar_product[n];
              c[m_max]=0;
              m_max=n;
            }
            else 
            { //if(scalar_product[n]<max){
             c[n]=0;
            }
          }
        }
      }
    }
  }

  // WE FILL INDEX VECTOR
  int m=0;
  for(int Np=index.size(), q=0; q<Np; q++)
  {
    if(c[index[q]]>0) 
      index[m++]=index[q];
  }
  index.resize(m);
  return;
}

/**
 * \fn subpixel_image_contours::build_index()
 * \brief build index vector
 * \author Luis Alvarez
 */
void subpixel_image_contours::build_index()
{
  int m=0,size=width*height;
  index.resize(size);
  for(int k=0; k<size; k++)
  {
    if(c[k]>0) 
      index[m++]=k;
  }
  index.resize(m);
}
