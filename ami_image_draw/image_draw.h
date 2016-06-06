/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */

#ifndef AMI_DLL_H
  #define AMI_DLL_H
#endif

/**
 * \file image_draw.h
 * \brief image_draw class AMI_DLL_H definition
 * \author Carlos Falcon \n \n
*/


#ifndef image_draw_H
#define image_draw_H

#ifndef NOMINMAX
#define NOMINMAX
#endif


#include <vector>

#include "../ami_image/image.h"


#define image_draw_DEBUG

using namespace std;

namespace ami
{


class AMI_DLL_H image_draw
{
  public :
  template <class  T>
  void draw_cercle(ami::image<T> &imagen, double x0,double y0,float radius, unsigned char red_arrow=255,unsigned char green_arrow=0,unsigned char blue_arrow=0);

};

/**
 * \fn template <class  T> void image_draw::draw_cercle(   ami::image<T> &imagen,   double x0,double y0,   float radius,
   unsigned char red_arrow,unsigned char green_arrow,unsigned char blue_arrow)
 * \brief Function to draw a cercle in an unsigned char image channel
 * \param  imagen INPUT/OUTPUT IMAGE WHERE THE CERCLE WILL BE DRAWN
 * \param [in] width, height IMAGE DIMENSIONS
 * \param [in] x0,y0 LOCATION OF THE CERCLE CENTRE
 * \param [in] radius RADIUS OF THE CERCLE
 * \param [in] red_arrow,green_arrow,blue_arrow COLOR TO FILL THE CIRCLE
 * \author Carlos Falcon
 */
template <class  T>
void image_draw::draw_cercle(
  ami::image<T> &imagen,
  double x0,double y0,
  float radius,
  unsigned char red_arrow,unsigned char green_arrow,unsigned char blue_arrow)

{
  int width=imagen.width();
  int height=imagen.height();
  int size=width*height;
  int i,j,cont,cont2;
  float m,n,paso=(float)0.06;
  // WE DRAW THE CERCLE
  for (int k=0;k<imagen.nChannels();k++)
  {
    unsigned char cercle_color = 0;
    switch(k)
    {
      case 0:
        cercle_color= red_arrow;
        break;
      case 1:
        cercle_color= green_arrow;
        break;
      case 2:
        cercle_color= blue_arrow;
        break;
    }
    if( (x0+radius)<0 || (y0+radius)<0 ) return;
    if( (x0-radius)>width || (y0-radius)>height) return;
    for(i=(int) (y0-radius-1);i<=(y0+radius);i++){
      for(j=(int) (x0-radius-1);j<=(x0+radius);j++){
        if(i<0 || i>=height  || j<0 || j>=width) continue;
        cont=0;
        if( ((i-y0)*(i-y0)+(j-x0)*(j-x0))<(radius*radius) ) cont++;
        if( ((i+1-y0)*(i+1-y0)+(j-x0)*(j-x0))<(radius*radius) ) cont++;
        if( ((i-y0)*(i-y0)+(j+1-x0)*(j+1-x0))<(radius*radius) ) cont++;
        if( ((i+1-y0)*(i+1-y0)+(j+1-x0)*(j+1-x0))<(radius*radius) ) cont++;
        if( cont==4 ) imagen[i*width+j+k*size]=cercle_color;
        else if ( cont>0 ){
          cont=0; cont2=0;
          for(n=(float)i;n<=(i+1);n+=paso){
            for(m=(float)j;m<=(j+1);m+=paso){
            cont++;
            if( ((n-y0)*(n-y0)+(m-x0)*(m-x0))<radius*radius ) cont2++;
            }
          }
          m=(float) cont2/cont;
          imagen[i*width+j+k*size]=(unsigned char) (cercle_color*m+imagen[i*width+j+k*size]*(1-m));
        }
      }
    }
  }
}
}
#endif
