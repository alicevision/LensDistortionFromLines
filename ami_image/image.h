/*
 * Copyright (c) 2010-2011, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under, at your option, the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version, or
 * the terms of the simplified BSD license.
 *
 * You should have received a copy of these licenses along this
 * program. If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */


/**
 * \file image.h
 * \brief Class to store multichannel image
 * \author Luis Alvarez, Pedro Henríquez \n \n
*/
#ifndef image_H
#define image_H

#include <vector>
#include <iostream>
#include <istream>
#include <ostream>
#include <typeinfo>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "io_png/io_png.h"

namespace ami {

// ITERATOR ALIAS

typedef std::vector<unsigned char>::iterator usIt;
typedef std::vector<int>::iterator iIt;
typedef std::vector<short>::iterator sIt;
typedef std::vector<float>::iterator fIt;
typedef std::vector<double>::iterator dIt;

typedef std::vector<unsigned char>::const_iterator uscIt;
typedef std::vector<int>::const_iterator icIt;
typedef std::vector<short>::const_iterator scIt;
typedef std::vector<float>::const_iterator fcIt;
typedef std::vector<double>::const_iterator dcIt;

}


namespace ami {

////////////////////////////////////////////////////////////////////////////////
/**
 * \class image
 * \brief Class to store multiChannel images and basic methods
 * \author Luis Alvarez
 */
template <class T>
class image{
  std::vector <T> image_ /** std vector to allocate image "channel 1 + channel 2
                              + ..." */;
  int width_ /** image width */;
  int height_ /** image height */;
  int nChannels_ /** number of image channels */;
  std::vector <int> roi_ /** ROI : Region Of Interest =(x0,x1,y0,y1,c0,c1)
                            subwindow corners used (if defined) to operate only
                            in a subwindow area */;

  std::vector <float> origin_ /** world coordinate origin */;
  std::vector <float> pixel_size_ /** pixel size */;

  public:

  // CONSTRUCTOR - DESTRUCTOR METHODS
  image();
  image(const int width,const int height);
  image(const T *vector,const int width,const int height,const int channels);
  image(const int width,const int height, const int nChannels, const T &a);
  image(const int width,const int height, const T &a, const T &b, const T &c);
  image(const int width,const int height, const int nChannels);
  image(std::string name /** INPUT FILE NAME */ );


  image(const image<bool> &img_in, bool isbool);
  image(const image<T> &img_in);      //----------------------------------------
  ~image();

  // CLASS ELEMENT ACCESS METHODS
  inline std::vector <T> *get_image(){return &image_;}
  inline int width() const {return width_;}
  inline int height() const {return height_;}
  inline int nChannels() const {return nChannels_;}
  inline int size() const {return image_.size();}
  inline T& operator[](const int &i);
  inline T& operator()(const int &x,const int &y,const int &channel);
  inline const T& operator[](const int &i) const;

  // BASIC OPERATION METHODS
  int write(std::string name);
  int write_bool(std::string name);
  int read(std::string name);
  image & operator=(const image &image2);
  void clear();
  void init(const int width,const int height, const int nChannels);
  void set_nchannels(int nchannel1);
  void set_size(int width,int height);
  void imageMirrored();
  image<T> const resize(const int width, const int height) const;
  image<T> const resize_no_omp(const int width, const int height) const ;
  image<T> const linear_transform(const double H[3][3]) const ;
  image<T> const sampling() const ;

  // BASIC METHODS TO MANIPULATE SUBWINDOWS (ROI : REGION OF INTEREST )
  std::vector <int> get_roi() const{return roi_;}
  void roi_clear(){roi_.clear();}
  void set_roi(const int x0,const int x1,const int y0,const int y1);
  void set_roi(const int x0,const int x1,const int y0,const int y1,
               const int c0,const int c1);
  void set_roi(const std::vector<int> roi);
  template <class U> void  get_roi_image(image<U> &image2) const;
  image<T>  const get_roi_image(const std::vector <int> &roi2) const;
  template <class U> void  set_roi_image(const image<U> &image2);

  // BASIC METHODS TO MANIPULATE IMAGE WITH ITERATORS
  typedef typename std::vector <T>::iterator iterator;
  inline iterator begin() {return(image_.begin());}
  inline iterator end() {return(image_.end());}

  //CHANGE CHANNEL
  void rgb_to_hsv(image<unsigned char> &H,image<unsigned char> &S,
                  image<unsigned char> &V);
  template <class U>
  void hsv_to_rgb(image<U> &image);

  // ACCESS VECTOR BOOL
  T get_value(const int x);
  void set_value(const int x,T value);

  // HISTOGRAM
  template<class U> void get_histogram(std::vector<U> &histograma,int channel);


};

////////////////////////////////////////////////////////////////////////////////
/**
* \fn template<class T> image<T> const image<T>::sampling() const
* \brief method to downsample the image by a factor of 2. We filter the image
          before downsampling
*/
template<class T>
image<T> const image<T>::sampling() const
{
  // WE DEFINE SAMPLED IMAGE
  image<T> sampled_img;
  int width_s=width_/2;
  int height_s=height_/2;
  sampled_img=image<T>(width_s,height_s,nChannels_);
  // WE FILTER THE IMAGE AND DOWNSAMPLING IT
  for(int n=0;n<nChannels_;n++){
    int n0=n*width_*height_;
    int n_s=n*width_s*height_s;
    int i_end=height_s-1;
		#ifdef _OPENMP
    #pragma omp parallel for shared(i_end,n0,n_s,sampled_img,width_s)
		#endif
    for(int i=1;i<i_end;i++){
      int m0=n0+2*i*width_;
      int m_s=n_s+i*width_s;
      int j_end=width_s-1;
      for(int j=1;j<j_end;j++){
        int m=m0+2*j;
        sampled_img[m_s+j]=0.25*image_[m]+0.125*(image_[m+1]+image_[m-1]+
                                                 image_[m+width_]+
                                                 image_[m-width_])+
                                          0.0625*(image_[m+width_+1]+
                                                  image_[m+width_-1]+
                                                  image_[m-width_+1]+
                                                  image_[m-width_-1]);
      }
    }
  }

  // WE SAMPLE IMAGE BORDER VERTICAL LINES
  for(int n=0;n<nChannels_;n++){
    int n0=n*width_*height_;
    int n_s=n*width_s*height_s;
    for(int i=0;i<height_s;i++){
      int i0=2*i;
      int m0=n0+i0*width_;
      int m_s=n_s+i*width_s;
      for(int j=0;j<width_s;j+=width_s-1){
        int j0=2*j;
        int m=m0+j0;
        int mN=m,mS=m,mW=m,mE=m,mNE=m,mNW=m,mSE=m,mSW=m;
        if(i0>0){ mS-=width_; mSE-=width_; mSW=m-width_; }
        if(i0<(height_-1)){ mN+=width_; mNE+=width_; mNW+=width_; }
        if(j0>0) {mW-=1; mNW-=1; mSW-=1;}
        if(j0<(width_-1)){ mE+=1; mNE+=1; mSE+=1;}

        sampled_img[m_s+j]=0.25*image_[m]+0.125*(image_[mN]+image_[mS]+
                                                 image_[mW]+image_[mE])+
                                          0.0625*(image_[mNW]+image_[mNE]+
                                                  image_[mSW]+image_[mSE]);
      }
    }
  }

  // WE SAMPLE IMAGE BORDER HORIZONTAL LINES
  for(int n=0;n<nChannels_;n++){
    int n0=n*width_*height_;
    int n_s=n*width_s*height_s;
    for(int i=0;i<height_s;i+=height_s-1){
      int i0=2*i;
      int m0=n0+i0*width_;
      int m_s=n_s+i*width_s;
      for(int j=1,j_end=width_s-1;j<j_end;j++){
        int j0=2*j;
        int m=m0+j0;
        int mN=m,mS=m,mW=m,mE=m,mNE=m,mNW=m,mSE=m,mSW=m;
        if(i0>0){ mS-=width_; mSE-=width_; mSW=m-width_; }
        if(i0<(height_-1)){ mN+=width_; mNE+=width_; mNW+=width_; }
        if(j0>0) {mW-=1; mNW-=1; mSW-=1;}
        if(j0<(width_-1)){ mE+=1; mNE+=1; mSE+=1;}

        sampled_img[m_s+j]=0.25*image_[m]+0.125*(image_[mN]+image_[mS]+
                                                 image_[mW]+image_[mE])+
                                          0.0625*(image_[mNW]+image_[mNE]+
                                                  image_[mSW]+image_[mSE]);
      }
    }
  }
  return(sampled_img);
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> template <class U> void image<T>::
    get_histogram(std::vector<U> &histograma, int channel)
 * \brief Get histogram of a channel of image
 * \author Carlos Falcon
*/
template <class T> template <class U>
void image<T>::get_histogram(std::vector<U> &histograma/** histograma de
                                                           dimension 256*/,
                              int channel/**channel to obtain the histogram*/)
{
  if (nChannels_<channel)
    return;
  if (channel <= 0)
    return;
  int size=width_*height_;
  histograma.resize(256);
	#ifdef _OPENMP
  #pragma omp parallel for
	#endif
  for(int i=0;i<256;i++) histograma[i]=0;

  int i_end=size*channel;
  //#pragma omp parallel for default(none) shared(size, channel,i_end)
  for(int i=size*(channel-1);i<i_end;i++)
          histograma[image_[i]]++;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> void image<T>::init(const int width,const int height,
                                              const int nChannels)
 * \brief initializes the image taking memory
 * \author Carlos Falcon
*/
template <class T>
void image<T>::init(const int width,const int height, const int nChannels)
{
  width_=width;
  height_=height;
  nChannels_=nChannels;
  int size=width_*height_*nChannels_;
  image_.resize(size);
}


////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> void image<T>::set_value(const int x,T value)
 * \brief set value to x position of vector
 * \author Carlos Falcon
*/
template <class T>
void image<T>::set_value(const int x,T value)
{
  if (x<=width_*height_*nChannels_)
    image_[x]=value;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> T image<T>::get_value(const int x)
 * \brief get value to x position of vector
 * \author Carlos Falcon
*/
template <class T>
T image<T>::get_value(const int x)
{
  T value=0;
  if (x<=width_*height_*nChannels_)
    return image_[x];
  return value;
}


////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T>::image(const T *vector,const int width,
                                          const int height,const int channels)
 * \brief image from T vector
 * \author Carlos Falcon
*/

template <class T>
image<T>::image(const T *vector,const int width,const int height,
                const int channels)
{
  width_=width ;
  height_=height;
  nChannels_=channels;
  int size=width_*height_*nChannels_;
  image_.resize(size);
  size=width_*height_;
	#ifdef _OPENMP
  #pragma omp parallel \
    shared(size, vector)
	#endif
    for(int n=0;n<nChannels();n++){
      int n2_=n*size;
      #ifdef _OPENMP
			#pragma omp for nowait
			#endif
      for(int i=0;i<size;i++)
      {
        image_[i+n2_]=vector[i+n2_];
      }
    }

//  for (int i=0;i< size;i++)
//    image_[i]=vector[i];
}
////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T>::image(const image<T> &img_in)
 * \brief image copy constructor
 * \author Carlos Falcon
*/

template <class T>
image<T>::image(const image<T> &img_in)
{
  width_=img_in.width() ;
  height_=img_in.height();
  nChannels_=img_in.nChannels();
  int size=width_*height_*nChannels_;
  image_.resize(size);
  size=width_*height_;
	#ifdef _OPENMP
  #pragma omp parallel \
    shared(size)
	#endif
    for(int n=0;n<img_in.nChannels();n++){
      int n2_=n*size;
      #ifdef _OPENMP
			#pragma omp for nowait
			#endif
      for(int i=0;i<size;i++)
      {
        image_[i+n2_]=img_in[i+n2_];
      }
    }

    if(img_in.get_roi().size() == 6){
        roi_.resize(6);
        roi_.at(0)=img_in.get_roi().at(0);
        roi_.at(1)=img_in.get_roi().at(1);
        roi_.at(2)=img_in.get_roi().at(2);
        roi_.at(3)=img_in.get_roi().at(3);
        roi_.at(4)=img_in.get_roi().at(4);
        roi_.at(5)=img_in.get_roi().at(5);
    }
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> template<class U> void image<T>::
        hsv_to_rgb(image<U> &image)
 * \brief function to get colour model RGB from HSV (unsigned char version)
 * \author Carlos Falcon
*/
template <class T> template<class U>
void image<T>::hsv_to_rgb(image<U> &image)
{
  if (nChannels_!=3) return;
  int size=width_*height_;
  for (int i=0; i<size;i++)
  {
    U hue=  image_[i];
    U value= image_[i+size];
    U saturation= image_[i+2*size];
    U red=0,blue=0,green=0;

    if(hue<=43){
    red=image[i]=value;
    image[i+size*2]=value-saturation*value/255;
    blue=image[i+size]=hue*(red-blue)/43+blue;
    }
    if(hue>43 && hue<=85){
    green=image[i+size]=value;
    blue=image[i+size*2]=value-saturation*value/255;
    image[i]=-((hue-85)*(green-blue)/43-blue);
    }
    if(hue>85 && hue<=128){
    red=image[i+size]=value;
    green=image[i]=value-saturation*value/255;
    image[i+size*2]=(hue-85)*(green-red)/43+red;
    }
    if(hue>128 && hue<=171){
    blue=image[i+2*size]=value;
    red=image[i]=value-saturation*value/255;
    image[i+size]=((hue-171)*(blue-red)/43+red);
    }
    if(hue>171 && hue<=214){
    blue=image[i+2*size]=value;
    green=image[i+size]=value-saturation*value/255;
    image[i]=-(hue-171)*(blue-green)/43+green;
    }
    if(hue>214){
    red=image[i]=value;
    green=image[i+size]=value-saturation*value/255;
    image[i+2*size]=((hue-214)*(red-green)/43+green);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> void image<T>::rgb_to_hsv(image<unsigned char> &H,
                                                    image<unsigned char> &S,
                                                    image<unsigned char> &V)
 * \brief function to get colour model HSV from RGB (unsigned char version)
 * \author Carlos Falcon
*/
template <class T>
void image<T>::rgb_to_hsv(image<unsigned char> &H,image<unsigned char> &S,
                          image<unsigned char> &V)
{
  if (nChannels_!=3) return;
  if (H.width()!=width_ || S.width()!=width_ || V.width()!=width_)return;

  int size_=width_*height_;
  for (int i=0; i<size_;i++)
  {
    unsigned char red, green, blue;
    red= (unsigned char) image_[i];
    green=(unsigned char) image_[i+size_];
    blue=(unsigned char) image_[i+2*size_];
    unsigned char rgb_min,rgb_max;
    rgb_min=rgb_max=red;
    if(rgb_min>green) rgb_min=green;
    if(rgb_min>blue) rgb_min=blue;
    if(rgb_max<green) rgb_max=green;
    if(rgb_max<blue) rgb_max=blue;
    if(rgb_max == 0){
      H[i]=0;
      S[i]=0;
      V[i]=0;  //HSV
      continue;
    }
    V[i]=(unsigned char) rgb_max;   //V
    S[i]=(unsigned char) ((int) 255*(rgb_max - rgb_min)/rgb_max);  //S
    if(S[i]==0){
      H[i]=0; //H
      continue;
    }
      //H
    if(rgb_max==red){
      H[i]=(unsigned char) ((int) 0+43*((int) green - blue)/
                            (rgb_max - rgb_min));
    }
    else if (rgb_max == green) {
      H[i]= (unsigned char) ((int) 85+ 43*((int) blue-red)/(rgb_max - rgb_min));
    }
    else{
      H[i]= (unsigned char) ((int) 171 + 43*((int) red - green)/
                             (rgb_max - rgb_min));
    }
   }
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> void image<T>::set_roi(const int x0,const int x1,
                                                 const int y0,const int y1)
 * \brief function to set a subwindow border
 * \author Luis Alvarez
*/
template <class T>
void image<T>::set_roi(const int x0,const int x1,const int y0,const int y1)
{
  if(x0>x1 || y0>y1 || x0>width_ || y0>height_) return;

  roi_.resize(6);
  roi_.at(0)=x0;
  roi_.at(1)=x1>width_?width_:x1;
  roi_.at(2)=y0;
  roi_.at(3)=y1>height_?height_:y1;
  roi_.at(4)=0;
  roi_.at(5)=nChannels_;

}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn void image::set_roi(const int x0,const int x1,const int y0,
   const int y1,const int c0,const int c1)
 * \brief function to set a subwindow border
 * \author Luis Alvarez
*/
template <class T>
void image<T>::set_roi(const int x0,const int x1,const int y0,
                          const int y1,const int c0, const int c1)
{
  if(x0>x1 || y0>y1 || x0>width_ || y0>height_ || c0>nChannels_ || c0>c1 )
    return;

  roi_.resize(6);
  roi_.at(0)=x0;
  roi_.at(1)=x1>width_?width_:x1;
  roi_.at(2)=y0;
  roi_.at(3)=y1>height_?height_:y1;
  roi_.at(4)=c0;
  roi_.at(5)=c1;

}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> void image<T>::set_roi(const std::vector<int> roi)
 * \brief function to set a subwindow border
 * \author Luis Alvarez
*/
template <class T>
void image<T>::set_roi(const std::vector<int> roi){  roi_=roi; }

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> template <class U> void  image<T>::
        set_roi_image(const image<U> &image2)
 * \brief function to fill the image with a subwindow image
 * \author Luis Alvarez
*/

template <class T> template <class U> void  image<T>::
  set_roi_image(const image<U> &image2)
{

  if(roi_.size()<6) return;
  int x0=roi_.at(0);
  int x1=roi_.at(1);
  int y0=roi_.at(2);
  int y1=roi_.at(3);
  int c0=roi_.at(4);
  int c1=roi_.at(5);

  if( x1<=width_ &&  y1<=height_ && x0<=x1 && y0<=y1 && c1<=nChannels_ &&
     c0<=c1){
    int width=x1-x0;
    int height=y1-y0;
    int size=width*height;

    if(image2.width()!=width || image2.height()!=height ||
       image2.nChannels()!=nChannels_) return;

    int i,j,n,i2,i2_,n2,n2_,size_=width_*height_;

		#ifdef _OPENMP
    #pragma omp parallel \
    shared(x0,x1,y0,y1,c0,c1,width,size,size_) \
    private(n,i,j,n2,n2_,i2,i2_)
		#endif
    for(n=c0;n<c1;n++){
      n2=(n-c0)*size;
      n2_=n*size_;
			#ifdef _OPENMP
      #pragma omp for nowait
			#endif
      for(i=y0;i<y1;++i){
        i2_=n2_+i*width_;
        i2=n2+(i-y0)*width-x0;
        for(j=x0;j<x1;++j){
          image_[i2_+j]=(T) image2[i2+j] ;
        }
      }
    }
    return;
  }
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T>  const image<T>::
        get_roi_image(const std::vector <int> &roi2) const
 * \brief function to get a subwindow image
 * \author Luis Alvarez
*/
template <class T>
image<T>  const image<T>::get_roi_image(const std::vector <int> &roi2) const
{
  int x0,x1,y0,y1,c0,c1;

  // IF THERE IS NOT ROI DEFINED WE RETURN AN EMPTY IMAGE
  if(roi2.size()<6){
     return(image<T>());
  }
  else{
    x0=roi2.at(0);
    x1=roi2.at(1);
    y0=roi2.at(2);
    y1=roi2.at(3);
    c0=roi2.at(4);
    c1=roi2.at(5);
  }

  int width,height;
  if( x1<= width_ &&  y1<= height_ && x0<x1 && y0<y1 &&  c1<= nChannels_ &&
     c0<c1){
    width=x1-x0;
    height=y1-y0;
  }
  else{
    return(image<T>());
  }

  image<T> image2(width,height,c1-c0);

  int size_=width_*height_;
  int size=width*height;
  int i,n,j,i2,i2_,n2,n2_;

	#ifdef _OPENMP
  #pragma omp parallel \
  shared(image2,x0,x1,y0,y1,c0,c1,width,size,size_) \
  private(n,i,j,n2,n2_,i2,i2_)
	#endif
  for(n=c0;n<c1;n++){
    n2=(n-c0)*size;
    n2_=n*size_;
    #ifdef _OPENMP
		#pragma omp for nowait
		#endif
    for(i=y0;i<y1;++i){
      i2_=n2_+i*width_;
      i2=n2+(i-y0)*width-x0;
      for(j=x0;j<x1;++j){
        image2[i2+j]= image_[i2_+j] ;
      }
    }
  }
  return(image2);
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> template <class U> void  image<T>::
        get_roi_image(image<U> &image2) const
 * \brief function to fill a subwindow image
 * \author Luis Alvarez
*/

template <class T> template <class U> void  image<T>::
  get_roi_image(image<U> &image2) const
{
  int x0,x1,y0,y1,c0,c1;

  // IF THERE IS NOT ROI DEFINED WE TAKE THE WHOLE IMAGE
  if(roi_.size()<6){
    x0=0;
    x1=width_;
    y0=0;
    y1=height_;
    c0=0;
    c1=nChannels_;
  }
  else{
    x0=roi_.at(0);
    x1=roi_.at(1);
    y0=roi_.at(2);
    y1=roi_.at(3);
    c0=roi_.at(4);
    c1=roi_.at(5);
  }

  if( x1<= width_ &&  y1<= height_ && x0<=x1 && y0<=y1 &&  c1<= nChannels_ &&
     c0<=c1){
    int width=x1-x0;
    int height=y1-y0;

    // WE ALLOCATE MEMORY IF NEEDED
    if(image2.width()!=width || image2.height()!=height ||
       image2.nChannels()!=nChannels_){
      image2=image<U>(width,height,c1-c0);
    }

    int size_=width_*height_;
    int size=width*height;
    int i,n,j,i2,i2_,n2,n2_;

		#ifdef _OPENMP
    #pragma omp parallel \
    shared(x0,x1,y0,y1,c0,c1,width,size,size_) \
    private(n,i,j,n2,n2_,i2,i2_)
		#endif
    for(n=c0;n<c1;n++){
      n2=(n-c0)*size;
      n2_=n*size_;
      #ifdef _OPENMP
			#pragma omp for nowait
			#endif
      for(i=y0;i<y1;++i){
        i2_=n2_+i*width_;
        i2=n2+(i-y0)*width-x0;
        for(j=x0;j<x1;++j){
          image2[i2+j]=(U) image_[i2_+j] ;
        }
      }
    }
    return;
  }
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn image<T> & image<T>::operator=(const image &image2)
 * \brief operator = (equality of images of different size is not allowed )
 * \author Luis Alvarez
*/
template <class T>
image<T> & image<T>::operator=(const image &image2)
{
  if( width_==image2.width() && height_==image2.height() &&
     nChannels_==3 && image2.nChannels()==1){
     if(image_.size()!=(unsigned int) (3*width_*height_) ) image_.resize(3*width_*height_);
     int size_=width_*height_;
     int size_2=2*size_;
     for(int m=0;m<size_;m++){
        image_[m]         = image2[m];
		image_[m+size_]   = image2[m];
		image_[m+size_2]  = image2[m];
     }
     return *this;
  }


  if( width_!=image2.width() || height_!=image2.height() ||
     nChannels_!=image2.nChannels()){
    width_=image2.width();
    height_=image2.height();
    nChannels_=image2.nChannels();
    image_.resize(image2.size());
  }

  int k,k_end=image2.size();
	#ifdef _OPENMP
  #pragma omp parallel for shared(k_end) private(k)
	#endif
  for(k=0;k<k_end;++k)
    image_[k]=  image2[k];

  return *this;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn image<T> const image<T>::linear_transform(const double H[3][3]) const ;
 * \brief function to apply a linear transform given by an homography H to the
          image. H goes to the output image to the input image
 * \author Luis Alvarez
*/
template <class T>
image<T> const image<T>::linear_transform(const double H[3][3]) const
{
  image<T> image2(width_,height_,nChannels_,255);

  int width_1=width_-1;
  int height_1=height_-1;
  for(int n=0;n<nChannels_;n++){
    int n0=n*width_*height_;
		#ifdef _OPENMP
    #pragma omp parallel for shared(H,image2,n,n0,width_1,height_1)
		#endif
    for(int i=0;i<height_;i++){
      for(int j=0;j<width_;j++){
        //WE COMPUTE THE HOMOGRAPHY TRANSFORMATION
        double z0=H[2][0]*j+H[2][1]*i+H[2][2];
        double x0=H[0][0]*j+H[0][1]*i+H[0][2];
        double y0=H[1][0]*j+H[1][1]*i+H[1][2];
        x0/=z0;
        y0/=z0;
        if(x0<0 || x0>=width_ || y0<0 || y0>=height_) continue;
        int y0i=y0;
        double dy0=y0-y0i;
        double dy0_1=1-dy0;
        int m0=n0+y0i*width_;
        int m2=n0+i*width_;
        int x0i=x0;
        double dx0=x0-x0i;
        double dx0_1=1-dx0;
        if( x0i<width_1){
          if( y0i<height_1){
            image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+width_+x0i]+
              dy0*dx0*image_[m0+width_+x0i+1]+
              dy0_1*dx0*image_[m0+x0i+1]
            );
          }
          else{
              image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+x0i]+
              dy0*dx0*image_[m0+x0i+1]+
              dy0_1*dx0*image_[m0+x0i+1]
            );
          }
        }
        else{
          if( y0i<height_1){
            image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+width_+x0i]+
              dy0*dx0*image_[m0+width_+x0i]+
              dy0_1*dx0*image_[m0+x0i]
            );
          }
          else{
              image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+x0i]+
              dy0*dx0*image_[m0+x0i]+
              dy0_1*dx0*image_[m0+x0i]
            );
          }
        }
      }
    }
  }
  return image2;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T> const image<T>::
        resize(const int width, const int height) const
 * \brief function to resize an image
 * \author Luis Alvarez
*/
template <class T>
image<T> const image<T>::resize(
const int width /** new image width */,
const int height /** new image height */) const
{
  printf("width=%d width_=%d height=%d height_=%d \n",width,width_,height,
         height_);

  if(width==width_ && height==height_) return (*this);

  if(width<=width_/2 && height<=height_/2){
     return(((*this).sampling()).resize(width,height));
  }

  if(height<100) return(resize_no_omp(width,height));

  if(width==0 || height==0) return(image<T> ());

  image<T> image2(width,height,nChannels_);


  // WE PROCESS IN THE CASE WE DO NOT NEED INTERPOLATION
  if( width_>width && width_%width==0  && height_>height && height_%height==0 ){
    int scale_x=width_/width;
    int scale_y=height_/height;
    for(int n=0;n<nChannels_;n++){
      int n0=n*width_*height_;
      int n2=n*width*height;
			#ifdef _OPENMP
      #pragma omp parallel for shared(image2,n,n0,n2,scale_x,scale_y)
			#endif
      for(int i=0;i<height;i++){
        int m0=n0+scale_y*i*width_;
        int m2=n2+i*width;
        for(int j=0;j<width;j++){
          image2[m2+j]=image_[m0+scale_x*j];
        }
      }
    }
    return image2;
  }

  // WE PROCESS THE CASE WHERE WE NEED A INTERPOLATION PROCEDURE
  double scale_x=(double) width_/width;
  double scale_y=(double) height_/height;

  int width_1=width_-1;
  int height_1=height_-1;
  for(int n=0;n<nChannels_;n++){
    int n0=n*width_*height_;
    int n2=n*width*height;
		#ifdef _OPENMP
    #pragma omp parallel for shared(image2,n,n0,n2,scale_x,scale_y,width_1,\
                                    height_1)
		#endif
    for(int i=0;i<height;i++){
      double y0=scale_y*i;
      int y0i=y0;
      double dy0=y0-y0i;
      double dy0_1=1-dy0;
      int m0=n0+y0i*width_;
      int m2=n2+i*width;
      for(int j=0;j<width;j++){
        double x0=scale_x*j;
        int x0i=x0;
        double dx0=x0-x0i;
        double dx0_1=1-dx0;
        if( x0i<width_1){
          if( y0i<height_1){
            image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+width_+x0i]+
              dy0*dx0*image_[m0+width_+x0i+1]+
              dy0_1*dx0*image_[m0+x0i+1]
            );
          }
          else{
              image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+x0i]+
              dy0*dx0*image_[m0+x0i+1]+
              dy0_1*dx0*image_[m0+x0i+1]
            );
          }
        }
        else{
          if( y0i<height_1){
            image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+width_+x0i]+
              dy0*dx0*image_[m0+width_+x0i]+
              dy0_1*dx0*image_[m0+x0i]
            );
          }
          else{
              image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+x0i]+
              dy0*dx0*image_[m0+x0i]+
              dy0_1*dx0*image_[m0+x0i]
            );
          }
        }
      }
    }
  }

  return image2;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T> const image<T>::
        resize_no_omp(const int width, const int height) const
 * \brief function to resize an image without omp
 * \author Luis Alvarez
*/
template <class T>
image<T> const image<T>::resize_no_omp(const int width, const int height) const
{

  if(width==0 || height==0) return(image<T> ());

  if(width==width_ && height==height_) return (*this);

  if(width<=width_/2 && height<=height_/2)
    return(((*this).sampling()).resize_no_omp(width,height));

  image<T> image2(width,height,nChannels_);


  // WE PROCESS IN THE CASE WE DO NOT NEED INTERPOLATION
  if( width_>width && width_%width==0  && height_>height && height_%height==0 ){
    int scale_x=width_/width;
    int scale_y=height_/height;
    for(int n=0;n<nChannels_;n++){
      int n0=n*width_*height_;
      int n2=n*width*height;
      for(int i=0;i<height;i++){
        int m0=n0+scale_y*i*width_;
        int m2=n2+i*width;
        for(int j=0;j<width;j++){
          image2[m2+j]=image_[m0+scale_x*j];
        }
      }
    }
    return image2;
  }

  // WE PROCESS THE CASE WHERE WE NEED A INTERPOLATION PROCEDURE
  double scale_x=(double) width_/width;
  double scale_y=(double) height_/height;

  int width_1=width_-1;
  int height_1=height_-1;
  for(int n=0;n<nChannels_;n++){
    int n0=n*width_*height_;
    int n2=n*width*height;
    for(int i=0;i<height;i++){
      double y0=scale_y*i;
      int y0i=y0;
      double dy0=y0-y0i;
      double dy0_1=1-dy0;
      if(y0i>=height_) y0i=height_-1;
      int m0=n0+y0i*width_;
      int m2=n2+i*width;
      for(int j=0;j<width;j++){
        double x0=scale_x*j;
        int x0i=x0;
        double dx0=x0-x0i;
        double dx0_1=1-dx0;
        if( x0i<width_1){
          if( y0i<height_1){
            image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+width_+x0i]+
              dy0*dx0*image_[m0+width_+x0i+1]+
              dy0_1*dx0*image_[m0+x0i+1]
            );
          }
          else{
              image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+x0i]+
              dy0*dx0*image_[m0+x0i+1]+
              dy0_1*dx0*image_[m0+x0i+1]
            );
          }
        }
        else{
          if( y0i<height_1){
            image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+width_+x0i]+
              dy0*dx0*image_[m0+width_+x0i]+
              dy0_1*dx0*image_[m0+x0i]
            );
          }
          else{
              image2[m2+j]=(T) (
              dy0_1*dx0_1*image_[m0+x0i]+
              dy0*dx0_1*image_[m0+x0i]+
              dy0*dx0*image_[m0+x0i]+
              dy0_1*dx0*image_[m0+x0i]
            );
          }
        }
      }
    }
  }
  return image2;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> int image<T>::write_bool(std::string name )
 * \brief function to write to disk a 8/16 bit tif image from boolean image
 * \author Carlos Falcon
*/
template <class T>
int image<T>::write_bool(std::string name /** image file name */)
{
  ami::image<unsigned char> aux(width_,height_,nChannels_,0);
  int size=width_*height_;
	#ifdef _OPENMP
  #pragma omp parallel for shared(size,aux)
	#endif
  for (int i=0; i<size;i++)
  {
    if (get_value(i))
      aux[i]=255;
  }
  return aux.write(name);
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> int image<T>::write( std::string name )
 * \brief Function to write an image to disk.
 * \param [in] name Filename of the image.
 * \author Pedro Henriquez and Daniel Santana-Cedrés
*/
template <class T>
int image<T>::write(std::string name /** image file name */)
{
  int pos=name.find_last_of('.');
  int size_ = width_*height_;
  unsigned char *red = new unsigned char[size_],
                *green = new unsigned char[size_],
                *blue = new unsigned char[size_];

  if(pos == (int)std::string::npos) return -1;

  for(int i=0; i<size_; i++)
  {
    red[i]   = image_[i];
    green[i] = image_[i+size_];
    blue[i]  = image_[i+size_*2];
  }

  std::string extension=name.substr(pos+1);
  if( (extension == std::string("png")) || (extension == std::string("PNG")))
  {
    int output_value = ami_write_png(strdup(name.c_str()),red,green,blue,width_,
                                     height_);
    delete []red;
    delete []green;
    delete []blue;
    return output_value;
  }
  delete []red;
  delete []green;
  delete []blue;
  printf("WRITE::Unrecognized image format\n");
  return -1;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> T& image<T>::operator()(const int &x,const int &y,
                                                  const int &channel)
 * \brief operator () to acces image value
 * \author Javier Martin
*/
template <class T>
T& image<T>::operator()(const int &x,const int &y,const int &channel){
  #ifdef IMAGE_DEBUG
    if((x+y*width_+channel*width_*height_)>=(image_.size())){
      printf("image<T>: bounds error vector access\n");
      printf("image size()=%d index to accces the vector =%d\n",
                                        image_.size(),x+y*width_+
                                        channel*width_*height_);
      int j; scanf("%d",&j);
      exit(0);
    }
  #endif
  return image_[x+y*width_+channel*width_*height_];
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> const T& image<T>::operator[](const int &i) const
 * \brief operator [] to acces image value
 * \author Luis Alvarez
*/
template <class T>
const T& image<T>::operator[](const int &i) const
{
  #ifdef IMAGE_DEBUG
    if(i>=(int)(image_.size())){
      printf("image<T>: bounds error vector access\n");
      printf("image size()=%d index to accces the vector =%d\n",
             (int)image_.size(),i);
      int j; scanf("%d",&j);
      exit(0);
    }
  #endif
  return image_[i];
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> T& image<T>::operator[](const int &i)
 * \brief operator [] to acces image value
 * \author Luis Alvarez
*/
template <class T>
T& image<T>::operator[](const int &i)
{
  #ifdef IMAGE_DEBUG
    if(i>=(int)(image_.size())){
      printf("image<T>: bounds error vector access\n");
      printf("image size()=%d index to accces the vector =%d\n",
             (int)image_.size(),i);
      int j; scanf("%d",&j);
      exit(0);
    }
  #endif
  return  (image_.at(i));
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T>::image(int width,int height)
 * \brief image constructor taking memory
 * \author Luis Alvarez
*/
template <class T>
image<T>::image(int width,int height)
{
  width_=width ;
  height_=height;
  nChannels_=1;
  image_.resize(width*height);
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn image<T>::image(const int width,const int height, const T &a, const T &b,
                       const T &c)
 * \brief image constructor taking memory and initialiting the channels: R with
          a value, G with b value and B with c value
 * \author Carlos Falcon
*/
template <class T>
image<T>::image(const int width,const int height, const T &a, const T &b,
                const T &c)
{
  width_=width ;
  height_=height;
  nChannels_=3;
  int size=width_*height_;
  image_.resize(size*3);
  int i;
	#ifdef _OPENMP
  #pragma omp parallel for shared(size) private(i)
	#endif
  for(i=0;i<size;i++)
    image_.at(i)=a;
  #ifdef _OPENMP
	#pragma omp parallel for shared(size) private(i)
	#endif
  for(i=size;i<2*size;i++)
    image_.at(i)=b;
  #ifdef _OPENMP
	#pragma omp parallel for shared(size) private(i)
	#endif
  for(i=2*size;i<3*size;i++)
    image_.at(i)=c;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn image<T>::image(const int width,const int height, const int nChannels,
                       const T &a)
 * \brief image constructor taking memory and initialiting with a value
 * \author Luis Alvarez
*/
template <class T>
image<T>::image(const int width,const int height,const int nChannels,const T &a)
{
  width_=width ;
  height_=height;
  nChannels_=nChannels;
  int size=width_*height_*nChannels_;
  image_.resize(size);
  int i;
	#ifdef _OPENMP
  #pragma omp parallel for shared(size) private(i)
	#endif
  for(i=0;i<size;i++) image_.at(i)=a;
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class T> image<T>::image(const int width,const int height,
                                          const int nChannels)
 * \brief image constructor taking memory and initialiting with a value
 * \author Luis Alvarez
*/
template <class T>
image<T>::image(const int width,const int height, const int nChannels)
{
  width_=width;
  height_=height;
  nChannels_=nChannels;
  int size=width_*height_*nChannels_;
  image_.resize(size);
}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn image<T>::~image()
 * \brief basic constructor
 * \author Luis Alvarez
*/
template <class T>
image<T>::image() : width_(0), height_(0), nChannels_(0) {}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn image<T>::image()
 * \brief basic destructor
 * \author Luis Alvarez
*/
template <class T>
image<T>::~image(){}

////////////////////////////////////////////////////////////////////////////////
/**
 * \fn  template <class T> image<T>::image(std::string name)
 * \brief Constructor from an image file
 * \param [in] name Input file name.
 * \author Pedro Henriquez
*/
template <class T>
image<T>::image(std::string name /** INPUT FILE NAME */ )
{
  read(name);
}

////////////////////////////////////////////////////////////////////////////////
/**
* \fn void image<T>::clear()
* \brief Function to put the image to 0 and clear its vectors.
* \author Pedro Henríquez
*/
template <class T>
void image<T>::clear()
{
  width_=0;
  height_=0;
  nChannels_=0;
  roi_clear();
  image_.clear();
}

////////////////////////////////////////////////////////////////////////////////
/**
* \fn void image<T>::set_nchannels(int nchannel1)
* \brief Change number of channels
* \param [in] nchannel1 New number of channels.
* \author Pedro Henríquez
*/
template <class T>
void image<T>::set_nchannels(int nchannel1)
{
  nChannels_=nchannel1;
  image_.resize(nChannels_*width_*height_);
}

////////////////////////////////////////////////////////////////////////////////
/**
* \fn void image<T>::set_size(int width,int height)
* \brief Change the image dimensions
* \param [in] width,height New image dimensions.
* \author Pedro Henríquez
*/
template <class T>
void image<T>::set_size(int width,int height)
{
  width_=width;
  height_=height;
  image_.resize(nChannels_*width_*height_);
}

////////////////////////////////////////////////////////////////////////////////
/**
* \fn image<T>::imageMirrored()
* \brief mirrored the image
* \author Javier Martín Abasolo
*/
template <class T>
void image<T>::imageMirrored()
{
  //#pragma omp parallel for
  // all of channels
  for(int k=0; k<nChannels_;k++){
    // reading for columns
    for(int i=0; i< height_;i++){
      // reading for rows
      for(int j=0; j<(width_/2);j++){
        T aux = image_[(k*height_*width_)+i*width_+j];
        // swapping mirrored values in the row
        image_[(k*height_*width_)+i*width_+j] = image_[(k*height_*width_)+
                                                        i*width_+width_-j-1];
        image_[(k*height_*width_)+i*width_+width_-j-1] = aux;
      }
    }
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
/**
* \fn int image<T>::read (std::string name)
* \brief Read an image selecting the library depend on the image format, returns
          0 when it can't load the image.
* \author Pedro Henríquez and Daniel Santana-Cedrés
*/
template <class T>
int image<T>::read (std::string name)
{
  unsigned char *red, *green, *blue;
  int width, height;
  // READ THE IMAGE WITH IO_PNG
  int output_value = ami_read_png(strdup(name.c_str()),&red,&green,&blue,
                                  &width,&height);
  if(output_value!=0)
  {
    return output_value;
  }
  // INITIALIZE ATTRIBUTES AND COUNTERS
  width_    = width;
  height_   = height;
  nChannels_= 3;
  int size_ = width_ * height_;
  image_.resize(size_*nChannels_);
  // FILL THE IMAGE VECTOR WITH THE INFORMATION FROM IO_PNG
  for(int i=0; i<size_; i++)
  {
    image_[i]         = (T)red[i];
    image_[i+size_]   = (T)green[i];
    image_[i+size_*2] = (T)blue[i];
  }
  // FREE MEMORY
  free(red);
  free(green);
  free(blue);
  return 0;
}

}

#endif
