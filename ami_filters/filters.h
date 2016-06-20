/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


 #ifndef AMI_DLL_H
  #define AMI_DLL_H
#endif

/**
 * \file filters.h
 * \brief filters class AMI_DLL_H definition
 * \author Luis Alvarez \n \n
*/


#ifndef filters_H
#define filters_H

#include <math.h>
#include <vector>
#include "../ami_image/image.h"
#include "../ami_utilities/utilities.h"
#include "../ami_primitives/subpixel_image_contours.h"
using namespace std;

/**
 * \fn template <class  T,class U> void ami::filters<T,U>::
        gauss_conv(ami::image<T> &img,ami::image<U> &img_conv,double sigma_x,
                   double sigma_y,double precision)
 * \brief Gaussian convolution filter (using Alvarez-Mazorra implementation)
 * \author Luis Alvarez
*/
template <class  T,class U>
void gauss_conv(
ami::image<T> &img /** INPUT IMAGE (IT CAN BE = OUTPUT IMAGE */,
ami::image<U> &img_conv /** OUTPUT IMAGE */ ,
double sigma_x /** STANDARD DEVIATION IN THE x VARIABLE */,
double sigma_y /** STANDARD DEVIATION IN THE y VARIABLE */,
double precision /** PRECISION TO COMPUTE THE GAUSSIAN CONVOLUTION BIGGER
IS THE PRECISION, MORE ITERATIONS TO COMPUTE THE CONVOLUTION */)
{

  if(img.get_roi().size()==6)
  {
    img.get_roi_image(img_conv);
  }
  else
  {
    if(img.width()!=img_conv.width() ||  img.height()!=img_conv.height() ||
       img.nChannels()!=img_conv.nChannels())
    {
      img_conv=ami::image<U>(img.width(),img.height(),img.nChannels());
    }
    int k, k_end=img.width()*img.height()*img.nChannels();
		#ifdef _OPENMP
    #pragma omp parallel for shared(k_end) private(k)
		#endif
    for(k=0; k<k_end; k++)
      img_conv[k] = (U) img[k];
  }

  unsigned int width=img_conv.width();
  unsigned int height=img_conv.height();
  unsigned int nChannels=img_conv.nChannels();
  unsigned int image_size=width*height;
  unsigned int image_size_total=image_size*nChannels;

  if(precision < 0)
    precision = 0;

  unsigned int Nc_x=(unsigned int) (precision*sigma_x>1 ? precision*sigma_x : 1);
  unsigned int Nc_y=(unsigned int) (precision*sigma_y>1 ? precision*sigma_y : 1);
  unsigned int Nc_a=(Nc_x>Nc_y) ? Nc_x : Nc_y;
  Nc_x = Nc_y = Nc_a;
  double t_y = sigma_x*sigma_x/(2*Nc_x);
  double t_x = sigma_y*sigma_y/(2*Nc_y);
  double l_x, v_x, l_y, v_y, l_x_1, l_y_1;
	l_x = v_x = l_y = v_y = l_x_1 = l_y_1 = 0.;

  // WE CHECK IS STANDARD DEVIATIONS ARE 0
  if(t_x==0 && t_y==0) 
    return;

  // PARAMETERS X AND Y
  if(t_x > 0)
  {
     l_x=(float)(1.+2.*t_x-sqrt((double) 4*t_x+1))/(2*t_x);
     v_x=l_x/t_x;
     l_x_1=(1-l_x);
  }
  if(t_y > 0)
  {
     l_y=(float)(1.+2.*t_y-sqrt((double) 4*t_y+1))/(2*t_y);
     v_y=l_y/t_y;
     l_y_1=(1-l_y);
  }

  //MAIN LOOP
  int cont;
  int c, y, x, x_end, m, y_end;
  for(cont=0; cont<(int)Nc_a; cont++)
  {
    if(t_y > 0)
    {
			#ifdef _OPENMP
      #pragma omp parallel \
      shared(width,image_size_total,image_size,l_y,v_y,l_y_1) \
      private(c,x,y,m,x_end)
			#endif
      for(c=0; c<(int)image_size_total; c+=image_size)
      {
        #ifdef _OPENMP
				#pragma omp for nowait
				#endif
        for(y=0; y<(int)image_size; y+=width)
        {
          m=c+y;
          vector <U> paso(width);
          paso[0]= img_conv[m]/l_y_1;
          x_end=m+width;
          unsigned int n=1;
          for(x=m+1; x<x_end; x++,n++)
          {
            paso[n]= img_conv[x]+l_y*paso[n-1];
          }
          img_conv[x_end-1]=paso[width-1]/l_y_1;
          n=width-2;
          for(x=x_end-2; x>m; x--,n--)
          {
            img_conv[x]= paso[n]+l_y*img_conv[x+1];
          }
          img_conv[x]= paso[n]+l_y*img_conv[x+1];
          for(x=m; x<x_end; x++)
          {
           img_conv[x]*=v_y;
          }
        }
      }
    }
    if(t_y > 0)
    {
		 #ifdef _OPENMP
     #pragma omp parallel \
     shared(width,height,image_size_total,image_size,l_x,v_x,l_x_1) \
     private(c,x,y,m,y_end)
		 #endif
      for(c=0; c<(int)image_size_total; c+=image_size)
      {
        #ifdef _OPENMP
				#pragma omp for nowait
				#endif
        for(x=0; x<(int)width; x++)
        {
          m=c+x;
          vector <U> paso(height);
          paso[0]= img_conv[m]/l_x_1;
          y_end=m+image_size;
          unsigned int n=1;
          for(y=m+width; y<y_end; y+=width,n++)
          {
            paso[n]= img_conv[y]+l_x*paso[n-1];
          }
          img_conv[y_end-width]=paso[height-1]/l_x_1;
          n=height-2;
          for(y=y_end-2*width; y>m; y-=width,n--)
          {
            img_conv[y]= paso[n]+l_x*img_conv[y+width];
          }
          img_conv[y]= paso[n]+l_x*img_conv[y+width];
          for(y=m; y<y_end; y+=width)
          {
           img_conv[y]*=v_x;
          }
        }
      }
    }
  }
}


/**
 * \fn template <class  T,class U> void ami::filters<T,U>::
        grad(const ami::image<T> &img,ami::image<U> &grad_x,
             ami::image<U> &grad_y, const bool NeigborhoodType)
 * \brief Gradient Computation
 * \author Luis Alvarez
*/
template <class  T,class U>
void grad(
  const ami::image<T> &img /** INPUT IMAGE */,
  ami::image<U> &grad_x /** OUTPUT x-GRADIENT IMAGE */,
  ami::image<U> &grad_y /** OUTPUT y-GRADIENT IMAGE */,
  const bool NeigborhoodType /** =0 means 5 size neighborhood,
                                 =9  means 9 size neighborhood*/)
{
  if(img.get_roi().size() == 6)
  {
    ami::image<T> img2;
    img.get_roi_image(img2);
    grad(img2,grad_x,grad_y,NeigborhoodType);
    return;
  }


  unsigned int width=img.width();
  unsigned int height=img.height();
  unsigned int nChannels=img.nChannels();
  unsigned int size_image=width*height;
  unsigned int size_image_width=size_image-width;
  unsigned int total_image_size=size_image*nChannels;
  unsigned int width_1=width-1;

  if((int)width!=grad_x.width() ||  (int)height!=grad_x.height() ||
     (int)nChannels!=grad_x.nChannels())
  {
    grad_x=ami::image<U>(width,height,nChannels);
  }

  if((int)width!=grad_y.width() ||  (int)height!=grad_y.height() ||
     (int)nChannels!=grad_y.nChannels())
  {
    grad_y=ami::image<U>(width,height,nChannels);
  }

  if(NeigborhoodType == 0)
  {
    //ami::utilities u;
    vector < vector <unsigned int> > b=boundary_neighborhood_5n(width,height);
    int m, c, y, k, k_end, b_size=b.size();

		#ifdef _OPENMP
    #pragma omp parallel \
    shared(width,width_1,total_image_size,size_image,size_image_width,b,b_size)\
    private(c,y,m,k,k_end)
		#endif
    for(c=0; c<(int)total_image_size; c+=size_image)
    {
      #ifdef _OPENMP
			#pragma omp for nowait
			#endif
      for(y=width; y<(int)size_image_width; y+=width)
      {
        m=c+y;
        k_end=m+width_1;
        for(k=m+1; k<k_end; k++){
          grad_y[k]=(U) img[k+width]- (U) img[k-width];
          grad_x[k]=(U) img[k-1]- (U) img[k+1];
        }
      }
			#ifdef _OPENMP
      #pragma omp for nowait
			#endif
      for(int k=0; k<b_size; k++)
      { // IMAGE BOUNDARY GRADIENT ESTIMATION
        grad_y[b[k][0]+c]=(U) img[b[k][2]+c]- (U) img[b[k][1]+c];
        grad_x[b[k][0]+c]=(U) img[b[k][4]+c]- (U) img[b[k][3]+c];
      }
    }
  }
  else
  {
    vector < vector <unsigned int> > b=boundary_neighborhood_9n(width,height);
    double coef1,coef2,c1,d1;
    coef1=sqrt((double) 2.);
    coef2=0.25*(2.-coef1);
    coef1=0.5*(coef1-1);
    int m, c, y, k, l, k_end, b_size=b.size();
    #ifdef _OPENMP
		#pragma omp parallel \
    shared(width,width_1,total_image_size,size_image,size_image_width,b,b_size,\
           coef1,coef2) private(c,y,m,k,k_end,c1,d1)
		#endif
    for(c=0; c<(int)total_image_size; c+=size_image)
    {
      #ifdef _OPENMP
			#pragma omp for nowait
			#endif
      for(y=width; y<(int)size_image_width; y+=width)
      {
        m=c+y;
        k_end=m+width_1;
        for(k=m+1; k<k_end; k++)
        {
          c1=img[k+width+1]-img[k-width-1];
          d1=img[k-width+1]-img[k+width-1];
          grad_y[k]=(U)(coef1*((U) img[k+width]-(U)img[k-width])+coef2*(c1-d1));
          grad_x[k]=(U)(-(coef1*((U) img[k+1]- (U) img[k-1])+coef2*(c1+d1)));
        }
      }
			#ifdef _OPENMP
      #pragma omp for nowait
			#endif
      for(l=0; l<b_size; l++)
      { // IMAGE BOUNDARY GRADIENT ESTIMATION
        c1=img[c+b[l][6]]-img[c+b[l][7]];
        d1=img[c+b[l][5]]-img[c+b[l][8]];
        grad_y[c+b[l][0]]=(U)(coef1*((U) img[c+b[l][2]]- (U) img[c+b[l][1]])+
                              coef2*(c1-d1));
        grad_x[c+b[l][0]]=(U)(-(coef1*((U) img[c+b[l][3]]- (U) img[c+b[l][4]])+
                                coef2*(c1+d1)));
      }
    }
  }
}

/**
 * \fn float ami_median_float(int k, int n, float *x)
 * \brief FUNCTION TO COMPUTE THE MEDIAN OF A VECTOR IN FLOAT PRECISION
 *  \param[in]  k : the estimated median separates the k lower values of the data sample.
 *  \param[in]  n : vector size
 *  \param[in]  x : pointer to the vector
 *  \return median of x separating the k lower values.
 * \author Luis Alvarez
 */
float ami_median_float(int k, int n, float *x)
{
	int i, ir, j, l, mid;
	float a, *y, paso;
	ami_malloc1d(y,float,n);
	for(int mm=0; mm<n; mm++)
    y[mm]=x[mm];

	l=0;
	ir=n-1;
	for (;;)
  {
		if (ir <= l+1)
    {
			if (ir == l+1 && y[ir] < y[l])
      {
				paso=y[l];
        y[l]=y[ir];
        y[ir]=paso;
			}
            a=y[k];
			free(y);
			return a;
		} 
    else 
    {
			mid=(l+ir) >> 1;
			paso=y[mid]; y[mid]=y[l+1]; y[l+1]=paso;
			if (y[l] > y[ir])
      {
				paso=y[l]; y[l]=y[ir]; y[ir]=paso;
			}
			if (y[l+1] > y[ir])
      {
				paso=y[l+1]; y[l+1]=y[ir]; y[ir]=paso;
			}
			if (y[l] > y[l+1])
      {
				paso=y[l]; y[l]=y[l+1]; y[l+1]=paso;
			}
			i=l+1;
			j=ir;
			a=y[l+1];
			for (;;) 
      {
				do i++;
        while (y[i] < a);
        
				do j--;
        while (y[j] > a);
        
				if (j < i)
          break;
        
				paso=y[i]; y[i]=y[j]; y[j]=paso;
			}
			y[l+1]=y[j];
			y[j]=a;
			if (j >= k)
        ir=j-1;
			if (j <= k)
        l=i;
		}
	}

}

//THIS PROCEDURE FILLS THE MAP RECURSIVELY (USED BY CANNY)
static void fill_imap(ami::image<int> &imap,
                     float *gradien_norm,
                     double low_threshold,
                     int k)
{
    if(imap[k]==1) return;
    int width = imap.width();
    imap[k] = 3;
    //NEIGHBOURS = O, E, N, S, NO, NE, SO, SE
    int nv[8] = {k-1, k+1, k-width, k+width,
                 k-width-1, k-width+1, k+width-1, k+width+1};
    //EXPLORE THE NEIGHBOURS OF THE PIXEL
    for(int iv=0; iv<8; iv++)
    {
        int nk = nv[iv];
        if(imap[nk]==3) continue;
        if((imap[nk]==0) && (gradien_norm[nk]>low_threshold))
            fill_imap(imap, gradien_norm, low_threshold, nk);
        else
            imap[nk] = 1;
    }
}

/**
 * \fn template<class T> void ami::filters::
        canny(ami::image<T> input, ami::image<T> &output, vector<double> seno,
            vector<double> coseno, double low_threshold, double high_threshold)
 * \brief Computes the edges with the Canny method
 * \author Luis Alvarez and Daniel Santana-Cedrés
*/
template<class T>
void canny(ami::image<T> input /**INPUT IMAGE (GRAY SCALE)*/,
			 ami::image<T> &output /**OUTPUT IMAGE WITH THE EDGES*/,
			 float *seno /**SINUS OF THE ORIENTATION*/,
			 float *coseno /**COSINUS OF THE ORIENTATION*/,
			 int *x /**COORDINATE X OF THE POSITION*/,
			 int *y /**COORDINATE Y OF THE POSITION*/,
			 const float per_low /**PERCENTAGE FOR THE LOW THRESHOLD (BETWEEN 0 AND 1)*/,
			 const float per_high /**PERCENTAGE FOR THE HIGH THRESHOLD (BETWEEN 0 AND 1)*/)
{
	//cout << "AMI::FILTERS::CANNY STARTS" << endl;
	bool neighborhoodtype = 9; //9 NEIGHBOURS FOR THE GRADIENT
	int width  = input.width();
	int height = input.height();
	int size   = width*height;
	//THRESHOLDS
	float low_threshold,high_threshold;
	//MAP
	ami::image<int> imap(width,height,1,0);
	//NORM OF THE GRADIENT
	float *gradient_norm = new float[size];
	//GRADIENT COMPONENTS
	ami::image<float> x_grad(width,height,input.nChannels(),0);
	ami::image<float> y_grad(width,height,input.nChannels(),0);
	//IMAGE FOR THE RESULT OF THE CONVOLUTION WITH THE GAUSSIAN
	ami::image<float> blurred(width,height,input.nChannels(),0);

	//WE APPLY A GAUSSIAN CONVOLUTION TO THE INPUT IMAGE (NOISE REDUCTION)
	gauss_conv(input,blurred,2.,2.,2.);

	//WE COMPUTE THE GRADIENT OF THE BLURRED IMAGE
	grad(blurred,x_grad,y_grad,neighborhoodtype);
	for(int k=0; k<size; k++)
			gradient_norm[k] = sqrt((x_grad[k]*x_grad[k]) + (y_grad[k]*y_grad[k]));

	//WE COMPUTE THE THRESHOLDS
	low_threshold  = ami_median_float(per_low*size, size, gradient_norm);
	high_threshold = ami_median_float(per_high*size, size, gradient_norm);
	//WE IMPOSE A MINIMUM VALUE FOR THE THRESHOLDS
	low_threshold=(low_threshold<2)?2:low_threshold;
	high_threshold=(high_threshold<4)?4:high_threshold;


	for(int i=0; i<height; i++)
	{
		for(int j=0; j<width; j++)
		{
			int pos     = i*width+j;
			//WE PUT THE MAP = 2 FOR VALUES > HIGH_THRESHOLD
			if(gradient_norm[pos] > high_threshold)
			{
				imap[pos] = 2;
			}
			//LIMITS OF THE IMAGE
			if((i==0) || (i==height-1) || (j==0) || (j==width-1))
			{
				imap[pos] = 1;
			}
		}
	}

	//NON-MAXIMUM SUPPRESSION
	float c = 1/(sqrt(2.));
	for(int i=1; i<height-1; i++)
	{
		for(int j=1; j<width-1; j++)
		{
			int k    = i*width+j;
			//ORIENTATIONS
			float o1 = fabs(x_grad[k]);
			float o2 = fabs(y_grad[k]);
			float o3 = fabs(c*(x_grad[k]+y_grad[k]));
			float o4 = fabs(c*(x_grad[k]-y_grad[k]));
			//NEIGHBOURS
			float neigh1=0;
			float neigh2=0;
			//WE SELECT THE MAXIMUM OF THE ORIENTATIONS
			float candidate = max(o1,max(o2,max(o3,o4)));
			//WE TAKE THE VALUE OF THE NEIGHBOURS BASED ON MAXIMUM ORIENTATION
			//0º
			if(candidate == o1)
			{
				//O - E
				neigh1 = gradient_norm[k-1];
				neigh2 = gradient_norm[k+1];
			}
			//90º
			if(candidate == o2)
			{
				//N - S
				neigh1 = gradient_norm[k-width];
				neigh2 = gradient_norm[k+width];
			}
			//45º
			if(candidate == o3)
			{
				//NE - SO
				neigh1 = gradient_norm[k-width+1];
				neigh2 = gradient_norm[k+width-1];
			}
			//135º
			if(candidate == o4)
			{
				//NO - SE
				neigh1 = gradient_norm[k-width-1];
				neigh2 = gradient_norm[k+width+1];
			}
			//WE UPDATE THE RESULT IN FUNCTION OF THE MAXIMUM VALUE OF
			//THE GRADIENT MODULE BETWEEN THE POINT AND ITS NEIGHBOURS
			float p_grad = gradient_norm[k]; //GRADIENT VALUE IN THE POINT
			if((p_grad>neigh1) && (p_grad>neigh2))
			{
				if(gradient_norm[k] >= high_threshold)
					imap[k] = 2;
				else
					imap[k] = 0;
			}
			else
			{
				imap[k] = 1;
			}
		}
	}

	//HISTERESIS PROCESS OF THE CANNY METHOD
	for(int k=0; k<size; k++)
	{
		if(imap[k]==2)
		{
			//NEIGHBOURS = O, E, N, S, NO, NE, SO, SE
			int nv[8] = {k-1, k+1, k-width, k+width,
									 k-width-1, k-width+1, k+width-1, k+width+1};
			for(int iv=0; iv<8; iv++)
			{
				int nk = nv[iv];
				if(imap[nk]>1) continue;
				if((imap[nk]==0) && (gradient_norm[nk]>low_threshold))
				{
					fill_imap(imap, gradient_norm, low_threshold, nk);
				}
				else
					imap[nk] = 1;
			}
		}
	}

	/*float maxi=0.;
	for(int i=0; i<size; i++)
	{
			maxi = max(maxi, gradient_norm[i]);
	}*/

	//FILL THE OUTPUT IMAGE TAKING INTO ACCOUNT THE IMAP IMAGE VALUES
	//THE EDGES ARE THE PIXELS WITH MAP > 1
	for(int i=1; i<height-1; i++)
	{
		for(int j=1; j<width-1; j++)
		{
			int pos = i*width+j;
			if(imap[pos]>1)
			{
				output[pos] = 255;
				//VALUE OF THE EDGES (NORMALIZED TO THE RANGE OF 0..255)
				/*output[pos]        = (unsigned char)(255*
																						 gradient_norm[pos]/maxi);
				output[pos+size]   = (unsigned char)(255*
																						 gradient_norm[pos]/maxi);
				output[pos+size*2] = (unsigned char)(255*
																						 gradient_norm[pos]/maxi);*/

				//ORIENTATION OF THE EDGE
				coseno[pos] = x_grad[pos]/gradient_norm[pos];
				seno[pos]   = y_grad[pos]/gradient_norm[pos];

				//POSITION OF THE EDGE
				x[pos]      = j;
				y[pos]      = i;
			}
			else
			{
				//BACKGROUND
				output[pos]        = 0;
				/*output[pos+size]   = 0;
				output[pos+size*2] = 0;*/
			}
		}
	}
	delete []gradient_norm;
	//cout << "AMI::FILTERS::CANNY ENDS" << endl;
}

/**
 * \fn template<class T> ami::subpixel_image_contours ami::filters::
        canny(ami::image<T> input, ami::image<T> &edges, float canny_low_threshold, float canny_high_threshold)
 * \brief Computes the edges with the Canny method and return a subpixel_image_contours object
 * \author Luis Alvarez and Daniel Santana-Cedrés
*/
template<class T>
ami::subpixel_image_contours canny(
ami::image<T> input /**INPUT IMAGE (GRAY SCALE)*/,
ami::image<T> &edges /**OUTPUT IMAGE WITH THE EDGES*/,
const float canny_low_threshold /**PERCENTAGE FOR THE LOW THRESHOLD (BETWEEN 0 AND 1)*/,
const float canny_high_threshold /**PERCENTAGE FOR THE HIGH THRESHOLD (BETWEEN 0 AND 1)*/)
{
  int size_=input.width()*input.height();
  float *seno   = new float[size_];//edge point orientation
  float *coseno = new float[size_];//edge point orientation
  int *x_pos    = new int[size_];//edge point location
  int *y_pos    = new int[size_];//edge point location
	vector<int> index(size_);
	int m=0;
  // WE COMPUTE THE EDGES
  canny(input,edges,seno,coseno,x_pos,y_pos,canny_low_threshold,canny_high_threshold);

  //Filling subpixel_image_contours object to call Hough transform
  ami::subpixel_image_contours contours(input.width(),input.height());
  for(int i=0; i<size_; i++)
  {
    if(edges[i]>0)
    { //edge point condition
      contours.get_c()[i]      = true;
      contours.get_x()[i]      = (float)x_pos[i];
      contours.get_y()[i]      = (float)y_pos[i];
      contours.get_coseno()[i] = seno[i];
      contours.get_seno()[i]   = coseno[i];
			index[m++]=i;
    }
    else{
      contours.get_c()[i] = false;
    }
  }
	index.resize(m);
  contours.set_index(index);
  delete []seno;
  delete []coseno;
  delete []x_pos;
  delete []y_pos;

  invert(edges,edges);//inverting black-white colors

  return(contours);

}
#endif
