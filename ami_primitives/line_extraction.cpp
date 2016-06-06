/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#ifndef AMI_DLL_CPP
  #define AMI_DLL_CPP
#endif

#include "line_extraction.h"
#include "../ami_lens_distortion/lens_distortion_procedures.h"
#include "../ami_utilities/utilities.h"
#include <omp.h>

//THIS FUNCTION CORRECTS THE ORIENTATION OF THE EDGE USING THE DIVISION MODEL.
//PARAMETERS: THE EDGE POSITION, THE ORIENTATION (SIN,COS) AND THE LENS
//            DISTORTION MODEL
vector<float> orientation_update(point2d<double> p, float seno, float coseno,
                                  lens_distortion_model ldm)
{
	vector<float> corrected_orientation;
	double a,b,norma;
	//POINT PLUS THE ORIENTATION
	point2d<double> p_ori(p.x + coseno, p.y + seno);
	//WE APPLY THE MODEL TO THE ORIGINAL POINT
	point2d<double> p_prime(ldm.evaluation(p));
	//WE APPLY THE MODEL TO THE POINT PLUS THE ORIENTATION
	point2d<double> p_ori_prime(ldm.evaluation(p_ori));

	//WE COMPUTE THE NEW ORIENTATION AND THE NORM
	a = p_ori_prime.x - p_prime.x;
	b = p_ori_prime.y - p_prime.y;
	norma = sqrt(a*a + b*b);

	if((a*a + b*b) <= 0)
	{
		corrected_orientation.push_back(coseno);
		corrected_orientation.push_back(seno);
		return corrected_orientation;
	}
	else
	{
		corrected_orientation.push_back(a/norma);
		corrected_orientation.push_back(b/norma);
		return corrected_orientation;
	}
}

/**
 * \fn double line_extraction::
        line_equation_distortion_extraction_improved_hough_quotient(
                              const subpixel_image_contours &subpixel_contours,
                              image_primitives &image_primitive,
                              float distance_point_line_max, int nlineas=13,
                              float angle_resolution=0.25,
                              float distance_resolution=2.,
                              float initial_distortion_parameter=-0.1,
                              float final_distortion_parameter=0.1,
                              float distortion_parameter_resolution=0.01,
                              float angle_point_orientation_max_difference=15.,
                              bool lens_distortion_estimation=true)
 * \brief Computation of lines using an improved version of Hough which includes
          1 parameter lens distortion division model
 * \author Luis Alvarez
 */
double line_equation_distortion_extraction_improved_hough(
  const subpixel_image_contours &subpixel_contours /** INPUT CONTOUR SUBPIXEL INFORMATION */,
  image_primitives &image_primitive /** OUTPUT IMAGE PRIMITIVES WHERE LINES AND DISTORTION MODEL IS DEFINED */,
  float distance_point_line_max /** INPUT MAX DISTANCE ALLOWED BETWEEN POINTS AND ASSOCIATED LINES */,
  int nlineas /** INPUT NUMBER OF LINES TO RETURN */,
  float angle_resolution /** INPUT ANGLE RESOLUTION */,
  float distance_resolution  /** INPUT DISTANCE RESOLUTION */,
  float initial_distortion_parameter /** INPUT INITIAL VALUE OF NORMALIZED DISTORTION PARAMETER */,
  float final_distortion_parameter /** INPUT FINAL VALUE NORMALIZED DISTORTION PARAMETER */,
  float distortion_parameter_resolution /** INPUT DISTORTION PARAMETER DISCRETIZATION STEP */,
  float angle_point_orientation_max_difference /** MAX DIFFERENCE (IN DEGREES) OF THE POINT ORIENTATION ANGLE
    AND THE LINE ANGLE */,
  bool lens_distortion_estimation /** BOOL TO CONTROL IF WE ESTIMATE THE LENS DISTORTION MODEL */,
  lens_distortion_model ini_ldm /** INITIAL DISTORTION MODEL */)
{
  int i,l,k;
  int width_score,height_score,depth_score;
  int width=subpixel_contours.get_width();
  int height=subpixel_contours.get_height();
  double x2,xc=(double) width/2.;
  double y2,yc=(double) height/2.;
  double best_distortion_parameter=0;
  double *seno,*coseno;
  double ami_pi=acos((double) -1.);
  double max_norm=(double) width*width+height*height;
  int nlineas_plus=nlineas;
  if( lens_distortion_estimation==true) nlineas_plus=nlineas;
  point2d<double> dc(xc,yc);
  double dmi=0.;

  // WE CHECK DISCRETIZATION PARAMETER VALUES
  if(distance_resolution<=0 || angle_resolution<=0 ){
    cout << "Problems in function line_equation_distortion_extraction_improved_hough()" << endl;
    cout << "Discrezation parameters lower or equal to 0 " << endl;
    cout << "Introduce a number to continue" << endl;
    int val = scanf("%d",&i);
		if(val==EOF)
			cout << "Read error" << endl;
    return -1;
  }
  // WE DEFINE SCORE VOLUME SIZE
  width_score=(int) (1.2*(2.*sqrt((double) width*width+height*height)/distance_resolution+2));
  height_score=(int) (180./angle_resolution);
  depth_score=distortion_parameter_resolution>0?(int)(1.+(final_distortion_parameter-initial_distortion_parameter)/distortion_parameter_resolution):1;
  // WE DEFINE ORIENTATION DISCRETIZATION VECTOR
  seno   =(double*)malloc(sizeof(double)*height_score);
  coseno =(double*)malloc(sizeof(double)*height_score);
  double paso_angle=angle_resolution*ami_pi/180.;
  for(l=0;l<height_score;l++) {
    coseno[l]=cos((double) l*paso_angle);
    seno[l]=sin((double) l*paso_angle);
  }
  subpixel_image_contours sp_contours_modifed(width, height);
  
	// UPDATE THE SUBPIXEL_CONTOURS OBJECT WITH THE INITIAL LENS DISTORTION MODEL
  vector<int> index = subpixel_contours.get_index();
  if(ini_ldm.get_d().size()<1)
  {
    // THE MODIFIED MODEL HAS THE DEFAULT VALUE. WE ONLY COPY
// 		#ifdef _OPENMP
//     #pragma omp parallel for
// 		#endif
    for(int i=0; i<(int)index.size(); i++)
    {
      sp_contours_modifed.get_x()[index[i]]      = subpixel_contours.get_x()[index[i]];
      sp_contours_modifed.get_y()[index[i]]      = subpixel_contours.get_y()[index[i]];
      sp_contours_modifed.get_coseno()[index[i]] = subpixel_contours.get_coseno()[index[i]];
      sp_contours_modifed.get_seno()[index[i]]   = subpixel_contours.get_seno()[index[i]];
      if(sp_contours_modifed.get_seno()[index[i]]<0){
        sp_contours_modifed.get_coseno()[index[i]] *=-1;
        sp_contours_modifed.get_seno()[index[i]]*=-1;
      }
      sp_contours_modifed.get_c()[index[i]]      = subpixel_contours.get_c()[index[i]];
    }
    if(ini_ldm.get_type()==DIVISION)
    {
      dmi = dc.norm();
      dmi *= dmi;
    }
  }
  else
  {
    // ELSE, WE UPDATE THE SUBPIXEL_CONTOURS OBJECT USING THE PROVIDED MODEL
// 		#ifdef _OPENMP
//     #pragma omp parallel for
// 		#endif
    for(int i=0; i<(int)index.size(); i++)
    {
      // WE UPDATE THE CONTOUR COMPONENT OF SUBPIXEL_CONTOURS OBJECT
      sp_contours_modifed.get_c()[index[i]] = subpixel_contours.get_c()[index[i]];
      // CORRECT THE POSITION
      point2d<double> ori(subpixel_contours.get_x()[index[i]],
                          subpixel_contours.get_y()[index[i]]);
      point2d<double> corrected(ini_ldm.evaluation(ori));
      sp_contours_modifed.get_x()[index[i]] = corrected.x;
      sp_contours_modifed.get_y()[index[i]] = corrected.y;

      // CORRECT THE ORIENTATION
      vector<float> corrected_orientation;
      double a,b,norma;
      // Point plus the orientation
      point2d<double> p_ori(ori.x + subpixel_contours.get_coseno()[index[i]],
                            ori.y + subpixel_contours.get_seno()[index[i]]);
      // Model applied to the original point
      point2d<double> p_prime(ini_ldm.evaluation(ori));
      // Model applied to the point plus the orientation
      point2d<double> p_ori_prime(ini_ldm.evaluation(p_ori));

      // We compute the new orientation and the norm
      a = p_ori_prime.x - p_prime.x;
      b = p_ori_prime.y - p_prime.y;
      norma = sqrt(a*a + b*b);

      if((a*a + b*b) <= 0)
      {
        corrected_orientation.push_back(subpixel_contours.get_coseno()[index[i]]);
        corrected_orientation.push_back(subpixel_contours.get_seno()[index[i]]);
      }
      else
      {
        corrected_orientation.push_back(a/norma);
        corrected_orientation.push_back(b/norma);
      }
      // WE UPDATE THE VALUE OF THE SIN AND COS
      sp_contours_modifed.get_coseno()[index[i]] = corrected_orientation[0];
      sp_contours_modifed.get_seno()[index[i]]   = corrected_orientation[1];
      if(sp_contours_modifed.get_seno()[index[i]]<0){
        sp_contours_modifed.get_coseno()[index[i]] *=-1;
        sp_contours_modifed.get_seno()[index[i]]*=-1;
      }
    }
    // WE UPDATE THE VALUE OF THE DISTORTION CENTER
    xc = ini_ldm.get_distortion_center().x;
    yc = ini_ldm.get_distortion_center().y;
    dc = point2d<double>(xc,yc);
    if(ini_ldm.get_type()==DIVISION)
    {
      point2d<double> corner(0,0);
      dmi= (dc-corner).norm2();
      corner.y=height;
      double distance_corner2=(dc-corner).norm2();
      if(distance_corner2>dmi)
        dmi=distance_corner2;
      corner.x=width;
      distance_corner2=(dc-corner).norm2();
      if(distance_corner2>dmi)
        dmi=distance_corner2;
      corner.y=0;
      distance_corner2=(dc-corner).norm2();
      if(distance_corner2>dmi)
        dmi=distance_corner2;
    }
  }

  double votation_score, max_votation_score = -1;
  for(int k=0; k<depth_score; k++)
  {
    float **score_k;
    // WE TAKE THE MEMORY FOR THE SLICE AND INITIALIZE IT TO ZERO
    ami2_malloc2d(score_k,float,width_score,height_score);
// 		ami2_malloc2d(score_k,float,height_score,width_score);
    for(int hs=0; hs<height_score; hs++)
    {
      for(int ws=0; ws<width_score; ws++)
      {
        score_k[ws][hs]=0.;
      }
    }
    // THE ANGLE INTERVAL
    double kdistortion = 0.0;
    if(ini_ldm.get_type()==DIVISION)
    {
      double p = initial_distortion_parameter + (k*distortion_parameter_resolution);
      kdistortion = -p / (dmi + dmi*p);
    }
    else
    {
      kdistortion = (initial_distortion_parameter+
                              k*distortion_parameter_resolution)/max_norm;
    }


    // WE BUILD THE MODEL FOR THIS ITERATION
    lens_distortion_model it_ldm;
		it_ldm.set_distortion_center(point2d<double>(xc,yc));
		it_ldm.get_d().resize(2);
		it_ldm.get_d()[0] = 1.;
		it_ldm.get_d()[1] = kdistortion;
		it_ldm.set_type(ini_ldm.get_type());
		
// 		#ifdef _OPENMP
// 		#pragma omp parallel for
// 		#endif
		for(int ind=0; ind<(int)index.size(); ind++)
		{
			double distancia,x2,y2,x3,y3,norm,distortion_parameter;
			int q,d,angle;
			q = index[ind];
			// WE COMPUTE THE DISTORTION MODEL ACCORDING TO THE TYPE
			x3 = x2 = sp_contours_modifed.get_x()[q];
			y3 = y2 = sp_contours_modifed.get_y()[q];
			if(it_ldm.get_d()[1]!=0)
			{
				if(ini_ldm.get_type()==DIVISION)
				{
					point2d<double> res = it_ldm.evaluation(point2d<double>(x2,y2));
					x3 = res.x;
					y3 = res.y;
				}
				else
				{
					norm=((x2-xc)*(x2-xc)+(y2-yc)*(y2-yc))/max_norm;
					distortion_parameter = norm*(initial_distortion_parameter+
																			 k*distortion_parameter_resolution);
					x3 = x2+distortion_parameter*(x2-xc);
					y3 = y2+distortion_parameter*(y2-yc);
				}
			}
			// ORIENTATION CORRECTION
			vector<float> corrected_orientation = orientation_update(point2d<double>(x2,y2),
																									sp_contours_modifed.get_seno()[q],
																									sp_contours_modifed.get_coseno()[q],
																									it_ldm);
			// WE ESTIMATE DE ANGLE INTERVAL
			if(corrected_orientation[1]>=0)
			{
				angle = (int) ((ami_pi-atan2(corrected_orientation[1],
												corrected_orientation[0]))/paso_angle);
			}
			else
			{
				angle = (int) ((ami_pi-atan2(-corrected_orientation[1],
												-corrected_orientation[0]))/paso_angle);
			}
			if(angle==height_score) angle=height_score-1;

			int l_min=(int) (angle-angle_point_orientation_max_difference/angle_resolution);
			int l_max=(int) (angle+angle_point_orientation_max_difference/angle_resolution);
			int id = 2;
			for(int l=l_min;l<(int) angle;l++)
			{
				if(l<0)
				{
					distancia=-coseno[height_score+l]*y3-seno[height_score+l]*x3;
					for(int nd=-id;nd<=id;nd++)
					{
						d= width_score/2+((int)(distancia/distance_resolution+0.5))+nd;
						if(d>=0 && d<width_score)
						{
							double distancia_recta = fabs((distancia+coseno[height_score+l]*y3+seno[height_score+l]*x3)+nd*distance_resolution);;
							score_k[d][height_score+l]+=1./(1.+distancia_recta+(angle-l)*angle_resolution);
						}
					}
				}
				else
				{
					distancia=-coseno[l]*y3-seno[l]*x3;
					for(int nd=-id;nd<=id;nd++)
					{
						d= width_score/2+((int)(distancia/distance_resolution+0.5))+nd;
						if(d>=0 && d<width_score)
						{
							double distancia_recta = fabs((distancia+coseno[l]*y3+seno[l]*x3)+nd*distance_resolution);
							score_k[d][l]+=1./(1.+distancia_recta+(angle-l)*angle_resolution);
						}
					}
				}
			}
			for(int l=(int) angle;l<l_max;l++)
			{
				if(l>=height_score)
				{
					distancia=-coseno[l-height_score]*y3-seno[l-height_score]*x3;
					for(int nd=-id;nd<=id;nd++)
					{
						d= width_score/2+((int)(distancia/distance_resolution+0.5))+nd;
						if(d>=0 && d<width_score)
						{
							double distancia_recta = fabs((distancia+coseno[l-height_score]*y3-seno[l-height_score]*x3)+nd*distance_resolution);
							score_k[d][l-height_score]+=1./(1.+distancia_recta+(l-angle)*angle_resolution);
						}
					}
				}
				else
				{
					distancia=-coseno[l]*y3-seno[l]*x3;
					for(int nd=-id;nd<=id;nd++)
					{
						d= width_score/2+((int)(distancia/distance_resolution+0.5))+nd;
						if(d>=0 && d<width_score)
						{
							double distancia_recta = fabs((distancia+coseno[l]*y3+seno[l]*x3)+nd*distance_resolution);
							score_k[d][l]+=1./(1.+distancia_recta+(l-angle)*angle_resolution);
						}
					}
				}
			}
		} // END INDEX

		// WE CREATE TWO VECTORS WITH THE POSITIONS INSIDE THE SCORE MATRIX WITH A
		// SCORE HIGHER THAN 10
    vector<int> i_pos(width_score*height_score,0);
    vector<int> j_pos(width_score*height_score,0);
    int pind=0;
    for(int i=0;i<height_score;i++)
    {
      for(int j=0;j<width_score;j++)
      {
        if(score_k[j][i]>10)
        {
          i_pos[pind] = i;
          j_pos[pind] = j;
          pind++;
        }
      }
    }
    i_pos.resize(pind);
    j_pos.resize(pind);

    // WE SELECT THE MAXIMUM OF THE SLICE
    votation_score=0.;
    vector<int> m(nlineas_plus);
    vector<int> n(nlineas_plus);
    vector<int> m_max(nlineas_plus);
    vector<int> n_max(nlineas_plus);
    vector<line_points> lines(nlineas_plus);
    float first_max_score = 0;
    for(int l=0;l<nlineas_plus;l++)
    {
      double shared_max=0;
// 			#ifdef _OPENMP
//       int vti[omp_get_num_threads()],
//           vtj[omp_get_num_threads()];
// 			#else
			int vti[1]={0},vtj[1]={0};
// 			#endif

// 			#ifdef _OPENMP
//       #pragma omp parallel
// 			#endif
      {
        double max_score=0;
// 				#ifdef _OPENMP
//         #pragma omp for nowait
// 				#endif
        for(int ind_pos=0; ind_pos<(int)i_pos.size(); ind_pos++)
        {
          int i = i_pos[ind_pos];
          int j = j_pos[ind_pos];
          int iu,id,jl,jr;
					iu = id = jl = jr = 0;
          (i==0)              ? iu=height_score-1 : iu=i-1;
          (i==height_score-1) ? id=0              : id=i+1;
          (j==0)              ? jl=width_score-1  : jl=j-1;
          (j==width_score-1)  ? jr=0              : jr=j+1;
          if(score_k[j][i]>=max_score && score_k[j][i]!=first_max_score &&
               (score_k[j][i]>=score_k[jl][iu] || score_k[jl][iu]==first_max_score) &&
               (score_k[j][i]>=score_k[j][iu]  || score_k[j][iu]==first_max_score)  &&
               (score_k[j][i]>=score_k[jr][iu] || score_k[jr][iu]==first_max_score) &&
               (score_k[j][i]>=score_k[jl][i]  || score_k[jl][i]==first_max_score)  &&
               (score_k[j][i]>=score_k[jr][i]  || score_k[jr][i]==first_max_score)  &&
               (score_k[j][i]>=score_k[jl][id] || score_k[jl][id]==first_max_score) &&
               (score_k[j][i]>=score_k[j][id]  || score_k[j][id]==first_max_score)  &&
               (score_k[j][i]>=score_k[jr][id] || score_k[jr][id]==first_max_score))
            {
							max_score = score_k[j][i];
// 							#ifdef _OPENMP
// 							vti[omp_get_thread_num()] = i;
// 							vtj[omp_get_thread_num()] = j;
// 							#else
							vti[0] = i;
							vtj[0] = j;
// 							#endif
            }
        } // END IND FOR
				// EVERY THREAD UPDATES THE MAXIMUM VALUE AND ITS POSITION
// 				#ifdef _OPENMP
//         #pragma omp critical
// 				#endif
        {
          shared_max = max(shared_max,max_score);
          if(shared_max==max_score && shared_max!=0)
          {
// 						#ifdef _OPENMP
//             m[l] = vti[omp_get_thread_num()];
//             n[l] = vtj[omp_get_thread_num()];
// 						#else
						m[l] = vti[0];
						n[l] = vtj[0];
// 						#endif
          }
        }
      } // END PRAGMA OMP PARALLEL

      lines[l].set_a((float)seno[m[l]]);
      lines[l].set_b((float)coseno[m[l]]);
      lines[l].set_c((float)(n[l]-width_score/2)*distance_resolution);

      if(shared_max>=20) votation_score+=shared_max;

      //We take the first maximum value+0.01
      if(l==0)
        first_max_score = shared_max+0.01;

      //WE SET TO FIRST_MAX_SCORE A NEIGHBORHOOD OF MAXIMUM SCORE POINT
      int line_angle = m[l];
      int d0 = n[l];
      int alpha_min=(int) (line_angle-angle_point_orientation_max_difference/angle_resolution);
      int alpha_max=(int) (line_angle+angle_point_orientation_max_difference/angle_resolution);
      int Ri = distance_point_line_max+2;
      int d_min = d0-Ri;
      int d_max = d0+Ri;
      bool first_pos=true;
      int kang=0;
      score_k[n[l]][m[l]]=first_max_score;
      // FROM ANGLE TO ANGLE_MAX
      for(int a=line_angle; a<alpha_max; a++)
      {
        int ang = a;
        float suma=0., suma_l=0., suma_r=0.;
        float media=0., media_l=0., media_r=0.;
        int n_ele=0;
        if(a>=height_score)
        {
          ang = a-height_score;
          if(first_pos)
          {
            d0 = width_score/2 - (d0 - width_score/2);
            first_pos = false;
          }
        }
        // SUMS
        // LEFT (CENTERED ON d0-1)
        if(d0>1)
        {
          suma_l = score_k[d0-2][ang] + score_k[d0-1][ang] + score_k[d0][ang];
          n_ele = ((d0-1)+Ri)-((d0-1)-Ri)+1;
          if (n_ele!=0) media_l = suma_l/n_ele;
        }
        // RIGHT (CENTERED ON d0+1)
        if(d0<width_score-2)
        {
          suma_r = score_k[d0][ang] + score_k[d0+1][ang] + score_k[d0+2][ang];
          n_ele = ((d0+1)+Ri)-((d0+1)-Ri)+1;
          if(n_ele!=0) media_r = suma_r/n_ele;
        }
        // CENTER (CENTERED ON d0)
        if(d0>0 && d0<width_score-1)
        {
          suma = score_k[d0-1][ang] + score_k[d0][ang] + score_k[d0+1][ang];
          n_ele = (d0+Ri)-(d0-Ri)+1;
          if(n_ele!=0) media = suma/n_ele;
        }
        // WE UPDATE D0 ACCORDING TO THE HIGHEST VALUE
        float val_l = suma_l*media_l, val = suma*media, val_r = suma_r*media_r;
        float max_val = max(val_l,max(val,val_r));
        if(max_val == val_l)
          d0 = d0-1;
        if(max_val == val_r)
          d0 = d0+1;

        d_min = d0-Ri;
        d_max = d0+Ri;
        // WE CHECK THE LIMITS
        if(d_min<0) d_min = 0;
        if(d_max>=width_score) d_max=width_score-1;
        // WE SET THE INTERVAL D_MIN - D_MAX TO FIRST_MAX_SCORE
        for(int dist=d_min; dist<=d_max; dist++)
          score_k[dist][ang] = first_max_score;

        // UPDATING Ri
        int f=5;
        if((kang*angle_resolution*f)<(distance_point_line_max+2))
          Ri = distance_point_line_max+2;
        else
          Ri = kang*angle_resolution*f;
        kang++;
      }
      // WE RESET Ri D0 AND MAXIMUM AND MINIMUM DISTANCE
      Ri = distance_point_line_max+2;
      d0 = n[l];
      d_min = d0-Ri;
      d_max = d0+Ri;
      bool first_neg = true;
      kang = 0;
      // FROM ANGLE TO ANGLE_MIN
      for(int a = line_angle; a>=alpha_min; a--)
      {
        int ang = a;
        float suma=0., suma_l=0., suma_r=0.;
        float media=0., media_l=0., media_r=0.;
        int n_ele=0;
        if(a<0)
        {
          ang = height_score+a;
          if(first_neg)
          {
            d0 = width_score/2 + (width_score/2-d0);
            first_neg=false;
          }
        }
        // SUMS
        // LEFT (CENTERED ON d0-1)
        if(d0>1)
        {
          suma_l = score_k[d0-2][ang] + score_k[d0-1][ang] + score_k[d0][ang];
          n_ele = ((d0-1)+Ri)-((d0-1)-Ri)+1;
          if(n_ele!=0) media_l = suma_l/n_ele;
        }
        // RIGHT (CENTERED ON d0+1)
        if(d0<width_score-2)
        {
          suma_r = score_k[d0][ang] + score_k[d0+1][ang] + score_k[d0+2][ang];
          n_ele = ((d0+1)+Ri)-((d0+1)-Ri)+1;
          if(n_ele!=0) media_r = suma_r/n_ele;
        }
        // CENTER (CENTERED ON d0)
        if(d0>0 && d0<width_score-1)
        {
          suma = score_k[d0-1][ang] + score_k[d0][ang] + score_k[d0+1][ang];
          n_ele = (d0+Ri)-(d0-Ri)+1;
          if(n_ele!=0) media = suma/n_ele;
        }
        // WE UPDATE D0 ACCORDING TO THE HIGHEST VALUE
        float val_l = suma_l*media_l, val = suma*media, val_r = suma_r*media_r;
        float max_val = max(val_l,max(val,val_r));
        if(max_val == val_l)
          d0 = d0-1;
        if(max_val == val_r)
          d0 = d0+1;

        d_min = d0-Ri;
        d_max = d0+Ri;
        // WE CHECK THE LIMITS
        if(d_min<0) d_min = 0;
        if(d_max>=width_score) d_max=width_score-1;
        // WE SET THE INTERVAL D_MIN - D_MAX TO FIRST_MAX_SCORE
        for(int dist=d_min; dist<=d_max; dist++)
          score_k[dist][ang] = first_max_score;

        // UPDATING Ri
        int f=5;
        if((kang*angle_resolution*f)<(distance_point_line_max+2))
          Ri = distance_point_line_max+2;
        else
          Ri = kang*angle_resolution*f;
        kang++;
      }
    } // END LINES LOOP
    
    // WE SELECT THE MAXIMUN DISTORTION LEVEL
    if(votation_score>max_votation_score)
    {
      max_votation_score=votation_score;

      image_primitive.set_lines(lines);
      if(ini_ldm.get_type()==DIVISION)
      {
        double p = initial_distortion_parameter +
                   (k*distortion_parameter_resolution);
        best_distortion_parameter= -p / (dmi + dmi*p);
      }
      else
      {
        best_distortion_parameter=(initial_distortion_parameter+
              k*distortion_parameter_resolution)/max_norm;
      }
      
      for(int l=0;l<nlineas_plus;l++){n_max[l]=n[l]; m_max[l]=m[l];}
    }
    ami2_free2d(score_k);
  } // END K LOOP
  
  // IMAGE_PRIMITIVES OBJECTS FOR CORRECTED POINTS AND ORIGINAL POINTS
  image_primitives image_primitive_corrected, image_primitive_original;
  image_primitive_corrected.set_lines(image_primitive.get_lines());
  image_primitive_original.set_lines(image_primitive.get_lines());
  
	// WE FILL THE IMAGE PRIMITIVE LENS DISTORTION MODEL
  lens_distortion_model ld;
  if(best_distortion_parameter!=0.)
	{
    ld.set_distortion_center(point2d<double>(xc,yc));
    ld.get_d().resize(2);
    ld.get_d()[0]=1.;
    ld.get_d()[1]=best_distortion_parameter;
    ld.set_type(ini_ldm.get_type());
    image_primitive.set_distortion(ld);
    image_primitive_original.set_distortion(ld);
  }
  else
  {
    ld.set_distortion_center(point2d<double>(xc,yc));
    ld.get_d().resize(2);
    ld.get_d()[0]=1.;
    ld.get_d()[1]=0.;
    ld.set_type(ini_ldm.get_type());
    image_primitive.set_distortion(ld);
    image_primitive_original.set_distortion(ld);
  }

  // WE COMPUTE THE POINTS OF THE LINE POINTS
  // FOLLOWING THE DISTANCE OF THE POINTS TO THE LINE
  double dot_product_min=cos(ami_pi*angle_point_orientation_max_difference/180.);
  for(k=0;k<nlineas_plus;k++) image_primitive.get_lines()[k].get_points().clear();
	for(int ind=0; ind<(int)index.size(); ind++)
	{
		int q=index[ind];
		x2=sp_contours_modifed.get_x()[q];
		y2=sp_contours_modifed.get_y()[q];
		double x2ori = subpixel_contours.get_x()[q];
		double y2ori = subpixel_contours.get_y()[q];
		point2d<double> p2(x2,y2);
		point2d<double> p2ori(x2ori,y2ori);
		point2d<double> p2d=ld.evaluation(p2);
		vector<float> corrected_orientation = orientation_update(p2,
																							sp_contours_modifed.get_seno()[q],
																							sp_contours_modifed.get_coseno()[q],
																							image_primitive.get_distortion());
		for(k=0;k<nlineas_plus;k++)
		{
			if(fabs(image_primitive.get_lines()[k].get_b()*corrected_orientation[0]
					-image_primitive.get_lines()[k].get_a()*corrected_orientation[1] )<
					dot_product_min) continue;
			if(fabs(image_primitive.get_lines()[k].evaluation(p2d))<
			(distance_point_line_max+0.5*distance_resolution))
			{
				 image_primitive.get_lines()[k].get_points().push_back(p2);
				 image_primitive_corrected.get_lines()[k].get_points().push_back(p2d);
				 image_primitive_original.get_lines()[k].get_points().push_back(p2ori);
				 k=nlineas_plus;
			}
		}
	}

  // DEBUGGING OF THE PRIMITIVES THROUGH THE ORIENTATION OF THE POINTS
  for(int il=0; il<(int)image_primitive_original.get_lines().size(); il++)
  {
    int pos_count=0, neg_count=0;
    for(int ip=0; ip<(int)image_primitive_original.get_lines()[il].get_points().size(); ip++)
    {
      point2d<double> current_point = image_primitive_original.get_lines()[il].get_points()[ip];
      int pos  = width*round(current_point.y) + round(current_point.x);
      double a = image_primitive_original.get_lines()[il].get_a();
      double b = image_primitive_original.get_lines()[il].get_b();
      double orientation_sign = subpixel_contours.get_coseno()[pos]*b +
                                subpixel_contours.get_seno()[pos]*-a;
      if(orientation_sign>0)
        pos_count++;
      if(orientation_sign<0)
        neg_count++;
    }
    if(pos_count!=0 && neg_count!=0)
    {
      for(int ip=0; ip<(int)image_primitive_original.get_lines()[il].get_points().size(); ip++)
      {
        point2d<double> current_point = image_primitive_original.get_lines()[il].get_points()[ip];
        int pos  = width*current_point.y + current_point.x;
        double a = image_primitive_original.get_lines()[il].get_a();
        double b = image_primitive_original.get_lines()[il].get_b();
        double orientation_sign = subpixel_contours.get_coseno()[pos]*b +
                                  subpixel_contours.get_seno()[pos]*-a;
        if(pos_count>neg_count) // ERASE NEGATIVES
        {
          if(orientation_sign<0)
          {
            image_primitive_original.get_lines()[il].get_points().erase(image_primitive_original.get_lines()[il].get_points().begin()+ip);
            image_primitive.get_lines()[il].get_points().erase(image_primitive.get_lines()[il].get_points().begin()+ip);
            image_primitive_corrected.get_lines()[il].get_points().erase(image_primitive_corrected.get_lines()[il].get_points().begin()+ip);
            ip--;
          }
        }
        else // ERASE POSITIVES
        {
          if(orientation_sign>0)
          {
            image_primitive_original.get_lines()[il].get_points().erase(image_primitive_original.get_lines()[il].get_points().begin()+ip);
            image_primitive.get_lines()[il].get_points().erase(image_primitive.get_lines()[il].get_points().begin()+ip);
            image_primitive_corrected.get_lines()[il].get_points().erase(image_primitive_corrected.get_lines()[il].get_points().begin()+ip);
            ip--;
          }
        }
      }
    }
  }

  // WE RECOMPUTE THE LINE EQUATIONS
  if (lens_distortion_estimation==true) 
	{
    for(i=0;i<(int)image_primitive.get_lines().size();i++)
		{
      if(image_primitive.get_lines()[i].get_points().size()>2)
      {
        image_primitive_corrected.get_lines()[i].points_to_equation();
        double a,b,c;
        image_primitive_corrected.get_lines()[i].get_abc(a,b,c);
        image_primitive.get_lines()[i].set_abc(a,b,c);
        image_primitive_original.get_lines()[i].set_abc(a,b,c);
      }
    }
  }
  else
	{
    for(i=0;i<(int)image_primitive.get_lines().size();i++)
      if(image_primitive.get_lines()[i].get_points().size()>2)
      {
        image_primitive.get_lines()[i].points_to_equation();
        double a,b,c;
        image_primitive.get_lines()[i].get_abc(a,b,c);
        image_primitive_corrected.get_lines()[i].set_abc(a,b,c);
        image_primitive_original.get_lines()[i].set_abc(a,b,c);
      }

  }

  // WE REMOVE LINES WITH A SMALL NUMBER OF POINTS AND WE JOINT LINES
  // WHICH ARE TOO CLOSE
  dot_product_min=cos(2.*ami_pi*angle_point_orientation_max_difference/180.);
  if(dot_product_min<0.95) dot_product_min=0.95;
  dot_product_min=0.99;

  float min_line_points=0.05*image_primitive.get_lines()[0].get_points().size();
  if(min_line_points<20) min_line_points=20;
  float distance_line_line_min=2.*distance_point_line_max;
  if(distance_line_line_min>10) distance_line_line_min=10.;

  for(k=0;k<(int)image_primitive.get_lines().size();k++)
	{
    if(image_primitive.get_lines()[k].get_points().size()<min_line_points){
      image_primitive.get_lines().erase(image_primitive.get_lines().begin()+k);
      image_primitive_corrected.get_lines().erase(image_primitive_corrected.get_lines().begin()+k);
      image_primitive_original.get_lines().erase(image_primitive_original.get_lines().begin()+k);
      k--;
      continue;
    }
    for(l=k+1;l<(int)image_primitive.get_lines().size();l++)
		{
      // WE CHECK THE NUMBER OF POINTS
      if(!image_primitive.get_lines()[l].get_points().size()>0)
        continue;
      double paso=
         image_primitive.get_lines()[k].get_a()*image_primitive.get_lines()[l].get_a()+
         image_primitive.get_lines()[k].get_b()*image_primitive.get_lines()[l].get_b();
      if(fabs(paso)>dot_product_min)
			{
        // WE CHECK THE ORIENTATION OF THE FIRST POINT OF EACH PRIMITIVE
        double a1=image_primitive_corrected.get_lines()[k].get_a(),
               a2=image_primitive_corrected.get_lines()[l].get_a(),
               b1=image_primitive_corrected.get_lines()[k].get_b(),
               b2=image_primitive_corrected.get_lines()[l].get_b();
        point2d<double> point=image_primitive_original.get_lines()[k].get_points()[0];
        int pos = width*point.y + point.x;
        double sign1 = b1*subpixel_contours.get_coseno()[pos]-a1*subpixel_contours.get_seno()[pos];
        point = image_primitive_original.get_lines()[l].get_points()[0];
        pos = width*point.y + point.x;
        double sign2 = b2*subpixel_contours.get_coseno()[pos]-a2*subpixel_contours.get_seno()[pos];
        if(!(sign1>0&&sign2>0) && !(sign1<0&&sign2<0))
           continue;
        // WE CHECK THE AVERAGE DISTANCE OF THE POINTS TO THE LINE
        paso=0;
        int m;
        for(m=0;m<(int)image_primitive.get_lines()[l].get_points().size();m++)
				{
          point2d<double> p2d=image_primitive_corrected.get_lines()[l].get_points()[m];
          double dist=fabs(image_primitive.get_lines()[k].evaluation(p2d));
          paso+=dist;
          if(dist>(5.*distance_point_line_max)){
            m=image_primitive.get_lines()[l].get_points().size()+1;
            break;
          }
        }
        if(m==((int)image_primitive.get_lines()[l].get_points().size()+1)) continue;
        if((paso/(int)image_primitive.get_lines()[l].get_points().size())<20)
				{
          // WE CHECK THE DISTANCE OF THE POINTS OF A LINE TO THE POINTS OF THE
          // OTHER LINE. THIS DISTANCE SHOULD BE BIG TO AVOID JOIN PARALLEL LINES
          paso=0;
          for(m=0;m<(int)image_primitive.get_lines()[l].get_points().size();m++)
					{
            point2d<double> p2d=image_primitive_corrected.get_lines()[l].get_points()[m];
            double a,b,c;
            image_primitive.get_lines()[k].get_abc(a,b,c);
            double d = a*p2d.x + b*p2d.y + c;
            point2d<double> np2d(p2d.x-d*a, p2d.y-d*b);
            paso+=image_primitive_corrected.get_lines()[k].distance(np2d);
          }
          if( (paso/image_primitive_corrected.get_lines()[l].get_points().size())>10)
					{
            // WE ADD THE POINTS OF THE LINE TO THE LINE POINT STRUCTURE
            image_primitive.get_lines()[k].get_points().insert(
              image_primitive.get_lines()[k].get_points().end(),
              image_primitive.get_lines()[l].get_points().begin(),
              image_primitive.get_lines()[l].get_points().end());

            image_primitive_corrected.get_lines()[k].get_points().insert(
              image_primitive_corrected.get_lines()[k].get_points().end(),
              image_primitive_corrected.get_lines()[l].get_points().begin(),
              image_primitive_corrected.get_lines()[l].get_points().end());

            image_primitive_original.get_lines()[k].get_points().insert(
              image_primitive_original.get_lines()[k].get_points().end(),
              image_primitive_original.get_lines()[l].get_points().begin(),
              image_primitive_original.get_lines()[l].get_points().end());
            // WE REMOVE THE LINE POINTS STRUCTURE
            image_primitive.get_lines().erase(image_primitive.get_lines().begin()+l);
            image_primitive_corrected.get_lines().erase(image_primitive_corrected.get_lines().begin()+l);
            image_primitive_original.get_lines().erase(image_primitive_original.get_lines().begin()+l);
            l--;
          }
        }
      }
    }
  }

  if (lens_distortion_estimation==true) {
    for(i=0;i<(int)image_primitive.get_lines().size();i++){
			image_primitive_corrected.get_lines()[i].points_to_equation();
			double a,b,c;
			image_primitive_corrected.get_lines()[i].get_abc(a,b,c);
			image_primitive.get_lines()[i].set_abc(a,b,c);
			image_primitive_original.get_lines()[i].set_abc(a,b,c);
    }
  }
  else{
    for(i=0;i<(int)image_primitive.get_lines().size();i++){
			image_primitive.get_lines()[i].points_to_equation();
			double a,b,c;
			image_primitive.get_lines()[i].get_abc(a,b,c);
			image_primitive_corrected.get_lines()[i].set_abc(a,b,c);
			image_primitive_original.get_lines()[i].set_abc(a,b,c);
    }
  }
	
  // WE FREE THE MEMORY
  free(seno); free(coseno);
  // WE SET THE ORIGINAL IMAGE PRIMITIVES AND RETURN THE NUMBER OF VOTES
  image_primitive = image_primitive_original;
  return max_votation_score;
}

