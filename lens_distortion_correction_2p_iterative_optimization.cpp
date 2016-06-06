/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


/**
 * @file lens_distortion_correction_2p_iterative_optimization.cpp
 * @brief distortion correction using ....
 *
 * @author Luis Alvarez <lalvarez@dis.ulpgc.es> and Daniel Santana-Cedrés <dsantana@ctim.es>
 */


//Included libraries
#include "ami_image/image.h"
#include "ami_filters/filters.h"
#include "ami_primitives/subpixel_image_contours.h"
#include "ami_primitives/line_extraction.h"
#include "ami_primitives/image_primitives.h"
#include "ami_lens_distortion/lens_distortion_procedures.h"
#include "ami_utilities/utilities.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

//------------------------------------------------------------------------------
//Minimizing the energy according to the type of the lens distortion model
//and the center optimization
double energy_minimization(lens_distortion_model &ldm, image_primitives &ip,
                           int w, int h, bool opt_center)
{
  int model_type = ldm.get_type();
  vector<bool> vtf(4); vtf[0]=vtf[1]=true; vtf[2]=vtf[3]=false;
  vector<bool> vtt(4,true);
  double error=0., tf_error=0., tt_error=0.;
  lens_distortion_model tf_ldm, tt_ldm;
  if(!opt_center)
  {
    if(model_type==POLYNOMIAL)
    {
      error = model_center_estimation_2p_polynomial(ldm.get_distortion_center(),
                                                    ip.get_lines(), ldm, w, h,
                                                    vtf);
    }
    else
    {
      error = model_center_estimation_2p_quotient(ldm.get_distortion_center(),
                                                  ip.get_lines(), ldm, w, h,
                                                  vtf);
    }
  }
  else
  {
    if(model_type==POLYNOMIAL)
    {
      tf_error = model_center_estimation_2p_polynomial(ldm.get_distortion_center(),
                                                       ip.get_lines(), ldm, w,
                                                       h, vtf);
      tf_ldm = ldm;
      tt_error = model_center_estimation_2p_polynomial(ldm.get_distortion_center(),
                                                       ip.get_lines(), ldm, w,
                                                       h, vtt);
      tt_ldm = ldm;
    }
    else
    {
      tf_error = model_center_estimation_2p_quotient(ldm.get_distortion_center(),
                                                     ip.get_lines(), ldm, w, h,
                                                     vtf);
      tf_ldm = ldm;
      tt_error = model_center_estimation_2p_quotient(ldm.get_distortion_center(),
                                                    ip.get_lines(), ldm, w, h,
                                                    vtt);
      tt_ldm = ldm;
    }
    
    if((fabs(tf_ldm.get_distortion_center().x-
              tt_ldm.get_distortion_center().x)>0.2*w) ||
         (fabs(tf_ldm.get_distortion_center().y-
              tt_ldm.get_distortion_center().y)>0.2*h) ||
         (check_invertibility(tt_ldm, w, h) == false)){
      ldm = tf_ldm;
      error = tf_error;
    }
    else
    {
      ldm = tt_ldm;
      error = tt_error;
    }
  }
  
  return error;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Iterative optimization
double iterative_optimization(
  ami::subpixel_image_contours &contours,
  image_primitives &i_primitives,
  float distance_point_line_max_hough,
  int max_lines,
  float angle_resolution,
  float distance_resolution,
  float distortion_parameter_resolution,
  float angle_point_orientation_max_difference,
  bool opt_center,
  int width, int height
)
{
  double final_error=0.;
  
  //We initiallize the previous model, the best one and the previous set of primitives
  i_primitives.get_distortion().get_d().resize(3);
  i_primitives.get_distortion().get_d()[2] = 0.;
  lens_distortion_model previous_model = i_primitives.get_distortion();
  lens_distortion_model best_model = previous_model;
  image_primitives previous_ip = i_primitives;
  //Tolerance for convergence
  double TOL = 1e-2;
  //Number of converngence iterations
  int convergence_iterations = 0;
  //Fails counter
  int fail_count = 0;
  //Number of points: current, best and next
  int num_points = count_points(i_primitives);
  int best_num_points = num_points;
  int next_num_points = num_points+TOL;
  lens_distortion_model ldm_tf, ldm_tt;
  //We apply the process until the number of points is not significantly greater
  //or until the process fails three times 
  while((next_num_points >= (num_points*(1+TOL))) || (fail_count<3))
  {
    double error = energy_minimization(previous_model, i_primitives, width,
                                       height, opt_center);
    if(check_invertibility(previous_model, width, height) == false)
    {
      if(fail_count==3)
        break;
    }
      
    i_primitives.clear();
    //CALL TO IMPROVED HOUGH WITH THE MODEL COMPUTED BEFORE
    line_equation_distortion_extraction_improved_hough(
                contours,i_primitives,distance_point_line_max_hough,
                max_lines,angle_resolution,
                distance_resolution, 0., 0.,
                distortion_parameter_resolution,
                angle_point_orientation_max_difference,
                true, previous_model);
    
    int local_num_points = count_points(i_primitives);
    if(local_num_points>next_num_points)
    {
      //We update the primitives only if the result is better
      if(local_num_points>best_num_points)
      {
        previous_ip = i_primitives;
        best_num_points = local_num_points;
        best_model = previous_model;
        final_error = error;
      }
    }
    else
    {
      fail_count++;
      if(fail_count==3) //THE LIMIT OF FAILS IS 3
        break;
    }
    num_points = next_num_points;
    next_num_points = local_num_points;
    convergence_iterations++;
  }
  //We take the last and best image primitives object and model
  i_primitives = previous_ip;
  i_primitives.set_distortion(best_model);
  
  //We return the average error
  return (final_error/count_points(i_primitives));
}
//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  if(argc!=13){
    print_function_syntax_lens_distortion_correction_2p_iterative_optimization();
    exit(EXIT_FAILURE);
  }

  if(check_params_lens_distortion_correction_2p_iterative_optimization(argv) != 0){
    manage_failure(argv,0);
    exit(EXIT_SUCCESS);
  }
  
  //We read the input image and some variables are initialized
  ami::image<unsigned char> input(argv[1]); // input image
  int width = input.width(), height = input.height();//input image dimensions
  int size_ = width*height; // image size
  ami::image<unsigned char> gray(width,height,1,0);//gray-level image to call canny
  ami::image<unsigned char> edges(width,height,1,0);//image to store edge information
  float canny_high_threshold = atof(argv[5]);// high threshold for canny detector
  float initial_distortion_parameter=atof(argv[6]);//left side of allowed distortion parameter interval
  float final_distortion_parameter = atof(argv[7]);//Hough parameter
  float distance_point_line_max_hough = atof(argv[8]);//Hough parameter
  //maximum difference allowed (in degrees) between edge and line orientation
  float angle_point_orientation_max_difference = atof(argv[9]);
  string tmodel(argv[10]);
  string s_opt_c(argv[11]);
  bool opt_center = false;
  (s_opt_c == string("True")) ? opt_center = true : opt_center = false;

  //Converting the input image to gray level
  for(int i=0; i<size_; i++)
    gray[i] = 0.3*input[i] + 0.59*input[i+size_] + 0.11*input[i+size_*2];

  //ALGORITHM STAGE 1 : Detecting edges with Canny
  cout << "Detecting edges with Canny..." << endl;
  float canny_low_threshold = 0.7; //default value for canny lower threshold
  ami::subpixel_image_contours contours=canny(gray,edges,canny_low_threshold,
                                              canny_high_threshold);

  //We create writable 3 channel images for edges and gray level
  ami::image<unsigned char> edges3c(width, height, 3, 255),
                            gray3c(width, height, 3);
  gray3c=gray;
    
  //We clean the contours
  int neighborhood_radius=2; //radius of neighborhood to take into account
  int min_neighbor_points=4; //min number of contour points in a neighborhood
  double min_orientation_value=0.95; /*  min average scalar product of 
                                         neigborhood point orientation */
  int min_distance_point=1; //minimum distance between contour points
  contours.clean(
    neighborhood_radius,
    min_neighbor_points,
    min_orientation_value,
    min_distance_point
  );
  vector<int> index = contours.get_index();
  for(int i =0; i<(int)index.size(); i++)
  {
    edges3c[index[i]]=0;
    edges3c[index[i]+size_]=0;
    edges3c[index[i]+2*size_]=0;
  }
  //Writing Canny detector output after the cleaning process
  edges3c.write(argv[2]);
  cout << "...edges detected" << endl;

  //ALGORITHM STAGE 2 : Detecting lines with improved_hough_quotient
  cout << "Detecting lines with improved Hough and " <<  tmodel <<" model..." << endl;
  image_primitives i_primitives;//object to store output edge line structure
  int max_lines=100; //maximun number of lines estimated
  float angle_resolution=0.1; // angle discretization step (in degrees)
  float distance_resolution=1.; // line distance discretization step
  float distortion_parameter_resolution=0.1;//distortion parameter discretization step
  lens_distortion_model ini_ldm;
  if(tmodel == string("pol"))
    ini_ldm.set_type(POLYNOMIAL);
  else
    ini_ldm.set_type(DIVISION);
  //we call 3D Hough line extraction
  line_equation_distortion_extraction_improved_hough(
    contours,
    i_primitives,
    distance_point_line_max_hough,
    max_lines,
    angle_resolution,
    distance_resolution,
    initial_distortion_parameter,
    final_distortion_parameter,
    distortion_parameter_resolution,
    angle_point_orientation_max_difference,
    true,
    ini_ldm
  );
  
  //ALGORITHM STAGE 3 : We apply the iterative optimization process
  double final_error = iterative_optimization(
                          contours,
                          i_primitives,
                          distance_point_line_max_hough,
                          max_lines,
                          angle_resolution,
                          distance_resolution,
                          distortion_parameter_resolution,
                          angle_point_orientation_max_difference,
                          opt_center, width, height);
                          
  //We check if the iterative optimization process finishes properly
  if(i_primitives.get_lines().size()==0){
    manage_failure(argv,0);
    exit(EXIT_SUCCESS);
  }
  cout << "...lines detected: " << i_primitives.get_lines().size() <<
          " with " << count_points(i_primitives) << " points" << endl;
  
  //Drawing the detected lines on the original image to illustrate the results
  drawHoughLines(i_primitives,gray3c);
  gray3c.write(argv[3]);
  
  //ALGORITHM STAGE 4 : Correcting the image distortion using the estimated model
  if(i_primitives.get_distortion().get_d().size()>0){
    cout << "Correcting the distortion..." << endl;
    if(i_primitives.get_distortion().get_type()==DIVISION)
    {
      ami::image<unsigned char> undistorted =undistort_quotient_image_inverse(
        input, // input image
        i_primitives.get_distortion(), // lens distortion model
        3.0 // integer index to fix the way the corrected image is scaled to fit input size image
      );
      //Writing the distortion corrected image
      undistorted.write(argv[4]);
    }
    else
    {
      lens_distortion_model ldm = i_primitives.get_distortion();
      int vs=0;
      (ldm.get_d().size()==2) ? vs=3 : vs=5;
      double *a = new double[vs];
      for(int i=0, ldmind=0; i<vs; i++)
      {
        if(i%2==0)
        {
          a[i] = ldm.get_d()[ldmind];
          ldmind++;
        }
        else
        {
          a[i] = 0.;
        }
      }
      ami::image<unsigned char> undistorted = undistort_image_inverse_fast(
        input,
        vs-1,
        a,
        i_primitives.get_distortion().get_distortion_center(),
        2.0
      );
      delete []a;
      //Writing the distortion corrected image
      undistorted.write(argv[4]);
    }
    cout << "...distortion corrected." << endl;
  }
  
  // WRITING OUTPUT TEXT DOCUMENTS
  // writing in a file the lens distortion model and the lines and associated points
  i_primitives.write(argv[12]);
  // writing function parameters and basic outputs :
  ofstream fs("output.txt");//Output file
  fs << "Selected parameters:" << endl;
  fs << "\t High Canny's threshold: " << argv[5] << endl;
  fs << "\t Initial normalized distortion parameter: " << argv[6] << endl;
  fs << "\t Final normalized distortion parameter: " << argv[7] << endl;
  fs << "\t Maximum distance between points and line: " << argv[8] << endl;
  fs << "\t Maximum differente between edge point and line orientations: "
     << argv[9]  << endl;
  fs << "\t Model applied: " << argv[10] << endl;
  fs << "\t Center optimization: " << argv[11] << endl;
  fs << "-------------------------" << endl;
  fs << "Results: " << endl;
  fs << "\t Number of detected lines: " << i_primitives.get_lines().size()
     << endl;
  int count = count_points(i_primitives);
  fs << "\t Total amount of line points: " << count << endl;
  fs << "\t Distortion center: (" << i_primitives.get_distorsion_center().x <<
        ", " << i_primitives.get_distorsion_center().y << ")" << endl; 
  double p1=0., p2=0.;
  bool is_division=false;
  (tmodel == string("pol")) ? is_division = false : is_division = true;
  compute_ps(p1,p2,i_primitives.get_distortion(),width,height,is_division);
  fs << "\t Estimated normalized distortion parameters: p1 = " << p1
     << " p2 = " << p2 << endl;
  fs << "\t Average squared error distance in pixels between line and associated points = " << final_error << endl;

  fs.close();
  exit(EXIT_SUCCESS);
}
