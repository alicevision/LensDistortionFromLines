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
#include <dirent.h>
#include <stdexcept>
#include <limits>

using namespace std;

//------------------------------------------------------------------------------

/**
 * @brief Minimizing the energy according to the type of the lens distortion
 *        model and the center optimization
 * 
 * @param[in,out] ldm: lens distortion model
 * @param[in] ip: image primitives with 2D/3D lines
 * @param[in] w
 * @param[in] h
 * @param[in] opt_center: optimize the distortion center
 * @return final error
 */
double energy_minimization(lens_distortion_model &ldm, image_primitives &ip,
                           int w, int h, bool opt_center)
{
  int model_type = ldm.get_type();
  std::vector<bool> vtf(4);
  vtf[0] = vtf[1] = true;
  vtf[2] = vtf[3] = false;
  std::vector<bool> vtt(4,true);
  double error = 0.;
  double tf_error = 0.;
  double tt_error = 0.;
  lens_distortion_model tf_ldm, tt_ldm;
  if(!opt_center)
  {
    if(model_type == POLYNOMIAL)
    {
      error = model_center_estimation_2p_polynomial(ldm.get_distortion_center(), ip.get_lines(), ldm, w, h, vtf);
    }
    else
    {
      error = model_center_estimation_2p_quotient(ldm.get_distortion_center(), ip.get_lines(), ldm, w, h, vtf);
    }
  }
  else
  {
    if(model_type == POLYNOMIAL)
    {
      tf_error = model_center_estimation_2p_polynomial(ldm.get_distortion_center(), ip.get_lines(), ldm, w, h, vtf);
      tf_ldm = ldm;
      tt_error = model_center_estimation_2p_polynomial(ldm.get_distortion_center(), ip.get_lines(), ldm, w, h, vtt);
      tt_ldm = ldm;
    }
    else
    {
      tf_error = model_center_estimation_2p_quotient(ldm.get_distortion_center(), ip.get_lines(), ldm, w, h, vtf);
      tf_ldm = ldm;
      tt_error = model_center_estimation_2p_quotient(ldm.get_distortion_center(), ip.get_lines(), ldm, w, h, vtt);
      tt_ldm = ldm;
    }
    
    if((fabs(tf_ldm.get_distortion_center().x - tt_ldm.get_distortion_center().x)>0.2*w) ||
       (fabs(tf_ldm.get_distortion_center().y - tt_ldm.get_distortion_center().y)>0.2*h) ||
       (check_invertibility(tt_ldm, w, h) == false))
    {
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

//Iterative optimization
double iterative_optimization(
  const ami::subpixel_image_contours &contours,
  image_primitives &i_primitives,
  const float distance_point_line_max_hough,
  const int max_lines,
  const float angle_resolution,
  const float distance_resolution,
  const float distortion_parameter_resolution,
  const float angle_point_orientation_max_difference,
  const bool opt_center,
  const int width,
  const int height
)
{
  double final_error = 0.;
  
  //We initialize the previous model, the best one and the previous set of primitives
  i_primitives.get_distortion().get_d().resize(3);
  i_primitives.get_distortion().get_d()[2] = 0.;
  lens_distortion_model previous_model = i_primitives.get_distortion();
  lens_distortion_model best_model = previous_model;
  image_primitives previous_ip = i_primitives;
  //Tolerance for convergence
  double TOL = 1e-2;
  //Number of convergence iterations
  int convergence_iterations = 0;
  //Fails counter
  int fail_count = 0;
  //Number of points: current, best and next
  int num_points = count_points(i_primitives);
  int best_num_points = num_points;
  int next_num_points = num_points+TOL;
  lens_distortion_model ldm_tf;
  lens_distortion_model ldm_tt;
  //We apply the process until the number of points is not significantly greater
  //or until the process fails three times 
  while((next_num_points >= (num_points*(1+TOL))) || (fail_count<3))
  {
    double error = energy_minimization(previous_model, i_primitives, width, height, opt_center);
    if(check_invertibility(previous_model, width, height) == false)
    {
      if(fail_count==3)
        break;
    }
    
    i_primitives.clear();
    //CALL TO IMPROVED HOUGH WITH THE MODEL COMPUTED BEFORE
    line_equation_distortion_extraction_improved_hough(
      contours,
      i_primitives,
      distance_point_line_max_hough,
      max_lines,
      angle_resolution,
      distance_resolution,
      0., 
      0.,
      distortion_parameter_resolution,
      angle_point_orientation_max_difference, 
      true,
      previous_model
    );
    
    int local_num_points = count_points(i_primitives);
    if(local_num_points > next_num_points)
    {
      //We update the primitives only if the result is better
      if(local_num_points > best_num_points)
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
      if(fail_count == 3)
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
  return (final_error / count_points(i_primitives));
}

//------------------------------------------------------------------------------

// Load directory's files into a vector.
std::string get_filename(const std::string& filepath)
{
  const std::size_t found = filepath.rfind("/");
  if(found == std::string::npos)
    return filepath;
  
  return filepath.substr(found+1, filepath.size());
}

void split_filename(const std::string& filename, std::string& basename, std::string& extension)
{
  const std::size_t found = filename.rfind(".");
  if(found == std::string::npos)
  {
    basename = filename;
    extension = "";
    return;
  }

  basename = filename.substr(0, found);
  extension = filename.substr(found+1, filename.size());
}

void read_directory(std::vector<std::string>& out_filepaths, const std::string& inputFolder)
{
  DIR* rep = opendir(inputFolder.c_str());
  cout << "Reading the directory ..." << endl;

  cout << "DT_REG: " << (int)DT_REG << endl;
  cout << "DT_LNK: " << (int)DT_LNK << endl;
  while(struct dirent* ent = readdir(rep))
  {
    cout << "ent->d_type: " << (int)ent->d_type << " => " << ent->d_name << endl;
    // ignore if not a file or a link
    if(ent->d_type == DT_REG || ent->d_type == DT_LNK)
      out_filepaths.push_back(inputFolder + "/" + std::string(ent->d_name));
  }
  closedir(rep);
  cout << "... finished" << endl << endl;
  cout << "found " << out_filepaths.size() << " files." << endl;
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  if(argc != 11)
  {
    print_function_syntax_lens_distortion_correction_2p_iterative_optimization();
    exit(EXIT_FAILURE);
  }
  
  if(check_params_lens_distortion_correction_2p_iterative_optimization(argv) != 0)
  {
    manage_failure(argv,0);
    exit(EXIT_SUCCESS);
  }
  
  std::string input_folder(argv[1]);
  std::vector<std::string> input_files;

  read_directory(input_files, input_folder);
  
  if(input_files.empty())
    throw std::logic_error("No input image");

  for(const std::string& filepath: input_files)
    std::cout << "file: " << filepath << std::endl;
  
  //Load first image to read image size
  ami::image<unsigned char> inputReadImageSize(input_files[0]);
  int width = inputReadImageSize.width();
  int height = inputReadImageSize.height();
  inputReadImageSize.clear();

  int size_ = width*height; // image size
  float canny_high_threshold = atof(argv[3]); // high threshold for canny detector
  const float initial_distortion_parameter = atof(argv[4]); //left side of allowed distortion parameter interval
  const float final_distortion_parameter = atof(argv[5]); //Hough parameter
  const float distance_point_line_max_hough = atof(argv[6]); //Hough parameter
  //maximum difference allowed (in degrees) between edge and line orientation
  const float angle_point_orientation_max_difference = atof(argv[7]);
  const string tmodel(argv[8]);
  const string s_opt_c(argv[9]);
  const int max_lines = atof(argv[10]);
  const bool opt_center = (s_opt_c == string("True"));

  //const int max_lines = 100; //maximum number of lines estimated per image
  const float angle_resolution = 0.1; // angle discretization step (in degrees)
  const float distance_resolution = 1.; // line distance discretization step
  const float distortion_parameter_resolution = 0.1; //distortion parameter discretization step
  lens_distortion_model ini_ldm;
  if(tmodel == string("pol"))
    ini_ldm.set_type(POLYNOMIAL);
  else
    ini_ldm.set_type(DIVISION);

  int count_frames = 0;
  double best_score = 0;
  //object to store all the output edge lines structure and the best first estimation of the lens model distortion
  image_primitives all_primitives; 

  cout << "\nEstimation of the lens distortion model:" << endl;
  
  //We read all the input image of the directory
  for(const std::string& input_filepath: input_files)
  {
    if(all_primitives.get_lines().size() > max_lines)
      // We have enough data to compute the global distortion
      break;

    const std::string input_filename = get_filename(input_filepath);
    std::string input_basename, input_extension;
    split_filename(input_filename, input_basename, input_extension);
    cout << " - Estimation of the lens distortion model on the image " << count_frames << ": " << input_filename << endl;

    image_primitives i_primitives; //object to store output edge line structure
    ami::subpixel_image_contours contours;

    ami::image<unsigned char> input(input_filepath);
    ami::image<unsigned char> gray(width,height,1,0); //gray-level image to call canny

    //Converting the input image to gray level
    for(int i=0; i<size_; i++)
      gray[i] = 0.3 * input[i] + 0.59 * input[i+size_] + 0.11 * input[i+size_*2];
    input.clear();

    //ALGORITHM STAGE 1 : Detecting edges with Canny   
    ami::image<unsigned char> edges(width, height, 1, 0); //image to store edge information
    cout << "    Detecting edges with Canny" << endl;
    const float canny_low_threshold = 0.7; //default value for canny lower threshold
    contours = canny(gray,edges,canny_low_threshold, canny_high_threshold);
    edges.clear();
    gray.clear();

    //We clean the contours
    int neighborhood_radius = 2; //radius of neighborhood to take into account
    int min_neighbor_points = 4; //min number of contour points in a neighborhood
    double min_orientation_value = 0.95; //min average scalar product of neighborhood point orientation
    int min_distance_point = 1; //minimum distance between contour points
    contours.clean(
      neighborhood_radius,
      min_neighbor_points,
      min_orientation_value,
      min_distance_point
    );
    
    // Save Canny result as image
    {
      //We create writable 3 channel images for edges
      const vector<int>& index = contours.get_index();
      ami::image<unsigned char> edges3c(width, height, 3, 255);
      for(int i =0; i<(int)index.size(); i++)
      {
        edges3c[index[i]] = 0;
        edges3c[index[i]+size_] = 0;
        edges3c[index[i]+2*size_] = 0;
      }
      //Writing Canny detector output after the cleaning process
      edges3c.write( std::string(argv[2]) + input_basename + "_canny.png" );
      cout << "    edges detected" << endl;
    }
    
    //ALGORITHM STAGE 2 : Detecting lines with improved_hough_quotient  
    cout << "    Detecting lines with improved Hough and " <<  tmodel << " model for " << input_filename << endl;

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
    double image_error = iterative_optimization(
                            contours,
                            i_primitives,
                            distance_point_line_max_hough,
                            max_lines,
                            angle_resolution,
                            distance_resolution,
                            distortion_parameter_resolution,
                            angle_point_orientation_max_difference,
                            opt_center,
                            width,
                            height
                          );
    
    cout << "        lines detected: " << i_primitives.get_lines().size() << std::endl;
    cout << "        nb points: " << count_points(i_primitives) << std::endl;
    cout << "        image error: " << image_error << std::endl;
    
    //We check if the iterative optimization process finishes properly
    if(i_primitives.get_lines().size() == 0)
    {
      manage_failure(argv,0);
      continue;
    }
    ++count_frames; // one new frame successfully computed
    
    const int nbPoints = count_points(i_primitives);
    const double score = nbPoints / std::max(0.5, image_error);
    
    //Selection of the best estimation of the lens distortion model
    all_primitives.get_lines().insert(all_primitives.get_lines().end(), i_primitives.get_lines().begin(), i_primitives.get_lines().end());
    if(score > best_score)
    {
      all_primitives.get_distortion() = i_primitives.get_distortion();
      best_score = score;
    }
  }

  //ALGORITHM STAGE 4 : We minimize the error on the whole group of images
  lens_distortion_model previous_model = all_primitives.get_distortion();
  double global_error = energy_minimization(previous_model, all_primitives, width, height, opt_center);

  for(const std::string& input_filepath: input_files)
  {
    const std::string input_filename = get_filename(input_filepath);
    std::string input_basename, input_extension;
    split_filename(input_filename, input_basename, input_extension);
    
    //Drawing the detected lines on the original image to illustrate the results
    {
      ami::image<unsigned char> gray3c(width, height, 3, 255);
      drawHoughLines(all_primitives, gray3c);
      gray3c.write(std::string(argv[2]) + input_basename + "_hough.png");
    }

    //ALGORITHM STAGE 4 : Correcting the image distortion using the estimated model
    if(all_primitives.get_distortion().get_d().size() > 0)
    {
      cout << "    Correcting the distortion for " << input_filename << endl;
      ami::image<unsigned char> inputImage(input_filepath);
      
      if(all_primitives.get_distortion().get_type() == DIVISION)
      {
        ami::image<unsigned char> undistorted = undistort_quotient_image_inverse(
          inputImage, // input image
          all_primitives.get_distortion(), // lens distortion model
          3.0 // integer index to fix the way the corrected image is scaled to fit input size image
        );
        
        //Writing the distortion corrected image
        undistorted.write(std::string(argv[2]) + input_basename + "_undistort.png");
      }
      else
      {
        lens_distortion_model ldm = all_primitives.get_distortion();
        int vs = (ldm.get_d().size() == 2) ? 3 : 5;
        double *a = new double[vs];
        for(int i=0, ldmind = 0; i < vs; i++)
        {
          if(i%2 == 0)
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
          inputImage,
          vs-1,
          a,
          all_primitives.get_distortion().get_distortion_center(),
          2.0
        );
        
        delete []a;
        //Writing the distortion corrected image
        undistorted.write(std::string(argv[2]) + input_basename + "_undistort.png");
      }
    }
  }

  // WRITING OUTPUT TEXT DOCUMENTS
  // writing in a file the lens distortion model and the lines and associated points
  const std::string calib_file(std::string(argv[2]) + "global.calib");
  cout << "Export final result: " << calib_file << endl;
  all_primitives.write(calib_file);
  // writing function parameters and basic outputs :
  ofstream fs(std::string(argv[2]) + "output.log"); // Output file
  fs << "Selected parameters:" << endl;
  fs << "\t High Canny's threshold: " << argv[3] << endl;
  fs << "\t Initial normalized distortion parameter: " << argv[4] << endl;
  fs << "\t Final normalized distortion parameter: " << argv[5] << endl;
  fs << "\t Maximum distance between points and line: " << argv[6] << endl;
  fs << "\t Maximum difference between edge point and line orientations: " << argv[7]  << endl;
  fs << "\t Model applied: " << argv[8] << endl;
  fs << "\t Center optimization: " << argv[9] << endl;
  fs << "\t Maximum number of lines for the lens distortion model estimation: " << argv[10] << endl;
  fs << "-------------------------" << endl;
  fs << "Results: " << endl;
  fs << "\t Number of detected lines: " << all_primitives.get_lines().size() << endl;

  int count = count_points(all_primitives);

  fs << "\t Total amount of line points: " << count << endl;
  fs << "\t Distortion center: (" << all_primitives.get_distorsion_center().x <<
        ", " << all_primitives.get_distorsion_center().y << ")" << endl; 

  double p1 = 0.; 
  double p2 = 0.;
  bool is_division = (tmodel != std::string("pol"));
  compute_ps(p1, p2, all_primitives.get_distortion(), width, height, is_division);

  fs << "\t Estimated normalized distortion parameters: p1 = " << p1 << " p2 = " << p2 << endl;
  fs << "\t Average squared error distance in pixels between line and associated points = " << global_error << endl;
  fs.close();
  
  exit(EXIT_SUCCESS);
}
