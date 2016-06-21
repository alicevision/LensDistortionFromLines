/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */

#ifndef AMI_DLL_H
  #define AMI_DLL_H
#endif

#ifndef LINE_EXTRACTION_H
#define LINE_EXTRACTION_H
#include "subpixel_image_contours.h"
#include "image_primitives.h"
#include <iostream>
#include <string>

double line_equation_distortion_extraction_improved_hough(
                  const subpixel_image_contours &subpixel_contours,
                  image_primitives &image_primitive,
                  const float distance_point_line_max,
                  const int nlineas = 100,
                  const float angle_resolution = 0.1,
                  const float distance_resolution = 1.,
                  const float initial_distortion_parameter=0.0,
                  const float final_distortion_parameter=1.0,
                  const float distortion_parameter_resolution=0.1,
                  const float angle_point_orientation_max_difference=2.0,
                  const bool lens_distortion_estimation=true,
                  const lens_distortion_model& ini_ldm = lens_distortion_model());

#endif
