% Automatic Lens Distortion Correction Using Two Parameter Polynomial and Division Models with iterative optimization

# ABOUT

* Author    : Miguel Alemán-Flores  <maleman@ctim.es>
              Luis Álvarez  <lalvarez@ctim.es>
              Luis Gómez <lgomez@ctim.es>
              Daniel Santana-Cedrés <dsantana@ctim.es>
* Copyright : (C) 2009-2014 IPOL Image Processing On Line http://www.ipol.im/
* License   : CC Creative Commons "Attribution-NonCommercial-ShareAlike" 
              see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en

# OVERVIEW

This source code provides an implementation of a lens distortion correction algorithm 
algorithm as described in IPOL

http://www.ipol.im/pub/algo/****

This program reads an input image and automatically estimates a 2 parameter 
polynomial or division lens distortion model, that it is used to correct image
distortion, after an iterative optimization process. 

The programs produces 4 outputs: 
   (1) The estimated distortion parameters. 
   (2) An output image with the results of Canny Edge Detector.
   (3) An output image with the estimated distorted lines.
   (4) An output image where the distortion model is applied to correct the 
       image distortion. 

# REQUIREMENTS

The code is written in ANSI C, and should compile on any system with
an ANSI C compiler.

The libpng header and libraries are required on the system for
compilation and execution. See http://www.libpng.org/pub/png/libpng.html

MacOSX already provides libpng in versions before Montain Lion.
However, if you have problems with png library, you can follow the next steps:
	- Install Hombrew (http://brew.hs) and run the commands:
		* brew doctor
		* brew update
		* brew install libpng
		* brew link libpng --force

# COMPILATION

We have checked our code on:
	- Fedora 11 (Leonidas) with GCC version 4.4.1-2
	- Ubuntu 14.04 LTS (trusty) with GCC version 4.8.2
	- Windows 7 with Cygwin 6.1 and GCC version 4.8.2
	- MacOSX 10.6.8 (Snow Leopard) Darwing Kernel version 10.8.0 wiht GCC version 4.2.1

Simply use the provided makefile, with the command `make`.
If you want to use the OpenMP library, please use `make OMP=1`.

Alternatively, you can manually compile
    g++ -Wall -Wextra -O3 lens_distortion_correction_2p_iterative_optimization.cpp 
		ami_primitives/subpixel_image_contours.cpp ami_lens_distortion/lens_distortion_procedures.cpp 
		ami_primitives/line_extraction.cpp ami_lens_distortion/lens_distortion_model.cpp 
		ami_primitives/line_points.cpp ami_image/io_png/io_png.cpp ami_lens_distortion/lens_distortion.cpp 
		ami_utilities/utilities.cpp ami_pol/ami_pol.cpp 
		-lpng -lm -o lens_distortion_correction_2p_iterative_optimization
		
		Or the following alternative compilation line for using the OpenMP library
		
		g++ -Wall -Wextra -O3 -fopenmp lens_distortion_correction_2p_iterative_optimization.cpp 
		ami_primitives/subpixel_image_contours.cpp ami_lens_distortion/lens_distortion_procedures.cpp 
		ami_primitives/line_extraction.cpp ami_lens_distortion/lens_distortion_model.cpp 
		ami_primitives/line_points.cpp ami_image/io_png/io_png.cpp ami_lens_distortion/lens_distortion.cpp 
		ami_utilities/utilities.cpp ami_pol/ami_pol.cpp 
		-lpng -lm -o lens_distortion_correction_2p_iterative_optimization

# USAGE

This program takes 10 parameters:
“exe_file  input_directory output_directory high_threshold_Canny initial_distortion_parameter final_distortion_parameter distance_point_line_max_hough angle_point_orientation_max_difference type_of_lens_distortion_model center_optimization primitives_file” 

* 'exe_file '                             : exe file (called ./lens_distortion_correction_2p_iterative_optimization) 
* 'input_directory'                       : input directory
* 'output_directory'                      : output directory
* 'high_threshold_Canny'                  : float value for the high threshold of the Canny method (between 0.7 and 1.0)
* 'initial_distortion_parameter'          : float value for the initial normalized distortion parameter (greater or equal to -0.5)
* 'final_distortion_parameter'            : float value for the final normalized distortion parameter (greater or equal to the initial value)
* 'distance_point_line_max_hough'         : maximum distance allowed between points and associated lines
* 'angle_point_orientation_max_difference': maximum difference (in degrees) of the point orientation angle and the line angle
* 'type_of_lens_distortion_model'         : type of the lens distortion model for the correction of the distortion (pol or div)
* 'center_optimization'                   : optimization of the center of the lens distortion model (True or False)
* 'max_lines'                            : int value for the maximal number of lines estimated per image (greater to 0)


# SOURCE CODE ORGANIZATION

The source code is organized in the following folders : 

* 'ami_pol' 		        : polynomial roots computation library.
* 'ami_filters' 		    : basic image filters including Gaussian convolution, gradient, 
					                Canny edge detector, etc.
* 'ami_image'			      : objects and basic methods to manage images.
* 'ami_image_draw'		  : basic procedure to draw primitives in an image.
* 'ami_lens_distortion'	: objects, methods and procedures to manage lens distortion models.
* 'ami_primitives'		  : basic objects and methods to manage primitives (lines and points)
* 'ami_utilities'		    : some auxilary macros and functions. 
* 'documentation'       : doxygen documentation of the source code.
* 'example'             : example input image and results.

# EXAMPLE

You can test the program with the provided test image (building.png) in the 
following way:

./lens_distortion_correction_2p_iterative_optimization /data/input/ /data/output/ 0.8 0.0 3.0 3.0 10.0 div True 1000

Furthermore, you can compare your results with the results present inside the folder 'example'

