#   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
#   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
#   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en


#Use OpenMP with make OMP=1
ifdef OMP
COPT	= -O3 -fopenmp
else
COPT	= -O3
endif

CFLAGS	= -Wall -Wextra $(COPT)
LDFLAGS	+= -lpng -lm

lens_distortion_correction_2p_iterative_optimization:\
						 lens_distortion_correction_2p_iterative_optimization.o\
						 subpixel_image_contours.o lens_distortion_procedures.o\
						 line_extraction.o lens_distortion_model.o line_points.o io_png.o\
						 lens_distortion.o ami_pol.o utilities.o
	$(CXX) $(CFLAGS) -o lens_distortion_correction_2p_iterative_optimization\
						lens_distortion_correction_2p_iterative_optimization.o\
						subpixel_image_contours.o lens_distortion_procedures.o\
						line_extraction.o lens_distortion_model.o\
						line_points.o io_png.o lens_distortion.o ami_pol.o\
						utilities.o $(LDFLAGS)
	
lens_distortion_correction_2p_iterative_optimization.o:\
	lens_distortion_correction_2p_iterative_optimization.cpp
	$(CXX) $(CFLAGS) -c lens_distortion_correction_2p_iterative_optimization.cpp

subpixel_image_contours.o: ami_primitives/subpixel_image_contours.cpp
	$(CXX) $(CFLAGS) -c ami_primitives/subpixel_image_contours.cpp

lens_distortion.o: ami_lens_distortion/lens_distortion.cpp
	$(CXX) $(CFLAGS) -c ami_lens_distortion/lens_distortion.cpp
	
lens_distortion_procedures.o: ami_lens_distortion/lens_distortion_procedures.cpp
	$(CXX) $(CFLAGS) -c ami_lens_distortion/lens_distortion_procedures.cpp

lens_distortion_model.o: ami_lens_distortion/lens_distortion_model.cpp
	$(CXX) $(CFLAGS) -c ami_lens_distortion/lens_distortion_model.cpp
	
line_extraction.o: ami_primitives/line_extraction.cpp
	$(CXX) $(CFLAGS) -c ami_primitives/line_extraction.cpp	
	
io_png.o: ami_image/io_png/io_png.cpp
	$(CXX) $(CFLAGS) -c ami_image/io_png/io_png.cpp

line_points.o: ami_primitives/line_points.cpp
	$(CXX) $(CFLAGS) -c ami_primitives/line_points.cpp
	
ami_pol.o: ami_pol/ami_pol.cpp
	$(CXX) $(CFLAGS) -c ami_pol/ami_pol.cpp
	
utilities.o: ami_utilities/utilities.cpp
	$(CXX) $(CFLAGS) -c ami_utilities/utilities.cpp

.PHONY: clean
clean:
	$(RM) -f *.o lens_distortion_correction_2p_iterative_optimization
	
