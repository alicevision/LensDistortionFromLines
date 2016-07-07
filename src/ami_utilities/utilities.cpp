/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


#include "utilities.h"
#include "../ami_primitives/point2d.h"

//------------------------------------------------------------------------------

/**
 * \fn std::vector < std::vector <unsigned int> >
        ami::utilities::boundary_neighborhood_9n(const unsigned int width,
                                                 const unsigned int height)
 * \brief function to compute the 9 neighborhood index of the image boundary
 * \author Luis Alvarez
 * \return b[k][l] where for a boundary point k : b[k][0]=m (image index point),
          b[k][1]=m-width, b[k][2]=m+width, b[k][3]=m+1, b[k][4]=m-1,
          b[k][5]=m-width+1, b[k][6]=m+width+1, b[k][7]=m-width-1,
          b[k][8]=m+width-1 (if point is outside the image we take the nearest
                             one inside the image)
*/

std::vector < std::vector <unsigned int> >
  boundary_neighborhood_9n(const unsigned int width, const unsigned int height)
{
  if( width<3 || height<3) return(std::vector< std::vector <unsigned int> >());
  unsigned int size=2*(width+height)-4;
  std::vector< std::vector <unsigned int> > neigh_vector(size,
                                                  std::vector<unsigned int>(9));

  unsigned int m,cont=0;

  // FIRST CORNER
  neigh_vector[cont][0]=0; // pixel position index
  neigh_vector[cont][1]=0; // North pixel position index
  neigh_vector[cont][2]=width; // South pixel position index
  neigh_vector[cont][3]=1; // Est pixel position index
  neigh_vector[cont][4]=0; // West pixel position index
  neigh_vector[cont][5]=1; // North-Est pixel position index
  neigh_vector[cont][6]=width+1; // South-Est pixel position index
  neigh_vector[cont][7]=0; // North-West pixel position index
  neigh_vector[cont][8]=0; // South-West pixel position index
  cont++;

  // BOUNDARY HORIZONTAL LINE
  for(unsigned int k=1,k_end=width-1;k!=k_end;++k,++cont){
    neigh_vector[cont][0]=k; // pixel position index
    neigh_vector[cont][1]=k; // North pixel position index
    neigh_vector[cont][2]=k+width; // South pixel position index
    neigh_vector[cont][3]=k+1; // Est pixel position index
    neigh_vector[cont][4]=k-1; // West pixel position index
    neigh_vector[cont][5]=k+1; // North-Est pixel position index
    neigh_vector[cont][6]=k+width+1; // South-Est pixel position index
    neigh_vector[cont][7]=k-1; // North-West pixel position index
    neigh_vector[cont][8]=k+width-1; // South-West pixel position index
  }

  // SECOND CORNER
  m=width-1;
  neigh_vector[cont][0]=m; // pixel position index
  neigh_vector[cont][1]=m; // North pixel position index
  neigh_vector[cont][2]=m+width; // South pixel position index
  neigh_vector[cont][3]=m; // Est pixel position index
  neigh_vector[cont][4]=m-1; // West pixel position index
  neigh_vector[cont][5]=m; // North-Est pixel position index
  neigh_vector[cont][6]=m+width; // South-Est pixel position index
  neigh_vector[cont][7]=m-1; // North-West pixel position index
  neigh_vector[cont][8]=m+width-1; // South-West pixel position index
  cont++;

  // BOUNDARY VERTICAL LINE
  for(unsigned int k=1,k_end=height-1;k!=k_end;++k,++cont){
    m=(k+1)*width-1;
    neigh_vector[cont][0]=m; // pixel position index
    neigh_vector[cont][1]=m-width; // North pixel position index
    neigh_vector[cont][2]=m+width; // South pixel position index
    neigh_vector[cont][3]=m; // Est pixel position index
    neigh_vector[cont][4]=m-1; // West pixel position index
    neigh_vector[cont][5]=m-width; // North-Est pixel position index
    neigh_vector[cont][6]=m+width; // South-Est pixel position index
    neigh_vector[cont][7]=m-width-1; // North-West pixel position index
    neigh_vector[cont][8]=m+width-1; // South-West pixel position index

  }

  // CORNER POINT
  m=height*width-1;
  neigh_vector[cont][0]=m; // pixel position index
  neigh_vector[cont][1]=m-width; // North pixel position index
  neigh_vector[cont][2]=m; // South pixel position index
  neigh_vector[cont][3]=m; // Est pixel position index
  neigh_vector[cont][4]=m-1; // West pixel position index
  neigh_vector[cont][5]=m-width; // North-Est pixel position index
  neigh_vector[cont][6]=m; // South-Est pixel position index
  neigh_vector[cont][7]=m-width-1; // North-West pixel position index
  neigh_vector[cont][8]=m-1; // South-West pixel position index
  cont++;

  // HORIZONTAL LINE
  for(unsigned int k=2,k_end=width;k!=k_end;++k,++cont){
    m=height*width-k;
    neigh_vector[cont][0]=m; // pixel position index
    neigh_vector[cont][1]=m-width; // North pixel position index
    neigh_vector[cont][2]=m; // South pixel position index
    neigh_vector[cont][3]=m+1; // Est pixel position index
    neigh_vector[cont][4]=m-1; // West pixel position index
    neigh_vector[cont][5]=m-width+1; // North-Est pixel position index
    neigh_vector[cont][6]=m+1; // South-Est pixel position index
    neigh_vector[cont][7]=m-width-1; // North-West pixel position index
    neigh_vector[cont][8]=m-1; // South-West pixel position index
  }

  // CORNER POINT
  m=(height-1)*width;
  neigh_vector[cont][0]=m; // pixel position index
  neigh_vector[cont][1]=m-width; // North pixel position index
  neigh_vector[cont][2]=m; // South pixel position index
  neigh_vector[cont][3]=m+1; // Est pixel position index
  neigh_vector[cont][4]=m; // West pixel position index
  neigh_vector[cont][5]=m-width+1; // North-Est pixel position index
  neigh_vector[cont][6]=m+1; // South-Est pixel position index
  neigh_vector[cont][7]=m-width; // North-West pixel position index
  neigh_vector[cont][8]=m; // South-West pixel position index
  cont++;

  // VERTICAL LINE
  for(unsigned int k=2,k_end=height;k!=k_end;++k,++cont){
    m=(height-k)*width;
    neigh_vector[cont][0]=m; // pixel position index
    neigh_vector[cont][1]=m-width; // North pixel position index
    neigh_vector[cont][2]=m+width; // South pixel position index
    neigh_vector[cont][3]=m+1; // Est pixel position index
    neigh_vector[cont][4]=m; // West pixel position index
    neigh_vector[cont][5]=m-width+1; // North-Est pixel position index
    neigh_vector[cont][6]=m+width+1; // South-Est pixel position index
    neigh_vector[cont][7]=m-width; // North-West pixel position index
    neigh_vector[cont][8]=m+width; // South-West pixel position index
  }

  return(neigh_vector);
}

//------------------------------------------------------------------------------

/**
 * \fn std::vector < std::vector <unsigned int> >
        ami::utilities::boundary_neighborhood_5n(const unsigned int width,
                                                 const unsigned int height)
 * \brief function to compute the 5 neighborhood index of the image boundary
 * \author Luis Alvarez
 * \return b[k][l] where for a boundary point k : b[k][0]=m (image index point),
            b[k][1]=m-width, b[k][2]=m+width, b[k][3]=m+1, b[k][4]=m-1,
            (if point is outside the image we take the nearest one inside the
             image)
*/
std::vector < std::vector <unsigned int> >
  boundary_neighborhood_5n(const unsigned int width, const unsigned int height)
{
  if( width<3 || height<3) return(std::vector< std::vector <unsigned int> >());
  unsigned int size=2*(width+height)-4;
  std::vector< std::vector <unsigned int> > neigh_vector(size,
                                                  std::vector<unsigned int>(5));

  unsigned int m,cont=0;

  // FIRST CORNER
  neigh_vector[cont][0]=0; // pixel position index
  neigh_vector[cont][1]=0; // North pixel position index
  neigh_vector[cont][2]=width; // South pixel position index
  neigh_vector[cont][3]=1; // Est pixel position index
  neigh_vector[cont][4]=0; // West pixel position index
  cont++;

  // BOUNDARY HORIZONTAL LINE
  for(unsigned int k=1,k_end=width-1;k!=k_end;++k,++cont){
    neigh_vector[cont][0]=k; // pixel position index
    neigh_vector[cont][1]=k; // North pixel position index
    neigh_vector[cont][2]=k+width; // South pixel position index
    neigh_vector[cont][3]=k+1; // Est pixel position index
    neigh_vector[cont][4]=k-1; // West pixel position index
  }

  // SECOND CORNER
  neigh_vector[cont][0]=width-1; // pixel position index
  neigh_vector[cont][1]=width-1; // North pixel position index
  neigh_vector[cont][2]=2*width-1; // South pixel position index
  neigh_vector[cont][3]=width-1; // Est pixel position index
  neigh_vector[cont][4]=width-2; // West pixel position index
  cont++;

  // BOUNDARY VERTICAL LINE
  for(unsigned int k=1,k_end=height-1;k!=k_end;++k,++cont){
    m=(k+1)*width-1;
    neigh_vector[cont][0]=m; // pixel position index
    neigh_vector[cont][1]=m-width; // North pixel position index
    neigh_vector[cont][2]=m+width; // South pixel position index
    neigh_vector[cont][3]=m; // Est pixel position index
    neigh_vector[cont][4]=m-1; // West pixel position index
  }

  // CORNER POINT
  m=height*width-1;
  neigh_vector[cont][0]=m; // pixel position index
  neigh_vector[cont][1]=m-width; // North pixel position index
  neigh_vector[cont][2]=m; // South pixel position index
  neigh_vector[cont][3]=m; // Est pixel position index
  neigh_vector[cont][4]=m-1; // West pixel position index
  cont++;

  // HORIZONTAL LINE
  for(unsigned int k=2,k_end=width;k!=k_end;++k,++cont){
    m=height*width-k;
    neigh_vector[cont][0]=m; // pixel position index
    neigh_vector[cont][1]=m-width; // North pixel position index
    neigh_vector[cont][2]=m; // South pixel position index
    neigh_vector[cont][3]=m+1; // Est pixel position index
    neigh_vector[cont][4]=m-1; // West pixel position index
  }

  // CORNER POINT
  m=(height-1)*width;
  neigh_vector[cont][0]=m; // pixel position index
  neigh_vector[cont][1]=m-width; // North pixel position index
  neigh_vector[cont][2]=m; // South pixel position index
  neigh_vector[cont][3]=m+1; // Est pixel position index
  neigh_vector[cont][4]=m; // West pixel position index
  cont++;

  // VERTICAL LINE
  for(unsigned int k=2,k_end=height;k!=k_end;++k,++cont){
    m=(height-k)*width;
    neigh_vector[cont][0]=m; // pixel position index
    neigh_vector[cont][1]=m-width; // North pixel position index
    neigh_vector[cont][2]=m+width; // South pixel position index
    neigh_vector[cont][3]=m+1; // Est pixel position index
    neigh_vector[cont][4]=m; // West pixel position index
  }

  return(neigh_vector);
}

//------------------------------------------------------------------------------

void drawHoughLines(image_primitives ip, ami::image<unsigned char> &bn)
{
	//DRAW THE LINES IN THE IMAGE
	ami::image_draw imgd;

	for(int i = 0; i<(int)ip.get_lines().size(); i++)
  {
    double a,b,c;
    ip.get_lines()[i].get_abc(a,b,c);
    int ia=round(a),ic=round(c);
    float point_size=3.;
    /* SEED TO INITIALIZE RANDOM FUNCTION */
    unsigned int graine =(unsigned int) (255*(fabs(ia)+1e-2)+ic);
    srand(graine);
    rand();
    unsigned char red= 64+(unsigned char) (172.*((double) rand()/RAND_MAX));
    unsigned char green=(unsigned char) (255.*((double) rand()/RAND_MAX));
    unsigned char blue= (unsigned char) (172.*((double) rand()/RAND_MAX));
    for(int j = 0; j<(int)ip.get_lines()[i].get_points().size(); j++)
    {
        point2d<double> ori = ip.get_lines()[i].get_points()[j];
        //DRAW THE ORIGINAL POINT IN r,g,b
        imgd.draw_cercle(bn,ori.x,ori.y,point_size,red,green,blue);
    }
  }
}

//------------------------------------------------------------------------------

void invert(ami::image<unsigned char> &input,
            ami::image<unsigned char> &output)
{
	int tam = input.width()*input.height()*input.nChannels();
	for(int i=0; i<tam; i++) output[i] = 255 - input[i];
}

//------------------------------------------------------------------------------

void print_function_syntax_lens_distortion_correction_2p_iterative_optimization()
{
	cout << "FUNCTION SYNTAX:" << endl;
	cout << "exe_file  /data/input/  /data/output/";
	cout << "initial_distortion_parameter_hough ";
	cout << "final_distortion_parameter_hough ";
	cout << "distance_point_line_max_hough ";
	cout << "angle_point_orientation_max_difference ";
	cout << "type_of_lens_distortion_model ";
	cout << "center_optimization ";
  cout << "max_lines ";
	cout << "primitives_file" << endl;
	
	cout << "Parameters description:" << endl;
	cout << "  exe_file: executable name, by default: ";
	cout << "lens_distortion_correction_2p_iterative_optimization"
			 << endl;
	cout << "  /data/input/: input directory" << endl;
	cout << "  /data/output/: output directory"<< endl;
	cout << "  high_threshold_Canny: float value for the high threshold of the ";
	cout << "Canny method (between 0.7 and 1.0)" << endl;
	cout << "  initial_distortion_parameter: float value for the initial ";
	cout << "normalized distortion parameter (greater or equal to -0.5)" << endl;
	cout << "  final_distortion_parameter: float value for the final normalized ";
	cout << "distortion parameter (greater or equal to the initial value)" << endl;
	cout << "  distance_point_line_max_hough: maximum distance allowed between ";
	cout << "points and associated lines." << endl;
	cout << "  angle_point_orientation_max_difference: maximum difference ";
	cout << "(in degrees) of the point orientation angle and the line angle" << endl;
	cout << "  type_of_lens_distortion_model: type of the lens distortion model for ";
	cout << "the correction of the distortion (pol or div)" << endl;
	cout << "  center_optimization: optimization of the center of the lens ";
	cout << "distortion model (True or False)" << endl;
  cout << "  max_lines: maximum number of lines estimated per image (greater than 0)" << endl << endl;
	cout << "Example command:" << endl;
	cout << " ./lens_distortion_correction_2p_iterative_optimization /data/input/";
	cout << " /data/output/ 0.8 ";
	cout << "0.0 3.0 3.0 10.0 div True 20" << endl;
}

//------------------------------------------------------------------------------

int check_params_lens_distortion_correction_2p_iterative_optimization(char *argv[])
{
	//Messages
	string error_message_head("The argument ");
	string error_message_image(" must be a png image.\n");
	string error_message_numgt(" must be greater than ");
	string error_message_numgte(" must be greater or equal to ");
	string error_message("");

//	//Input image
//	string str(argv[1]);
//	int pos = str.find_last_of('.');
//	string ext = str.substr(pos+1);
//	if(ext != "png" && ext != "PNG")
//		error_message += error_message_head + "1 (input image)" +
//										 error_message_image;

//	//Canny image
//	str = string(argv[2]);
//	pos = str.find_last_of('.');
//	ext = str.substr(pos+1);
//	if(ext != "png" && ext != "PNG")
//		error_message += error_message_head + "2 (output Canny image)" +
//										 error_message_image;
//
//	//Hough image
//	str = string(argv[3]);
//	pos = str.find_last_of('.');
//	ext = str.substr(pos+1);
//	if(ext != "png" && ext != "PNG")
//		error_message += error_message_head + "3 (output Hough lines image)" +
//										 error_message_image;
//
//	//Corrected image
//	str = string(argv[4]);
//	pos = str.find_last_of('.');
//	ext = str.substr(pos+1);
//	if(ext != "png" && ext != "PNG")
//		error_message += error_message_head + "4 (output corrected image)" +
//										 error_message_image;

	//High threshold Canny
	float num = atof(argv[3]);
	if(num <= 0.7)
		error_message += error_message_head + "3 (high threshold Canny)" +
										 error_message_numgt + "0.7\n";

	//Initial distortion parameter
	num = atof(argv[4]);
	if(num<-0.5)
		error_message += error_message_head + "4 (initial distortion parameter)" +
										 error_message_numgte + "-0.5\n";

	//Final distortion parameter
	num = atof(argv[5]);
	if(num < atof(argv[4]))
		error_message += error_message_head + "5 (final distortion parameter)" +
										 error_message_numgte + argv[4]+"\n";

	//Maximum distance between points and line
	num = atof(argv[6]);
	if(num < 0.0)
		error_message += error_message_head +
										 "6 (maximum distance between points and line)" +
										 error_message_numgte + "0.0\n";

	//Maximum difference between point angle and line angle
	num = atof(argv[7]);
	if(num < 0.0 || num > 45.0)
		error_message += error_message_head + "7 (maximum difference between point"+
									" angle and line angle) must be between 0.0 and 45.0\n";

	//Type of the lens distortion model
	string val(argv[8]);
	if(val!=string("pol") && val!=string("div"))
		error_message += error_message_head + "8 (type of the lens distortion model)."+
									" The allowed values are: pol or div\n";
									
	//Optimization of the center of the lens distortion model
	val.clear();
	val = string(argv[9]);
	if(val!=string("True") && val!=string("False"))
		error_message += error_message_head + "9 (optimization of the center of "+
									"the lens distortion model). The allowed values are: True or "+
									"False\n";
  
  //Maximum number of frames to export
	num = atof(argv[10]);
	if(num < 0)
		error_message += error_message_head +
										 "10 (maximum number of lines estimated per image)" +
										 error_message_numgt + " 0\n";
  
	if(error_message.length() > 0)
	{
		cout << error_message;
		return -1;
	}
	return 0;
}

//------------------------------------------------------------------------------

int count_points(image_primitives ip)
{
	int count = 0;
	for(int i=0; i<(int)ip.get_lines().size(); i++)
		count += ip.get_lines()[i].get_points().size();
	return count;
}

//------------------------------------------------------------------------------

void manage_failure(char *argv[], int code)
{
	ami::image<unsigned char> input(argv[1]);
	//Write the output images as a copy of the input
	input.write(argv[2]);
	//Write output file
	ofstream fs("output.txt");
	fs << "Selected parameters:" << endl;
    fs << "\t High Canny's threshold: " << argv[3] << endl;
    fs << "\t Initial normalized distortion parameter: " << argv[4] << endl;
    fs << "\t Final normalized distortion parameter: " << argv[5] << endl;
    fs << "\t Maximum distance between points and line: " << argv[6] << endl;
    fs << "\t Maximum different between edge point and line orientations: "
     << argv[7]	<< endl;
		fs << "\t Model applied: " << argv[8] << endl;
		fs << "\t Center optimization: " << argv[9] << endl;
    fs << "\t Maximum frames for estimate the lens distortion model: " << argv[10] << endl;
    fs << "-------------------------" << endl;
    fs << "Results: " << endl;
	fs << "\t Program failed:" << endl;
	switch(code)
	{
		case 0:
			fs << "\t\t Probably due to a bad parameter choice" << endl;
			break;
		case 1:
			fs << "\t\t The lens distortion model doesn't pass the invertivility test" 
			   << endl;
			break;
	}
	fs.close();
	//And write an empty output file of primitives information
	ofstream pf(argv[8]);
	pf.close();
}

//------------------------------------------------------------------------------

double update_rsqmax(point2d<double> dc, int w, int h)
{
    point2d<double> corner(0,0);
    double max_distance_corner2= (dc-corner).norm2();
    corner.y=h;
    double distance_corner2=(dc-corner).norm2();
    if(distance_corner2>max_distance_corner2)
        max_distance_corner2=distance_corner2;
    corner.x=w;
    distance_corner2=(dc-corner).norm2();
    if(distance_corner2>max_distance_corner2)
        max_distance_corner2=distance_corner2;
    corner.y=0;
    distance_corner2=(dc-corner).norm2();
    if(distance_corner2>max_distance_corner2)
        max_distance_corner2=distance_corner2;

    //Return the maximum r^2_max
    return max_distance_corner2;
}

//------------------------------------------------------------------------------

void compute_ps(double &p1, double &p2, const lens_distortion_model &ldm,
								int w, int h, bool quo)
{
	double dmi = update_rsqmax(ldm.get_distortion_center(),w,h);
	double r1  = sqrt(dmi);
	double r2  = r1/2;
	double k1,k2;
	if(ldm.get_d().size()<3)
  {
    k1 = ldm.get_d()[1];
    k2 = 0.;
  }
  else
  {
    k1 = ldm.get_d()[1];
    k2 = ldm.get_d()[2];
  }
	if(quo)
	{
		p1 = (1/(1+k1*r1*r1+k2*r1*r1*r1*r1))-1;
		p2 = (1/(1+k1*r2*r2+k2*r2*r2*r2*r2))-1;
	}
	else
	{
		p1 = k1*4*r2*r2 + k2*16*r2*r2*r2*r2;
		p2 = k1*r2*r2 + k2*r2*r2*r2*r2;
	}
}

//------------------------------------------------------------------------------

int ami2_gauss(double **A,double *b,int N)
/* RETURNS THE SOLUTION IN THE VECTOR b */
     /*double **A,*b; int N; */
{
  double **PASO,max,paso,mul;
  int i,j,i_max,k;
  PASO=(double **)malloc(sizeof(double*)*N);
  for(i=0;i<N;i++)
    PASO[i]=(double *)malloc(sizeof(double)*(N+1));

  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      PASO[i][j]=A[i][j];
    }
    PASO[i][N]=b[i];
  }

  for(i=0;i<N;i++){
    max=fabs(PASO[i][i]);
    i_max=i;
    for(j=i;j<N;j++){
       if(fabs(PASO[j][i])>max){
         i_max=j; max=fabs(PASO[j][i]);
       }
    }
    if(max<10e-30){
      printf("System does not have solution 0\n");
      for(i=0;i<N;i++)
        free(PASO[i]);
      free(PASO);
      return(-1);
    }
    if(i_max>i){
      for(k=0;k<=N;k++){
        paso=PASO[i][k];
        PASO[i][k]=PASO[i_max][k];
        PASO[i_max][k]=paso;
      }
    }
    for(j=i+1;j<N;j++){
      mul=-PASO[j][i]/PASO[i][i];
      for(k=i;k<=N;k++) PASO[j][k]+=mul*PASO[i][k];
    }
  }
  if(fabs(PASO[N-1][N-1])<10e-30){
      printf("System does not have solution 1\n");
      for(i=0;i<N;i++)
       free(PASO[i]);
      free(PASO);
      return(-1);
    }

  for(i=N-1;i>0;i--){
    for(j=i-1;j>=0;j--){
      mul=-PASO[j][i]/PASO[i][i];
      for(k=i;k<=N;k++) PASO[j][k]+=mul*PASO[i][k];
    }
  }
  for(i=0;i<N;i++)
      b[i]=PASO[i][N]/PASO[i][i];

  for(i=0;i<N;i++)
    free(PASO[i]);
  free(PASO);
  return(0);
}

//------------------------------------------------------------------------------

bool check_invertibility(lens_distortion_model &ldm, int w, int h)
{
	//double r1 = update_rsqmax(ldm.get_distortion_center(), w, h);
	double r1sq = update_rsqmax(ldm.get_distortion_center(), w, h);//r1*r1;
	double r1p4 = r1sq*r1sq;
	double k1 = ldm.get_d()[1];
	double k2 = ldm.get_d()[2];
	
	if(ldm.get_type()==POLYNOMIAL)
	{
		if( (((r1sq*k1) < (-2./3.)) && ((9.*r1p4*k1*k1 - 20.*r1p4*k2) < 0.0))
				||
				(((r1sq*k1) >=(-2./3.)) && ((5.*r1p4*k2 + 3.*r1sq*k1 + 1.) > 0.0))
		  )
			return true;
	}
	else
	{
    if((r1sq*k1>-2.) && (r1sq*k1<2.))
    {
      if((r1p4*k2>-1.-r1sq*k1) && (r1p4*k2<((1.-r1sq*k1)/3.)))
        return true;
    }
    else
    {
      if(r1sq*k1>=2.)
      {
        if((r1p4*k2>-1.-r1sq*k1) && (r1p4*k2<(-r1p4*k2/12.)))
          return true;
      }
    }
		/*if( (((r1sq*k1) < 2.) && ((r1sq*k1 + 3.*r1p4*k2 - 1.) < 0.0))
				||
				(((r1sq*k1) >=2.) && ((r1p4*k1*k1 + 12.*r1p4*k2) < 0.0))
			)
			return true;*/
	}
	return false;
}
