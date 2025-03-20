   typedef long long fixed;
   #define fixeddot 16

   #define VERBOSE 0
   #define BOOSTBLURFACTOR 90.0

   #include <stdio.h>
   #include <stdlib.h>
   #include <math.h>

   int read_pgm_image(char *infilename, unsigned char **image, int *rows,
      int *cols);
   int write_pgm_image(char *outfilename, unsigned char *image, int rows,
      int cols, char *comment, int maxval);

   void canny(unsigned char *image, int rows, int cols, float sigma,
            float tlow, float thigh, unsigned char **edge, char *fname);
   void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma,
         short int **smoothedim);
   void make_gaussian_kernel(float sigma, float **kernel, int *windowsize);
   void derrivative_x_y(short int *smoothedim, int rows, int cols,
         short int **delta_x, short int **delta_y);
   void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
         short int **magnitude);
   void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols,
         float tlow, float thigh, unsigned char *edge);
   void radian_direction(short int *delta_x, short int *delta_y, int rows,
      int cols, float **dir_radians, int xdirtag, int ydirtag);
   double angle_radians(double x, double y);

   int main(int argc, char *argv[])
   {
      char *infilename = NULL;  /* Name of the input image */
      char *dirfilename = NULL; /* Name of the output gradient direction image */
      char outfilename[128];    /* Name of the output "edge" image */
      char composedfname[128];  /* Name of the output "direction" image */
      unsigned char *image;     /* The input image */
      unsigned char *edge;      /* The output edge image */
      int rows, cols;           /* The dimensions of the image. */
      float sigma,              /* Standard deviation of the gaussian kernel. */
      tlow,               /* Fraction of the high threshold in hysteresis. */
      thigh;             
      if(argc < 5){
      fprintf(stderr,"\n<USAGE> %s image sigma tlow thigh [writedirim]\n",argv[0]);
         fprintf(stderr,"\n      image:      An image to process. Must be in ");
         fprintf(stderr,"PGM format.\n");
         fprintf(stderr,"      sigma:      Standard deviation of the gaussian");
         fprintf(stderr," blur kernel.\n");
         fprintf(stderr,"      tlow:       Fraction (0.0-1.0) of the high ");
         fprintf(stderr,"edge strength threshold.\n");
         fprintf(stderr,"      thigh:      Fraction (0.0-1.0) of the distribution");
         fprintf(stderr," of non-zero edge\n                  strengths for ");
         fprintf(stderr,"hysteresis. The fraction is used to compute\n");
         fprintf(stderr,"                  the high edge strength threshold.\n");
         fprintf(stderr,"      writedirim: Optional argument to output ");
         fprintf(stderr,"a floating point");
         fprintf(stderr," direction image.\n\n");
         exit(1);
      }

      infilename = argv[1];
      sigma = atof(argv[2]);
      tlow = atof(argv[3]);
      thigh = atof(argv[4]);

      if(argc == 6) dirfilename = infilename;
      else dirfilename = NULL;

      if(VERBOSE) printf("Reading the image %s.\n", infilename);
      if(read_pgm_image(infilename, &image, &rows, &cols) == 0){
         fprintf(stderr, "Error reading the input image, %s.\n", infilename);
         exit(1);
      }

      if(VERBOSE) printf("Starting Canny edge detection.\n");
      if(dirfilename != NULL){
         sprintf(composedfname, "%s_s_%3.2f_l_%3.2f_h_%3.2f.fim", infilename,
         sigma, tlow, thigh);
         dirfilename = composedfname;
      }

      ///////  
      canny(image, rows, cols, sigma, tlow, thigh, &edge, dirfilename);
      ///////
      sprintf(outfilename, "%s_s_%3.2f_l_%3.2f_h_%3.2f.pgm", infilename,
         sigma, tlow, thigh);
      if(VERBOSE) printf("Writing the edge iname in the file %s.\n", outfilename);
      if(write_pgm_image(outfilename, edge, rows, cols, "", 255) == 0){
         fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
         exit(1);
      }
      return 0;
   }

   void canny(unsigned char *image, int rows, int cols, float sigma,
            float tlow, float thigh, unsigned char **edge, char *fname)
   {
      FILE *fpdir=NULL;          /* File to write the gradient image to.     */
      unsigned char *nms;        /* Points that are local maximal magnitude. */
      short int *smoothedim,     /* The image after gaussian smoothing.      */
               *delta_x,        /* The first devivative image, x-direction. */
               *delta_y,        /* The first derivative image, y-direction. */
               *magnitude;      /* The magnitude of the gadient image.      */
      int r, c, pos;
      float *dir_radians=NULL;   /* Gradient direction image.                */

      if(VERBOSE) printf("Smoothing the image using a gaussian kernel.\n");
      gaussian_smooth(image, rows, cols, sigma, &smoothedim);

      if(VERBOSE) printf("Computing the X and Y first derivatives.\n");
      derrivative_x_y(smoothedim, rows, cols, &delta_x, &delta_y);

      if(fname != NULL)
      {
         radian_direction(delta_x, delta_y, rows, cols, &dir_radians, -1, -1);
         
         if((fpdir = fopen(fname, "wb")) == NULL)
      {
            fprintf(stderr, "Error opening the file %s for writing.\n", fname);
            exit(1);
         }
         fwrite(dir_radians, sizeof(float), rows*cols, fpdir);
         fclose(fpdir);
         free(dir_radians);
      }
      
      if(VERBOSE) printf("Computing the magnitude of the gradient.\n");
      magnitude_x_y(delta_x, delta_y, rows, cols, &magnitude);

      
      if(VERBOSE) printf("Doing the non-maximal suppression.\n");
      if((nms = (unsigned char *) malloc(rows*cols*sizeof(unsigned char)))==NULL)
      {
         fprintf(stderr, "Error allocating the nms image.\n");
         exit(1);
      }
      non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);

      if(VERBOSE) printf("Doing hysteresis thresholding.\n");
      if((*edge=(unsigned char *)malloc(rows*cols*sizeof(unsigned char))) ==NULL)
      {
         fprintf(stderr, "Error allocating the edge image.\n");
         exit(1);
      }
      apply_hysteresis(magnitude, nms, rows, cols, tlow, thigh, *edge);

      
      free(smoothedim);
      free(delta_x);
      free(delta_y);
      free(magnitude);
      free(nms);
   }


   void radian_direction(short int *delta_x, short int *delta_y, int rows,
      int cols, float **dir_radians, int xdirtag, int ydirtag)
   {
      int r, c, pos;
      float *dirim=NULL;
      double dx, dy;
      
      if((dirim = (float *) malloc(rows*cols* sizeof(float))) == NULL){
         fprintf(stderr, "Error allocating the gradient direction image.\n");
         exit(1);
      }
      *dir_radians = dirim;

      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++){
            dx = (double)delta_x[pos];
            dy = (double)delta_y[pos];

            if(xdirtag == 1) dx = -dx;
            if(ydirtag == -1) dy = -dy;

            dirim[pos] = (float)angle_radians(dx, dy);
         }
      }
   }

   double angle_radians(double x, double y)
   {
      double xu, yu, ang;

      xu = fabs(x);
      yu = fabs(y);

      if((xu == 0) && (yu == 0)) return(0);

      ang = atan(yu/xu);

      if(x >= 0){
         if(y >= 0) return(ang);
         else return(2*M_PI - ang);
      }
      else{
         if(y >= 0) return(M_PI - ang);
         else return(M_PI + ang);
      }
   }

   void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
         short int **magnitude)
   {
      int r, c, pos, sq1, sq2;

      if((*magnitude = (short *) malloc(rows*cols* sizeof(short))) == NULL){
         fprintf(stderr, "Error allocating the magnitude image.\n");
         exit(1);
      }

      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++){
            sq1 = (int)delta_x[pos] * (int)delta_x[pos];
            sq2 = (int)delta_y[pos] * (int)delta_y[pos];
            (*magnitude)[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
         }
      }

   }

   void derrivative_x_y(short int *smoothedim, int rows, int cols,
         short int **delta_x, short int **delta_y)
   {
      int r, c, pos;
      
      if(((*delta_x) = (short *) malloc(rows*cols* sizeof(short))) == NULL){
         fprintf(stderr, "Error allocating the delta_x image.\n");
         exit(1);
      }
      if(((*delta_y) = (short *) malloc(rows*cols* sizeof(short))) == NULL){
         fprintf(stderr, "Error allocating the delta_x image.\n");
         exit(1);
      }


      for(r=0;r<rows;r++){
         pos = r * cols;
         (*delta_x)[pos] = smoothedim[pos+1] - smoothedim[pos];
         pos++;
         for(c=1;c<(cols-1);c++,pos++){
            (*delta_x)[pos] = smoothedim[pos+1] - smoothedim[pos-1];
         }
         (*delta_x)[pos] = smoothedim[pos] - smoothedim[pos-1];
      }

      for(c=0;c<cols;c++){
         pos = c;
         (*delta_y)[pos] = smoothedim[pos+cols] - smoothedim[pos];
         pos += cols;
         for(r=1;r<(rows-1);r++,pos+=cols){
            (*delta_y)[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
         }
         (*delta_y)[pos] = smoothedim[pos] - smoothedim[pos-cols];
      }
   }


   void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma,
         short int **smoothedim)
   {
      int r, c, rr, cc,     /* Counter variables. */
         windowsize,        /* Dimension of the gaussian kernel. */
         center;            /* Half of the windowsize. */
      float *tempim,        /* Buffer for separable filter gaussian smoothing. */
            *kernel,        /* A one dimensional gaussian kernel. */
            dot,            /* Dot product summing variable. */
            sum;            /* Sum of the kernel weights variable. */

      if(VERBOSE) printf("   Computing the gaussian smoothing kernel.\n");
      make_gaussian_kernel(sigma, &kernel, &windowsize);
      center = windowsize / 2;

      if((tempim = (float *) malloc(rows*cols* sizeof(float))) == NULL)
      {
         fprintf(stderr, "Error allocating the buffer image.\n");
         exit(1);
      }
      if(((*smoothedim) = (short int *) malloc(rows*cols*sizeof(short int))) == NULL)
      {
         fprintf(stderr, "Error allocating the smoothed image.\n");
         exit(1);
      }

      
      if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
      for(r=0;r<rows;r++)
      {
         for(c=0;c<cols;c++)
      {
            dot = 0.0;
            sum = 0.0;
            for(cc=(-center);cc<=center;cc++)
         {
               if(((c+cc) >= 0) && ((c+cc) < cols))
            {
                  dot += (float)image[r*cols+(c+cc)] * kernel[center+cc];
                  sum += kernel[center+cc];
               }
            }
            tempim[r*cols+c] = dot/sum;
         }
      }

      
      if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
      for(c=0;c<cols;c++)
      {
         for(r=0;r<rows;r++)
      {
            sum = 0.0;
            dot = 0.0;
            for(rr=(-center);rr<=center;rr++)
         {
               if(((r+rr) >= 0) && ((r+rr) < rows))
            {
                  dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
                  sum += kernel[center+rr];
               }
            }
            (*smoothedim)[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
         }
      }

      free(tempim);
      free(kernel);
   }

   void make_gaussian_kernel(float sigma, float **kernel, int *windowsize)
   {
      int i, center;
      float x, fx, sum=0.0;

      *windowsize = 1 + 2 * ceil(2.5 * sigma);
      center = (*windowsize) / 2;

      if(VERBOSE) printf("      The kernel has %d elements.\n", *windowsize);
      if((*kernel = (float *) malloc((*windowsize)* sizeof(float))) == NULL){
         fprintf(stderr, "Error callocing the gaussian kernel array.\n");
         exit(1);
      }

      for(i=0;i<(*windowsize);i++){
         x = (float)(i - center);
         fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
         (*kernel)[i] = fx;
         sum += fx;
      }

      for(i=0;i<(*windowsize);i++) (*kernel)[i] /= sum;

      if(VERBOSE){
         printf("The filter coefficients are:\n");
         for(i=0;i<(*windowsize);i++)
            printf("kernel[%d] = %f\n", i, (*kernel)[i]);
      }
   }

