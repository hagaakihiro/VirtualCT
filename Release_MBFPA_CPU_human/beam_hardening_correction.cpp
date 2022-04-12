/////////////////////////////////////////////////////
//   Simple beam hardening correction algorithm code for herical projection (for CPU)
//   Written by Akihiro Haga
//   The Tokushima University
//   Email: haga@tokushima-u.ac.jp
//   Ref: Kai-Wen Li, et al., Physica Medica, v89, p182, 2021
//   -- Usage --
//   0. Prepare sinogram (projection image)
//      ex. reprojection_float.raw
//   1. Compile as gcc beam_hardening_correction.cpp
//   2. ./a.out "alpha" "beta"
//      ex. ./a.out 0.01 2.00 (for kV CT) or ./a.out 0.01 3.70 (for MV CT)
//   3. output image is producted (as "reprojection_float_cor.raw")
//   Note that "alpha" and "beta" could depend on CT geometory(SDD, SID, Det.size etc.)/protocol(angle interval etc.) as well as the photon energy and phantom size.
/////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char* argv[])
{
  ///////////////  Simple beam hardening correction ///////////////////

  if( argc != 3 ){
    fprintf(stderr,"\n[usage]$ no alpha and beta \n\n");
    exit(0);
  }
  double alpha = atof(argv[1]);
  double beta = atof(argv[2]);
  
  FILE*	fp;
  FILE* fp3;
  char filename[128];
  int num_det = 609;
  int num_angle = 800;
  int reprojection_size = num_det * num_angle;
  sprintf(filename,"reprojection_float.raw");
  float*   reprojection_float =  (float*)malloc(reprojection_size*sizeof(float) );
  if( (fp=fopen(filename,"rb")) == NULL )
    {
      fprintf(stderr,"data file: %s not found\n",filename);
      return -1;
    }  
  int readSize = fread( (void*)reprojection_float, sizeof(float), reprojection_size, fp );
  fclose(fp);
  double min_nii = 1000000;
  double min_org = 1000000;
  for(int i = 0; i < reprojection_size; i++)
    {
      if(min_org > reprojection_float[i]) min_org = reprojection_float[i];
      double yi = -log(reprojection_float[i]);
      //yi = yi - (alpha*yi+beta*yi*yi);
      //printf("%lf %lf %lf\n", yi, fabs(yi), pow(abs(yi),beta));
      yi = fabs(yi) + (alpha*pow(fabs(yi),beta));

      reprojection_float[i] = (float)exp(-yi);
      if(min_nii > reprojection_float[i]) min_nii = reprojection_float[i];
    }
  printf("%lf %lf \n", min_org, min_nii);
  sprintf(filename,"reprojection_float_cor.raw");
  fp3=fopen(filename,"wb");	      
  fwrite( reprojection_float,
	  sizeof(float), reprojection_size, fp3);
  fclose(fp3);
  free(reprojection_float);
  
  return 0;
}

