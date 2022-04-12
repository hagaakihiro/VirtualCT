/////////////////////////////////////////////////////
//   Material based forward projection algorithm code for herical projection (for CPU)
//   Written by Akihiro Haga 
//   The Tokushima University
//   Email: haga@tokushima-u.ac.jp
//   Ref: Kai-Wen Li, et al., Br J Radiol, 2021
//   -- Usage --
//   0. Prepare elementary density image and incident X-ray spectrum
//      ex. Weight_input_ICRU110.raw, Spectrum/PhaseSpace_photon_120kV_0cm.txt
//   1. make
//   2. mvct.exe (input.txt)
//      ex. mvct.exe IR_TOMO_input_virtual_projection_120kV.txt
//   3. output image is producted (as "reprojection_float.raw")
/////////////////////////////////////////////////////
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <sys/stat.h>//MHIR
#include <sys/types.h>//MHIR

#include "physParamsTOMO.h"
#include "mallocMD.h"
#include "misc_TOMO.h"
#include "projectionDataTOMO.h"
#include "filterSinogramTOMO.h"
#include "IR_ImageTOMO.h"
#include "lengthproj_TOMO.h"
#include "reprojection_TOMO.h"
//#include "projOnGPU_TOMO.h"
#include "virtual_projection.h"

//#define USE_GPU // Comment out when GPU is not used.

  // Uniform distribution
  double Uniform0( void ){
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
  }

  double rand_gauss0( double mu, double sigma ){
    double z = sqrt( -2.0*log(Uniform0()) ) * sin( 2.0*PI*Uniform0() );
    return mu + sigma*z;
  }

int gaussiannoise0(float* ProjImagefloat, int NN, double noise)
{
  double projection_min = 10000000;
  for(int i = 0; i < NN; i++)
    {
      if(projection_min > ProjImagefloat[i])
	{
	  projection_min = ProjImagefloat[i];
	}
    }
  printf("projection_min = %lf, used noise = %lf\n",projection_min, noise);
  double mu=0;
  for(int i = 0; i < NN; i++)
    {
      //double nstar = MAX_PIXVALUE * ProjImagefloat[i]; // nstar = number of photons
      double nstar = ProjImagefloat[i]; 
      //double sigma = projection_min/10; // noise > 10 dB
      double sigma = noise; // noise > 10 dB
      double psn = rand_gauss0(0,sigma);
      nstar = nstar + psn;
      ProjImagefloat[i] = nstar;
      if(nstar < 0)  ProjImagefloat[i] = projection_min/10;
      if(nstar > 1.0)  ProjImagefloat[i] = 1.0;
    }
  return 0;
}

void print_npf_float(float* npf_float, PROJECTION_DATA* projData)
{

  char filename[256];
  FILE*	fp2;
  sprintf(filename,"npf_float.raw");
  fp2=fopen(filename,"wb");
  
  fwrite( npf_float,
	  getByteSizeOfUnit(UNIT_FLOAT32),
	  projData->reconSize*projData->reconSize,fp2);
  fclose(fp2);
  
}
  
void initial_image_def(int reconZ, float* reconImagefloat_before, PROJECTION_DATA* projData, short* reconImageShort, int num_mat)
{

  // Initial image definition
  for(int n=0;n<num_mat;n++){
    for(int j=0;j<projData->reconSize;j++){ // Height
      for(int jj=0;jj<projData->reconSize;jj++){ // Width
	int jjjj =  n * projData->reconSize * projData->reconSize + j * projData->reconSize + jj;
	int jjjjj = reconZ * projData->reconSize * projData->reconSize + j*projData->reconSize + jj;
	//reconImagefloat_before[jjjj] = (float)projData->FBP_reconImageshort[jjjjj]*0.0716/1024+0.0716;
	reconImagefloat_before[jjjj] = (float)reconImageShort[jjjj]/1000;
	//reconImagefloat_before[jjjj] = (float)0.07;
      }//Width
    }//Height
  }
  
}

void calc_squdiff_quad(int ite, int num_ite, double& objective, double& quadratic, double objective_before, double quadratic_before, int& num_pro, int start_angle_num, int end_angle_num, int thinout, float* projection_float, float* reprojection_float, PROJECTION_DATA* projData)
{
  int thinstep = (int)((ite-1)/thinout);
  int thinnumb = (ite-1)-thinstep*thinout;  
  num_pro = thinnumb;
  int sinogram_thinout = int((end_angle_num - start_angle_num)/thinout);
  for(int j_thin = 0; j_thin < sinogram_thinout; j_thin++) // sinogram range for one slice;
    {
      for(int jj=Start_Channel;jj<End_Channel;jj++)	    // Width
	{
	  int jjjj = num_pro * projData->projImgWidth + jj;
	  double temp_for_obj = (double)(projection_float[jjjj] - reprojection_float[jjjj]);
	  objective += temp_for_obj * temp_for_obj;
	  if(ite == num_ite){
	    double temp_for_quad = log(projection_float[jjjj]/reprojection_float[jjjj]);
	    quadratic += 0.5 * projection_float[jjjj] * temp_for_quad * temp_for_quad;
	  }
	}
      num_pro = num_pro + thinout;
    }
  
  //if(ite == num_ite)
  if(ite==ite/100*100) printf("Square difference =  %lf, Quadratic = %lf \n", objective, quadratic);
}

void update_wbefore(int ite, double objective, double objective_before, double& wbefore)
{

  if(objective < objective_before){
    wbefore = wbefore + 0.0001;
  }
  else{
    if(ite > 1) wbefore = wbefore - 0.0002;
    //if(ite != 0) wbefore = wbefore - 0.00010;
  }
  if(wbefore <= 0.0) wbefore = 0.0001;
  //fprintf(stderr,"wbefore = %lf \n", wbefore);
  //wbefore = wbefore + ;
  
}

void calc_diff_copy_reconImflo_before(double& diff, int reconZ, float* reconImagefloat_before, float* reconImagefloat, PROJECTION_DATA* projData)
{

  for(int j=0;j<projData->reconSize;j++){ // Hight
    for(int jj=0;jj<projData->reconSize;jj++){ // Width
      int jjjj =  j * projData->reconSize + jj;
      int jjjjj = reconZ * projData->reconSize * projData->reconSize + j*projData->reconSize + jj;
      
      //reconImagefloat_before_3Dnew[jjjjj] = reconImagefloat[jjjj];
      diff+= (reconImagefloat_before[jjjj]-reconImagefloat[jjjj])
	*(reconImagefloat_before[jjjj]-reconImagefloat[jjjj]);
      reconImagefloat_before[jjjj] = reconImagefloat[jjjj];
    } // Width
  } // Height
  
}

void output_recon_image(int ite, int reconZ, float* reconImagefloat, PROJECTION_DATA* projData)
{
  //if(ite == ite/10000*10000){
  //if(ite == ite)
    {
    char filename[128];
    FILE*	fp3;
    sprintf(filename,"ReconImageRaw/Reconstruction_ite%04d_z%04d.raw", ite, reconZ);
    fp3=fopen(filename,"wb");
    
    fwrite( reconImagefloat,
	    getByteSizeOfUnit(UNIT_FLOAT32),
	    projData->reconSize*projData->reconSize,fp3);
    fclose(fp3);
  }

}

//void calc_reconImagefloat_min(double& min_quadra, double quadratic, float* reconImagefloat_min, float* reconImagefloat, PROJECTION_DATA* projData)
void calc_reconImagefloat_min(double& min_squdiff, double objective, float* reconImagefloat_min, float* reconImagefloat, PROJECTION_DATA* projData)
{

  if(min_squdiff > objective){
    min_squdiff = objective;
    //printf("min_squdiff = %f \n",min_squdiff);
    for(int j=0;j<projData->reconSize;j++){ // Hight
      for(int jj=0;jj<projData->reconSize;jj++){ // Width
	int jjjj =  j * projData->reconSize + jj;
	reconImagefloat_min[jjjj] = reconImagefloat[jjjj];
      }
    }
  }
  
//  if(min_quadra > quadratic){
//    min_quadra = quadratic;
//    for(int j=0;j<projData->reconSize;j++){ // Hight
//      for(int jj=0;jj<projData->reconSize;jj++){ // Width
//	int jjjj =  j * projData->reconSize + jj;
//	reconImagefloat_min[jjjj] = reconImagefloat[jjjj];
//      }
//    }
//  }
  
}

int main(int argc,char* argv[])
{
  // show usage unless command is OK
  if( argc != 2 ){
    fprintf(stderr,"\n[usage]$ mvct <info filename>\n\n");
    exit(0);
  }
  
  PROJECTION_DATA* projData = newProjectionData();
  
  // load data required for reconstuction (input data, prior image, etc.)
  if( loadData(argv[1],projData) == 0 ){
    deleteProjectionData( projData );
    exit(0);
  }
  
  
  FILE*	fp;
  if( (fp=fopen(projData->reconDataFilename,"rb")) == NULL ){
    fprintf(stderr,"No prior image \n");
    //exit(0);
  }
  else{
    fprintf(stderr,"Prior image: %s \n", projData->reconDataFilename);
  }
  
  int num_material = 6; // number of element; here, (H, C, N, O, P, Ca) See attenuation_coefficient
  int recon2Dsize = projData->reconSize*projData->reconSize;
  //int recon3Dsize = projData->reconSize*projData->reconSize*projData->reconSlices;
  int proje2Dsize = projData->projImgWidth*projData->projImgHeight;
  //int proje3Dsize = projData->projImgWidth*projData->projImgHeight*projData->usedProjNumber;
  int reconMDsize = projData->reconSize*projData->reconSize*num_material;
  // Data allocation
  float*   reprojection_float = (float*)new1DArray( proje2Dsize, UNIT_FLOAT32 );
  float*   projection_float = (float*)new1DArray( proje2Dsize, UNIT_FLOAT32 );
  float*   npf_float = (float*)new1DArray( recon2Dsize, UNIT_FLOAT32 );
  float*   npf_float_re = (float*)new1DArray( recon2Dsize, UNIT_FLOAT32 );
  float*   Prior_image_2D =  (float*)new1DArray(recon2Dsize, UNIT_FLOAT32 );
  float*   reconImagefloat_before = (float*)new1DArray(reconMDsize, UNIT_FLOAT32 );
  //float*   reconImagefloat_before_3D = (float*)new1DArray( recon3Dsize, UNIT_FLOAT32 );
  //float*   reconImagefloat_before_3Dnew = (float*)new1DArray( recon3Dsize, UNIT_FLOAT32 );
  short*   reconImageShort = (short*)new1DArray( reconMDsize, UNIT_SINT16 );
  float*   reconImagefloat = (float*)new1DArray( recon2Dsize, UNIT_FLOAT32 );
  float*   reconImagefloat_min = (float*)new1DArray( recon2Dsize, UNIT_FLOAT32 );

  float*   systemMatrixfloat = (float*)new1DArray( projData->projImgHeight*recon2Dsize, UNIT_FLOAT32 );
  float*   xraytracedetector = (float*)new1DArray( projData->projImgHeight*recon2Dsize, UNIT_FLOAT32 );
  float*   reprojection_systemMatrixfloat = (float*)new1DArray( proje2Dsize*projData->reconSize*2, UNIT_FLOAT32 );
  int*     reprojection_vnum = (int*)new1DArray( proje2Dsize*projData->reconSize*2, sizeof(int) );
  int*     reprojection_sum = (int*)new1DArray( proje2Dsize, sizeof(int) );
    
  double*  Corrected_Detector_Weight0 = (double*)new1DArray( End_Channel-Start_Channel,UNIT_FLOAT64 );
  double*  Corrected_Detector_Weight1 = (double*)new1DArray( End_Channel-Start_Channel,UNIT_FLOAT64 );
  double*  Corrected_Detector_Weight2 = (double*)new1DArray( End_Channel-Start_Channel,UNIT_FLOAT64 );

  char*	   reconMaskImage = prepareMaskData(projData->reconSize);
  double   *cosTable, *sinTable;
//#ifdef USE_GPU
//  //////////////by yuya nemoto///////////////////
//  float*   projection_floatGPU = (float*)new1DArray( (End_Channel-Start_Channel)*projData->projImgHeight, UNIT_FLOAT32 );
//#endif
  // Interpolation required due to the discrepancy between the detector curvature and the rotation one // filterSinogramTOMO.cpp
  prepareCorrectedDetectorInterpolation(End_Channel-Start_Channel, Corrected_Detector_Weight0, Corrected_Detector_Weight1, Corrected_Detector_Weight2);
  // cos and sin tables used in the reconstruction // calculated once here // misc_TOMO.cpp
  prepareTrigonometricTable( projData->projImgHeight, projData->anglesRad, &cosTable, &sinTable );
		
				
  fprintf(stderr,"usedProj: %d SngrmH: %d SngrmW: %d\n", projData->usedProjNumber,  projData->projImgHeight, projData->projImgWidth);		
  struct stat st;
  if(stat("ReconImageRaw", &st) != 0)
    {
      mkdir("ReconImageRaw", 0775);
    }

//#ifdef USE_GPU
//  initializeGPU( projData->reconSize, projData->reconPixSize,
//		 projData->projImgHeight, projData->projImgWidth, 
//		 cosTable, sinTable, npf_float,  npf_float_re, reprojection_float, num_material);
//  //		 cosTable, sinTable, npf_float,  npf_float_re, reprojection_float, projection_floatGPU);
//#endif
  {
    
    char	filename[256];
  
    double zstart = projData->reconOffset; // reconstruction start Z
    double zend = projData->reconOffset + (projData->reconSlices -1) * projData->reconStep; // reconstruction end Z
    double s_start = projData->StartCouchPosition + projData->CouchSpeed/(UsedReconAngle); // possible start limit for Z
    double s_end = projData->StartCouchPosition + (projData->usedProjNumber+1) 
      * projData->CouchSpeed - projData->CouchSpeed/(UsedReconAngle); // possible end limit for Z
    //if( zstart + 0.001 < s_start || zend - 0.001 > s_end){
    //  fprintf(stderr,"\n Wrong reconstruction range !!-- s_start %lf < zstart %lf ?  STOP\n\n", s_start, zstart);
    //  fprintf(stderr,"\n Wrong reconstruction range !!-- s_end %lf > zend %lf ?  STOP\n\n", s_end, zend);
    //  exit(0);
    //}
    fprintf(stderr,"Allowable reconstraction range = %lf to %lf \n", s_start, s_end);
    fprintf(stderr,"Reconstructing z = %lf to %lf, Step = %lf [cm] ... \n", zstart, zend, projData->reconStep);
    
    int thinout = projData->osemset;
    //fprintf(stderr, "*** Thin out = %d *** \n", thinout);
    
    double diff; // image difference by iteration
    
    char filename_recon[128];
    FILE*	fprecon;
    sprintf(filename_recon,"ReconImageRaw/Reconstruction_512x512x%d.raw", projData->reconSlices);
    fprecon=fopen(filename_recon,"wb");
    
    //z loop zstart
    //for(int reconZ=0;reconZ<projData->reconSlices; reconZ++){
    for(int reconZ=0;reconZ<1; reconZ++){
      
//      char Iteration_data_name[128];
//      sprintf(Iteration_data_name,"iteration_data_%d.txt",reconZ);
//      FILE*	fpp;
//      if( (fpp=fopen(Iteration_data_name,"w")) == NULL ){
//	fprintf(stderr,"No Record \n");
//	exit(0);
//      }
      
      clear_npf( recon2Dsize, npf_float);

      // Energy spectrum loading
      char spectrum_data_name[128];
      sprintf(spectrum_data_name,"%s",projData->SpectrumData);
      float* Attenuation = (float*)new1DArray( 1024*num_material, sizeof(float) );
      float* pdf_pe = (float*)new1DArray( 1024, sizeof(float) );
      int NE;
      int nthin=1;
      // Attenuation coefficient for each material and each energy -- Attenuation
      attenuation_coefficient( spectrum_data_name, num_material, Attenuation, pdf_pe, &NE, nthin);

      //projection used in recon Z
      double z_position = zstart + (double)reconZ*projData->reconStep; // reconstruction Z position
      int sngrmHeight=projData->projImgHeight; // sinogram height (MVCT uses 800 projections (angles))
      double detector_Z_pitch = (double)(projData->CouchSpeed/sngrmHeight); // Couch movement in [cm] for a projection
      //int start_angle_num = (int)(z_position/detector_Z_pitch-sngrmHeight/(UsedReconAngle)+0.00000001); // Sinogram No. at start
      int start_angle_num = 0;
      int end_angle_num = start_angle_num + sngrmHeight/(UsedReconAngle)*2; // Sinogram No. at end
      int usesinogram = int(sngrmHeight/(UsedReconAngle)*2); // Num. of sinogram used in the reconstruction // original: 800 (360 degrees)
      printf("%d \n",end_angle_num);
      
      // Sinogram data interpolation due to helical geometry
      int num_pro = 0;
      for(int j = start_angle_num; j < end_angle_num; j++){ // sinogram range used;
	int iphase = 0;
	if(num_pro < usesinogram/2) iphase = 1;
	else if( num_pro > usesinogram/2) iphase = -1;
	if( j+iphase*sngrmHeight  < sngrmHeight ) iphase = 0;
	if( j+iphase*sngrmHeight  > (projData->usedProjNumber-1)*sngrmHeight ) iphase = 0;
	double z_weight = 0.5*fabs(num_pro-(sngrmHeight/2.0))/(sngrmHeight/2.0);
	for(int jj = 0; jj < projData->projImgWidth; jj++){
	  //projection_float[num_pro*projData->projImgWidth + jj] = projData->projVolume[j*projData->projImgWidth + jj];
	  projection_float[num_pro*projData->projImgWidth + jj] = 
	    (1.0-z_weight)*projData->projVolume[j*projData->projImgWidth + jj] 
	    + z_weight*projData->projVolume[(j+iphase*sngrmHeight)*projData->projImgWidth+jj];
	}
	num_pro++;
      }
      /*
      FILE*	fpproj;
      sprintf(filename,"projection_float.raw"); // Sinogram used in the reconstruction at reconZ
      fpproj=fopen(filename,"wb");
      
      fwrite( projection_float,
	      getByteSizeOfUnit(UNIT_FLOAT32),
	      projData->projImgHeight*projData->projImgWidth,fpproj);
      fclose(fpproj);
      */
      double min_squdiff = 10000000.0;
      //double min_quadra = 10000000.0;
      double wbefore = projData->gradientdecentconst+0.0005;
      //  Iteration loop	    
      double objective = 10000000.0;
      double quadratic = 10000000.0;
      double objective_before;
      double quadratic_before;
      //for(int ite = 0; ite < projData->num_ite + 1; ite++){
      for(int ite = 0; ite < 2; ite++){
	if(ite==ite/100*100) fprintf(stderr,"Iteration loop, %d step,	\n",ite);
        
	if(ite == 0){ // npf creation and initial image definition
	  double xalpha[1];
//#ifdef USE_GPU
//	  //	  lengthprojOnGPU_TOMO( ite, projData->reconSize, npf_float, npf_float_re, projection_float, projection_floatGPU,
//	  lengthprojOnGPU_TOMO( ite, projData->reconSize, npf_float, npf_float_re, projection_float, 
//				projData->reconPixSize,
//				projData->projImgWidth, projData->projImgHeight,
//				z_position,
//				projData->CouchSpeed,
//				reconMaskImage, systemMatrixfloat, xraytracedetector, thinout);
//	  
//#else
	  //clock_t start = clock();
	  lengthproj_TOMO( ite, projData->reconSize, npf_float, projection_float,
			   projData->reconPixSize,
			   projData->projImgWidth, projData->projImgHeight,
			   cosTable, sinTable,
			   z_position,
			   projData->CouchSpeed,
			   reconMaskImage, systemMatrixfloat, xraytracedetector,
			   1);
	  //clock_t end = clock();     // 
	  //std::cout << "len_TOMO duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
//#endif
	  
	  print_npf_float(npf_float, projData);

	  FILE*	input_image;
	  sprintf(filename,"Weight_input_ICRU110.raw");
	  input_image=fopen(filename,"rb");
	  size_t a;
	  a = fread( reconImageShort,
		  getByteSizeOfUnit(UNIT_SINT16),
		  reconMDsize,input_image);
	  fclose(input_image);
	  initial_image_def(reconZ, reconImagefloat_before, projData, reconImageShort, num_material);
	  //projData->xfactor = 10000;
	}//if ite == 0
	
	else{ //if ite != 0
	  double xalpha[1];
	  
	  clearSinogramImage( proje2Dsize, reprojection_float);
	  // reprojection with reconImagefloat_before
//#ifdef USE_GPU
//	  reprojectionOnGPU_TOMO( ite, projData->reconSize, reconImagefloat_before, reprojection_float,
//				  projData->reconPixSize,
//				  projData->projImgWidth, projData->projImgHeight,
//				  z_position,
//				  projData->CouchSpeed,
//				  projData->xfactor, reprojection_systemMatrixfloat,
//				  reprojection_vnum, reprojection_sum, thinout,
//				  Attenuation, pdf_pe, num_material, NE);
//	  FILE*	reproj;
//	  sprintf(filename,"reprojection_float_gpu.raw");
//	  reproj=fopen(filename,"wb");
//	  float projection_max = 0;
//	  for(int i = 0; i < projData->projImgHeight*projData->projImgWidth; i++)
//	    {
//	      if(projection_max < reprojection_float[i])
//		{
//		  projection_max = reprojection_float[i];
//		}
//	    }
//	  for(int i = 0; i < projData->projImgHeight*projData->projImgWidth; i++)
//	    {
//	      reprojection_float[i] = reprojection_float[i]/projection_max;
//	    }
//	  // Gaussian Noise Input
//	  //gaussiannoise0(reprojection_float,projData->projImgHeight*projData->projImgWidth,projData->photonnoise );
//
//	  fwrite( reprojection_float,
//		  getByteSizeOfUnit(UNIT_FLOAT32),
//		  projData->projImgHeight*projData->projImgWidth,reproj);
//	  fclose(reproj);
//	  		
//	  
//#else
	  if(ite == 1){
	    //clock_t start = clock();
	    float*    reprojection_float_full = (float*)new1DArray( projData->projImgWidth*projData->projImgHeight, UNIT_FLOAT32 );
	    /*
	    //for normal ct projection
	    NE = 1;
	    num_material = 1;
	    Attenuation[0] = 1;
	    pdf_pe[0] = 1;
	    for(int mat = 1; mat < num_material; mat++)
	      {
		for(int i = 0; i < projData->reconSize*projData->reconSize; i++) reconImagefloat_before[i] += reconImagefloat_before[mat*projData->reconSize*projData->reconSize+i];
	      }
	    for(int i = 0; i < projData->reconSize*projData->reconSize; i++)
	      {
	    	reconImagefloat_before[i] = reconImagefloat_before[i]*2.1;
	    	if(reconImagefloat_before[i]>0.0001) printf("%f \n", reconImagefloat_before[i]);
	      }
	    //exit(0);
	    */
	    for(int ke = 0; ke < NE; ke++)
	      {
		reprojection_TOMO( ite, projData->reconSize, reconImagefloat_before, reprojection_float,
				   projData->reconPixSize,
				   projData->projImgWidth, projData->projImgHeight,
				   cosTable, sinTable,
				   z_position,
				   projData->CouchSpeed,
				   projData->xfactor, xalpha, reconMaskImage, reprojection_systemMatrixfloat,
				   reprojection_vnum, reprojection_sum,
				   1, Attenuation, pdf_pe, num_material, NE, ke);

		for(int i = 0; i < projData->projImgWidth*projData->projImgHeight; i++) reprojection_float_full[i] += reprojection_float[i];

	      }
	    
	    FILE*	reproj;
	    sprintf(filename,"reprojection_float.raw");
	    reproj=fopen(filename,"wb");
	    float projection_max = 0;
	    for(int i = 0; i < projData->projImgHeight*projData->projImgWidth; i++)
	    {
	      if(projection_max < reprojection_float_full[i])
		{
		  projection_max = reprojection_float_full[i];
		}
	    }
	    for(int i = 0; i < projData->projImgHeight*projData->projImgWidth; i++)
	    {
	      reprojection_float_full[i] = reprojection_float_full[i]/projection_max;
	    }
	    // Gaussian Noise Input
	    //gaussiannoise0(reprojection_float_full,projData->projImgHeight*projData->projImgWidth);
	    gaussiannoise0(reprojection_float_full,projData->projImgHeight*projData->projImgWidth,projData->photonnoise );

	    fwrite( reprojection_float_full,
		  getByteSizeOfUnit(UNIT_FLOAT32),
		  projData->projImgHeight*projData->projImgWidth,reproj);
	    fclose(reproj);
	    //clock_t end = clock();
	    //std::cout << "rep_TOMO duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
	  }
	  else{
	    //clock_t start = clock();
	    for(int ke = 0; ke < NE; ke++)
	      {
		reprojection_TOMO_2( ite, projData->reconSize, reconImagefloat_before, reprojection_float,
				 projData->reconPixSize,
				 projData->projImgWidth, projData->projImgHeight,
				 cosTable, sinTable,
				 z_position,
				 projData->CouchSpeed,
				 projData->xfactor, xalpha, reconMaskImage, reprojection_systemMatrixfloat,
				 reprojection_vnum, reprojection_sum,
				 thinout, Attenuation, pdf_pe, num_material, NE, ke);
	      }
	    //clock_t end = clock(); 
	    //std::cout << "rep_TOMO_2 duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
	  }
//#endif
	  
	  objective_before = objective;
	  quadratic_before = quadratic;
	  objective = 0.0;
	  quadratic = 0.0;

	  calc_squdiff_quad(ite, projData->num_ite, objective, quadratic, objective_before, quadratic_before, num_pro, start_angle_num, end_angle_num, thinout, projection_float, reprojection_float, projData);
	  clear_npf( recon2Dsize, npf_float_re);
	  
	  // Estimated npf
//#ifdef USE_GPU
//	  //	  lengthprojOnGPU_TOMO( ite, projData->reconSize, npf_float, npf_float_re, reprojection_float, projection_floatGPU,
//	  lengthprojOnGPU_TOMO( ite, projData->reconSize, npf_float, npf_float_re, reprojection_float, 
//				projData->reconPixSize,
//				projData->projImgWidth, projData->projImgHeight,
//				z_position,
//				projData->CouchSpeed,
//				reconMaskImage,
//				systemMatrixfloat, xraytracedetector, thinout );
//
//#else
	  //clock_t start = clock();
	  lengthproj_TOMO_2( ite, projData->reconSize, npf_float_re, reprojection_float,
			     projData->reconPixSize,
			     projData->projImgWidth, projData->projImgHeight, 
			     z_position,
			     projData->CouchSpeed,						       
			     reconMaskImage,
			     systemMatrixfloat, xraytracedetector,
			     thinout);
	  //clock_t end = clock();     // 
	  //std::cout << "len_TOMO_2 duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
//#endif
		
	  // Weight
	  double weight_tv = projData->totalvariation;//0.05, 0.5  // Weight for total variation
	  double weight_prior = projData->priorconst;     // Weight for prior image constraint
	  double weight_marko = projData->markovconst;

	  update_wbefore(ite, objective, objective_before, wbefore);

#ifdef USE_GPU
	  IR_ImageOnGPUTOMO( ite, projData->reconSize, reconImagefloat_before, reconImagefloat,
			     npf_float_re, wbefore,
			     weight_tv, weight_prior, thinout);
#else
	  IR_ImageTOMO( ite, projData->reconSize, reconImagefloat_before, reconImagefloat,
			Prior_image_2D, npf_float, npf_float_re, wbefore,
			weight_tv, weight_prior, reconMaskImage, thinout);
#endif		

	  calc_diff_copy_reconImflo_before(diff, reconZ, reconImagefloat_before, reconImagefloat, projData);

	  // Output of Reconstructed Image //
	  output_recon_image(ite, reconZ, reconImagefloat, projData);

	  calc_reconImagefloat_min(min_squdiff, objective, reconImagefloat_min,
				   reconImagefloat, projData);
//	  calc_reconImagefloat_min(min_quadra, quadratic, reconImagefloat_min,
//				   reconImagefloat, projData);
	  
	}//if ite !=0
	    
      } // iteration loop end

      for(int j=0;j<projData->reconSize*projData->reconSize;j++)        // Hight
	{
	  reconImageShort[j] = (short)(1024*(reconImagefloat_min[j]-0.0716)/0.0716);
	}
      if(fprecon)//MHIR	
	fwrite( reconImageShort,
		getByteSizeOfUnit(UNIT_SINT16),
		projData->reconSize*projData->reconSize,fprecon);
      //	fwrite( reconImagefloat_min,
      //		getByteSizeOfUnit(UNIT_FLOAT32),
      //		projData->reconSize*projData->reconSize,fprecon);
//      if(fpp){  //MHIR
//	fclose(fpp);    
//	fpp = NULL;
//      }
    } // reconZ end
    if(fprecon){ //MHIR
      fclose(fprecon);
      fprecon = NULL;
    }
    fprintf(stderr,"done.\n");
    if(fp){ //MHIR
      fclose(fp);
      fp = NULL;
    }
    /////////////////////////// 
    
    delete1DArray( (void*)reprojection_float );
    delete1DArray( (void*)projection_float );
    delete1DArray( (void*)npf_float );
    delete1DArray( (void*)npf_float_re );
    delete1DArray( (void*)Prior_image_2D );
    delete1DArray( (void*)reconImagefloat_before );
    //delete1DArray( (void*)reconImagefloat_before_3D );
    //delete1DArray( (void*)reconImagefloat_before_3Dnew );
    delete1DArray( (void*)reconImagefloat );
    delete1DArray( (void*)reconImagefloat_min );
    delete1DArray( (void*)systemMatrixfloat );
    delete1DArray( (void*)xraytracedetector );
    delete1DArray( (void*)reprojection_systemMatrixfloat );
    delete1DArray( (void*)reprojection_vnum );
    delete1DArray( (void*)reprojection_sum );
    delete1DArray( (void*)Corrected_Detector_Weight0);
    delete1DArray( (void*)Corrected_Detector_Weight1);
    delete1DArray( (void*)Corrected_Detector_Weight2);
    delete1DArray( (void*)reconMaskImage );
    delete1DArray( (void*)cosTable );
    delete1DArray( (void*)sinTable );

//#ifdef USE_GPU
//    delete1DArray( (void*)projection_floatGPU );
//#endif
//#ifdef USE_GPU
//    terminateGPU();
//#endif
  }
  
  deleteProjectionData( projData );
}


