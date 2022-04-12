/////////////////////////////////////////////////////
//   Filtered back projection algorithm code for CT reconstruction (for CPU)
//   Copyright (c) 2022 Akihiro Haga
//   The Tokushima University
//   Email: haga@tokushima-u.ac.jp
//   This software is released under the MIT License.
//   http://opensource.org/licenses/mit-license.php
//   Ref: Akihiro Haga, et al., Radiation Oncology, v9, 2014
//   -- Usage --
//   0. Prepare sinogram (projection image)
//      ex. reprojection_float.raw from "Release_MBFPA_CPU_human"
//   1. make
//   2. mvct.exe (input.txt)
//      ex. mvct.exe TOMO_input_virtual_projection_human.txt
//   3. output image is producted (as "FBP_virtual_projection_512x512_human.raw")
/////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "physParamsTOMO.h"
#include "mallocMD.h"
#include "misc.h"
#include "projectionDataTOMO.h"
#include "filterSinogramTOMO.h"
#include "reconstructImageTOMO.h"
//#include "reconstructImageOnGPU_TOMO.h"
//#define USE_GPU

int
main(int argc,char* argv[])
{
  // show usage unless command is OK
  if( argc != 2 )
    {
      fprintf(stderr,"\n[usage]$ cbct <info filename>\n\n");
      exit(0);
    }
  // all in one as data structure (projeVolume, angleRad, xshift, yshift)
  PROJECTION_DATA* projData = newProjectionData();
  
  // load data required for reconstuction (angle, projection images, etc.)
  if( loadData(argv[1],projData) == 0 )
    {
      deleteProjectionData( projData );
      exit(0);
    }
  
  
  FILE*	fp;
  
  if( (fp=fopen(projData->reconDataFilename,"wb")) == NULL )
    {
      fprintf(stderr,"ERROR: output file NOT writable\n");
      exit(0);
    }
  else
    {
      float**     fltdSinogram = (float**)new2DArray( projData->projImgWidth,
						      projData->projImgHeight*projData->usedProjNumber, UNIT_FLOAT32 );
      float***	reconImage = (float***)new3DArray( projData->reconSize,
						   projData->reconSize, projData->reconSlices, UNIT_FLOAT32 );
      short*** 	reconImageshort = (short***)new3DArray( projData->reconSize,
							projData->reconSize, projData->reconSlices, UNIT_SINT16 );
      double**	constMap = prepareConstMap( projData->projImgWidth, projData->ReconstructionFilter );
      
      double*     Corrected_Detector_Weight0 = (double*)new1DArray( End_Channel-Start_Channel,UNIT_FLOAT64 );
      double*     Corrected_Detector_Weight1 = (double*)new1DArray( End_Channel-Start_Channel,UNIT_FLOAT64 );
      double*     Corrected_Detector_Weight2 = (double*)new1DArray( End_Channel-Start_Channel,UNIT_FLOAT64 );
      prepareCorrectedDetectorInterpolation(End_Channel-Start_Channel, Corrected_Detector_Weight0, Corrected_Detector_Weight1, Corrected_Detector_Weight2);
      char**		reconMaskImage = prepareMaskData(projData->reconSize);
      double		*cosTable, *sinTable;
      //int**       intStartAngle = (int**)new2DArray( projData->reconSize, projData->reconSize, UNIT_FLOAT32 );
      //int**       intEndAngle = (int**)new2DArray( projData->reconSize, projData->reconSize, UNIT_FLOAT32 );
      //prepareAngles( projData->usedProjNumber, projData->anglesRad, projData->reconSize, intStartAngle, intEndAngle);
      
      prepareTrigonometricTable( projData->projImgHeight, projData->anglesRad, &cosTable, &sinTable );

      /*
	#ifdef USE_GPU
	initializeGPU( projData->reconSize, projData->reconSlices, 
	projData->projImgHeight, projData->projImgWidth, projData->usedProjNumber,
	projData->anglesRad, 
	cosTable, sinTable,
	projData->xshift,projData->yshift );
	#endif*/
      //printf("%d \n",projData->usedProjNumber);
      {
      clearReconstructedImage( projData->reconSize,
	projData->reconSlices, reconImage, reconImageshort);
      clearSinogramImage( projData->projImgWidth,
	projData->projImgHeight*projData->usedProjNumber, fltdSinogram);
      
      startTimer();
      
			
      /* #ifdef USE_GPU
	 filterSinogramOnGPU( projData->projImgHeight, projData->projImgWidth, projData->usedProjNumber,
	 projData->projVolume, fltdSinogram,
	 projData->xshift,
	 projData->yshift );
      */
      /*
	char	filename[256];
	sprintf(filename,"fltsinogram%dx%d.float32.raw",projData->projImgWidth,projData->usedProjNumber);
	FILE*	fp=fopen(filename,"wb");
	for(int th=0;th<projData->usedProjNumber;th++)
	for(int m=0;m<projData->projImgWidth;m++)
	{
	float flt=fltdSinogram[th][256][m];
	fwrite( &flt, getByteSizeOfUnit(UNIT_FLOAT32),1,fp);
	}
	fclose(fp);*/
      
      //			stopTimer();	showLapTime("filterSinogramOnGPU_TOMO");	startTimer();
      
      double zstart = projData->reconOffset;
      double zend = projData->reconOffset + (projData->reconSlices -1) * projData->reconStep;
      double s_start = projData->StartCouchPosition + projData->CouchSpeed/(UsedReconAngle) + projData->CouchSpeed - 0.0000001;
      double s_end = projData->StartCouchPosition + (projData->usedProjNumber) * projData->CouchSpeed - projData->CouchSpeed/(UsedReconAngle + 0.0000001);
      fprintf(stderr,"Allowable reconstraction range = %lf to %lf \n", s_start, s_end);
      fprintf(stderr,"Reconstructing z = %lf to %lf, Step = %f ... \n", zstart, zend, projData->reconStep);
      if( zstart < s_start || zend > s_end)
	{
      fprintf(stderr,"\n Wrong reconstruction range !!-- STOP\n\n");
      exit(0);
    }
      //int zslice_min = (int)((zstart - projData->StartCouchPosition)/projData->CouchSpeed+0.000001);
      //if(zslice_min < 1) zslice_min = 1;
      //int zslice_max = int((zend+projData->CouchSpeed/(UsedReconAngle)-projData->StartCouchPosition)/projData->CouchSpeed);
      //if(zslice_max > projData->usedProjNumber) zslice_max = projData->usedProjNumber;
      int zslice_min = 1;
      int zslice_max = 1;
      fprintf(stderr,"z slice num. = %d to %d, usedProjNum = %d\n", zslice_min, zslice_max, projData->usedProjNumber );
      //#else
      filterSinogram( projData->projImgHeight, zslice_min, zslice_max,
	projData->projVolume, fltdSinogram,
	constMap,
	Corrected_Detector_Weight0, Corrected_Detector_Weight1, Corrected_Detector_Weight2,
	projData->xshift, projData->yshift,
	projData->InPrClassString
	);
      
      stopTimer();	showLapTime("filterSinogram_TOMO");	startTimer();
      //#endif
      // clear and reconstruct image
      
      /*#ifdef USE_GPU
	fprintf(stderr,"%d ... %d \n",projData->xfactor,projData->projImgWidth);
			
	reconstructImageOnGPU( projData->reconSize, projData->reconOffset, projData->reconSlices, projData->reconStep, reconImage, reconMaskImage, projData->xfactor,
	projData->projImgHeight, projData->projImgWidth, projData->usedProjNumber, fltdSinogram, intStartAngle, intEndAngle );
	
	stopTimer();	showLapTime("reconstructImageOnGPU_TOMO");
	#else*/
            
      reconstructImage(projData->reconSize,
	reconImage, reconMaskImage,
	projData->projImgWidth, projData->projImgHeight,
	fltdSinogram,
	cosTable, sinTable,
	projData->reconOffset, projData->reconStep, projData->reconSlices,
	projData->StartCouchPosition, projData->CouchSpeed, projData->usedProjNumber,
	projData->xfactor, projData->reconPixSize,
	projData->InPrClassString
	);
      stopTimer();	showLapTime("reconstructImage_TOMO");
      
      //#endif
      /*			for(int k=0;k<projData->reconSlices;k++)
				for(int j=0;j<projData->reconSize;j++)
				for(int i=0;i<projData->reconSize;i++)
				{
				if(reconImage[k][j][i] < SHRT_MIN  || SHRT_MAX < reconImage[k][j][i]){
				printf("ERROR!!! overflow...\n");
				exit(0);
				}
				reconImageshort[k][j][i] = (short)reconImage[k][j][i];
				//printf("SHRT_MIN = %d, SHRT_MAX=%d\n", SHRT_MIN, SHRT_MAX);
				}
				fwrite( get1DArrayOf3DArray((void***)reconImageshort,UNIT_SINT16),
					getByteSizeOfUnit(UNIT_SINT16),
					projData->reconSize*projData->reconSize*projData->reconSlices,
			fp);
*/			
      fwrite( get1DArrayOf3DArray((void***)reconImage,UNIT_FLOAT32),
	getByteSizeOfUnit(UNIT_FLOAT32),
	projData->reconSize*projData->reconSize*projData->reconSlices,
	fp);
      
      fprintf(stderr,"done.\n");
    }
      fclose(fp);
      
      // Output of ReconVoxData //
      /*FILE*	fp1;
	
	if( (fp1=fopen("ReconVoxData.raw","wb")) == NULL )
	{
	fprintf(stderr,"ERROR: output file NOT writable\n");
	exit(0);
	}
	else
	{
	fwrite( get1DArrayOf3DArray((void***)reconVoxData,UNIT_FLOAT64),
	getByteSizeOfUnit(UNIT_FLOAT64),
	projData->reconSize*projData->reconSize*330,
	fp1);
	}
      */
      /////////////////////////// 
      delete2DArray( (void**)fltdSinogram, projData->projImgWidth,
		     projData->projImgHeight*projData->usedProjNumber, UNIT_FLOAT32 );
      delete3DArray( (void***)reconImage, projData->reconSize, projData->reconSize, projData->reconSlices, UNIT_FLOAT32 );
      delete3DArray( (void***)reconImageshort, projData->reconSize, projData->reconSize, projData->reconSlices, UNIT_SINT16 );
      delete2DArray( (void**)constMap, End_Channel-Start_Channel, End_Channel-Start_Channel, UNIT_FLOAT64 );
      delete2DArray( (void**)reconMaskImage, projData->reconSize, projData->reconSize, UNIT_UINT8 );
      delete1DArray( (void*)cosTable );
      delete1DArray( (void*)sinTable );
#ifdef USE_GPU
      terminateGPU();
#endif
}
  
  deleteProjectionData( projData );
}


