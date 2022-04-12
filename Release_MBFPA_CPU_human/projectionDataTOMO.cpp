#include <stdio.h>
#include <string.h>
#include <stdlib.h> // for MAC
#include <malloc.h> // for Ubuntu

#include "physParamsTOMO.h"
#include "mallocMD.h"
#include "projectionDataTOMO.h"


#define LINEBUF_SIZE	512

static char	lineBuf[LINEBUF_SIZE];

//////

static void
skipLines(FILE* fp,int nLines)
{
  char* cread;
  for(int i=0;i<nLines;i++)	cread = fgets(lineBuf,LINEBUF_SIZE,fp);
}

static void
loadInteger(FILE* fp,int* value)
{
  char* cread;
  skipLines(fp,2);	cread = fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%d",value);
}

static void
loadfloat(FILE* fp,float* value)
{
  char* cread;
  skipLines(fp,2);	cread = fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%f",value);
}

static void
loaddouble(FILE* fp,double* value)
{
  char* cread;
  skipLines(fp,2);	cread = fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%lf",value);
}

static void
loadDouble3(FILE* fp, double* value1, double* value2, double* value3)
{
  char* cread;
  skipLines(fp,2);	cread = fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%lf %lf %lf",value1,value2,value3);
}

static void
loadIntegerHex(FILE* fp,int* value)
{
  char* cread;
  skipLines(fp,2);	cread = fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%X",value);
}

static void
loadString(FILE* fp,char* string)
{
  char* cread;
  skipLines(fp,2);	cread = fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%s",string);
}

//////

int
loadData(char* infoFilename,PROJECTION_DATA* projData)
{
  FILE*	fp;

  if( (fp=fopen(infoFilename,"r")) == NULL )
    {
      fprintf(stderr,"info file: %s not found\n",infoFilename);
      return 0;
    }
  
  // skip comments (6lines)
  skipLines(fp,6);
	
  // load each element
  loadString( fp, projData->InPrClassString );		// In-treatment/Pre-treatment label
  loadInteger( fp, &projData->Prior_image_ONOFF );	// Prior image ON/OFF (0 : no prior, 1 : planCT, 2 : MVCT)
  loadInteger( fp, &projData->projImgWidth );           // sinogram image width
  loadInteger( fp, &projData->projImgHeight );		// sinogram image height
  loadString( fp, projData->projImgDataType );          // Data type
  loadInteger( fp, &projData->projImgFileOffset );	// projection image offset
  loadString( fp, projData->projImgFileName );		// projection image file name
  loadString( fp, projData->AirScanImgFileName );	// Air scan image file name
  loaddouble( fp, &projData->GantryRotationPeriod );	// Gantry rotation period [s]
  loaddouble( fp, &projData->DataSamplingRate );        // Data Sampling Rate [Hz]
  loaddouble( fp, &projData->StartViewAngle );          // Starting View Angle [degree]
  loaddouble( fp, &projData->CouchSpeed );              // Couch Speed [cm/rotation]
  loaddouble( fp, &projData->StartCouchPosition );      // Start couch position [cm]
  loadInteger( fp, &projData->reconSize );		// reconstruction size (square only)
  loaddouble( fp, &projData->reconPixSize );		// reconstruction pixel size in cm
  loadString( fp, projData->reconTgtClassString );	// reconstruction target class label
  loaddouble( fp, &projData->reconOffset );         	// reconstruction start z position in cm
  loaddouble( fp, &projData->reconStep );		// reconstruction z step in cm
  loadInteger( fp, &projData->reconSlices );            // number of reconstructed slices
  loadInteger( fp, &projData->xfactor );	        // factor for CT value
  loadString( fp, projData->reconDataFilename );	// reconstruction file name (prior image)
  loadInteger( fp, &projData->FBP_slice );		// number of reconstructed slices
  loadDouble3( fp, &projData->x_reg, &projData->y_reg, &projData->z_reg );// Prior image shift
  ////////// This is required in Delivery CT (in-treatment CT)
  double ScanListZvalue;
  loaddouble( fp, &ScanListZvalue );			// reconstruction start point for z
  projData->z_reg = ScanListZvalue - projData->CouchSpeed * 1.5 + projData->reconOffset + projData->z_reg;
  //////////
  loadInteger( fp, &projData->osemset );		// Thin out for OSEM
  loadDouble3( fp, &projData->totalvariation, &projData->priorconst, &projData->markovconst );// Constants for regularization terms weight
  loaddouble( fp, &projData->gradientdecentconst );	// reconstruction z step in cm
    
  loadInteger( fp, &projData->num_ite );		// number of Iteration

  loadInteger( fp, &projData->sino_fact );		// sinogram factor (sinogram used in reconstruction)

  loadString( fp, projData->SpectrumData );     	// Spectrum data file name
  
  loaddouble( fp, &projData->photonnoise );     	// Photon noise in detector
  
  fclose(fp);

  // Data sampling rate // Thinning out for "IN-treatment"
  projData->DataSamplingRateOrig = projData->DataSamplingRate;
  if(projData->DataSamplingRate >= 200)
    {
      projData->DataSamplingRate = projData->DataSamplingRate/5;
    }
    
  projData->totalProjNumber = (int)(projData->projImgHeight/
				    (projData->GantryRotationPeriod * projData->DataSamplingRateOrig));
  projData->projImgHeight = (int)(projData->GantryRotationPeriod * projData->DataSamplingRate);
    
  if(strcmp(projData->projImgDataType,"float")==0)
    {
        projData->datasize = sizeof(float);
        projData->datatype = UNIT_FLOAT32;
    }
	else if(strcmp(projData->projImgDataType,"unsigned short")==0)
    {
        projData->datasize = sizeof(unsigned short);
        projData->datatype = UNIT_UINT16;
    }
	else if(strcmp(projData->projImgDataType,"signed short")==0)
    {
        projData->datasize = sizeof(short);
        projData->datatype=UNIT_SINT16;
    }
    else if(strcmp(projData->projImgDataType,"double")==0)
    {
        projData->datasize = sizeof(double);
        projData->datatype=UNIT_FLOAT64;
    }
    
  int readSize;
  // Prior image input //
  short* Preinput = (short*)new1DArray(projData->reconSize*projData->reconSize,UNIT_SINT16);
  if(projData->Prior_image_ONOFF==1) // planning CT (DICOM format)
    {
      projData->FBP_reconImageshort = (short*)new1DArray( pCTsize*pCTsize*projData->FBP_slice,
						      UNIT_SINT16 ); //// pCT input with original size
      projData->Pri_reconImageshort = (short*)new1DArray( projData->reconSize*projData->reconSize*projData->reconSlices,
						      UNIT_SINT16 ); //// pCT input with modified size
      char CTname[128];
      for(int i = 0; i < projData->FBP_slice; i++)
        {
	  sprintf(CTname, "%s/CT%d.dcm", projData->reconDataFilename,i+1);
	  fprintf(stderr, "%s\n",CTname);
	  if( (fp=fopen(CTname,"rb")) == NULL )
            {
	      fprintf(stderr,"CT file: %s not found\n",CTname);
	      return 0;
            }
	  unsigned char eTag1, eTag2, gTag1, gTag2, buf[4];
	  int allendflag = 0;
	  while( allendflag == 0)
            {
	      readSize = fread( buf, sizeof(unsigned char), 1, fp);
	      if( readSize != 1 )
		{
		  fprintf(stderr,"WARNING: data size unmatch\n");
		}
	      eTag2 = eTag1;
	      eTag1 = gTag2;
	      gTag2 = gTag1;
	      gTag1 = *(unsigned char*) buf;
	      if(eTag2 == 0xE0 && eTag1 == 0x7F && gTag2 == 0x10 &&gTag1 == 0x00)
                {
		  readSize = fread(buf, sizeof(unsigned char), 4, fp);
		  if( readSize != 4 )
		    {
		      fprintf(stderr,"WARNING: data size unmatch\n");
		    }
		  int dLen = *(int*)buf;
		  //fprintf(stderr,"%d\n",dLen);
		  readSize = fread(Preinput,1,dLen,fp);
		  if( readSize != dLen )
		    {
		      fprintf(stderr,"WARNING: data size unmatch\n");
		    }
		  //fread( projData->FBP_reconImageshort, sizeof(short), pCTsize * pCTsize, fp );
		  allendflag = 1;
                }
            }
	  for(int k = 0; k < pCTsize*pCTsize; k++)
            {
	      projData->FBP_reconImageshort[k+i*pCTsize*pCTsize] = Preinput[k];
            }
	  
        }
      fclose(fp);
      // END; pCT read

      // pCT input with modified size (projData->Pri_reconImageshort)
      for(int jjj=0;jjj<projData->reconSlices;jjj++)
	for(int j=0;j<projData->reconSize;j++)              // Hight
	  for(int jj=0;jj<projData->reconSize;jj++)	    // Width
	    {
	      double x1 = (-projData->reconSize/2+jj)*projData->reconPixSize,
		y1 = (-projData->reconSize/2+j)*projData->reconPixSize,
		z1 = (-projData->reconSlices/2+jjj)*projData->reconStep;
	      
	      int ppx = (x1+projData->x_reg)/pCTpixelsize + pCTsize/2,
		ppy = (y1-projData->y_reg)/pCTpixelsize + pCTsize/2,
		ppz = (z1+projData->z_reg)/pCTw + projData->FBP_slice/2;
	      
	      int jjjj = jjj * projData->reconSize * projData->reconSize + j * projData->reconSize + jj;
	      int jjjjj = ppz*pCTsize*pCTsize + ppy*pCTsize + ppx;///////////////////////////////////
	      projData->Pri_reconImageshort[jjjj] = 0;
	      if(ppx >= 0 && ppx < pCTsize && ppy >= 0 && ppy < pCTsize && ppz >= 0 && ppz < projData->FBP_slice)
		{
		  projData->Pri_reconImageshort[jjjj] = projData->FBP_reconImageshort[jjjjj]+1024;
		  if(projData->Pri_reconImageshort[jjjj] <= 0) projData->Pri_reconImageshort[jjjj] = 1;
		}
	    }
      //fprintf(stderr,"%lf %d \n", pCTpixelsize, pCTsize);
      projData->x_reg = 0.0;
      projData->y_reg = 0.0;
      projData->z_reg = 0.0;
      if( (fp=fopen("PreCT.raw","wb")) == NULL )
        {
	  fprintf(stderr,"PreCT file: not found\n");
	  return 0;
        }
      //fwrite(projData->FBP_reconImageshort, sizeof(short), pCTsize * pCTsize * projData->FBP_slice, fp);
      fwrite(projData->Pri_reconImageshort, sizeof(short),
	     projData->reconSize * projData->reconSize * projData->reconSlices, fp);
      fclose(fp);
    }
  else if(projData->Prior_image_ONOFF==2) // MVCT (FBP reconstruction)
    {
      projData->FBP_reconImageshort = (short*)new1DArray( pMVCTsize*pMVCTsize*projData->FBP_slice,
						      UNIT_SINT16 ); //// pMVCT input
      projData->Pri_reconImageshort = (short*)new1DArray( projData->reconSize*projData->reconSize*projData->reconSlices,
						      UNIT_SINT16 ); //// pMVCT input
      char CTname[128];
      sprintf(CTname, "%s", projData->reconDataFilename);
      if( (fp=fopen(CTname,"r")) == NULL )
        {
	  fprintf(stderr,"prior input CT file: %s not found\n",CTname);
	  return 0;
        }
      readSize = fread( projData->FBP_reconImageshort, sizeof(short), pMVCTsize * pMVCTsize * projData->FBP_slice, fp );
      if( readSize != pMVCTsize * pMVCTsize * projData->FBP_slice )
	{
	  fprintf(stderr,"WARNING: FBP data size unmatch\n");
	}
      fclose(fp);
    }
  delete1DArray((void*)Preinput);
  return loadProjectionData( projData );
}



PROJECTION_DATA*
newProjectionData()
{
  PROJECTION_DATA* ret = (PROJECTION_DATA*)malloc(sizeof(PROJECTION_DATA));

  if( ret == NULL )
    {
      fprintf(stderr,"ERROR: not enough memory for newProjectionData()\n");
      return NULL;
    }

  ret->anglesRad = NULL;
  ret->xshift = NULL;
  ret->yshift = NULL;
  
  ret->projVolume = NULL;
  ret->mask_projVolume = NULL;
  ret->FBP_reconImageshort = NULL;
  ret->Pri_reconImageshort = NULL;
  
  return ret;
}

void
deleteProjectionData(PROJECTION_DATA* projData)
{
  if( projData == NULL )	return;
  
  delete1DArray( (void*)projData->anglesRad );
  delete1DArray( (void*)projData->xshift );
  delete1DArray( (void*)projData->yshift );
  
  int datatype;
  projData->usedProjNumber = projData->totalProjNumber;
  
  delete1DArray( (void*)projData->projVolume );
  delete1DArray( (void*)projData->mask_projVolume );
  delete1DArray( (void*)projData->FBP_reconImageshort );
  delete1DArray( (void*)projData->Pri_reconImageshort );
  free(projData);
}



// load angle, shift, pixels
int
loadProjectionData(PROJECTION_DATA* projData)
{
  int readSize;

  FILE*	 fp;
  FILE*  fp1;
  int	 tgtCount=0;
  double angle, x=0, y=0;
  int    nskip = projData->sino_fact;                     ////// 2017.11.11
  projData->projImgHeight = projData->projImgHeight/nskip;////// 2017.11.11

  projData->usedProjNumber = projData->totalProjNumber;

    
  projData->anglesRad = (double*)new1DArray(projData->projImgHeight,UNIT_FLOAT64);
  projData->xshift = (double*)new1DArray(projData->projImgHeight,UNIT_FLOAT64);
  projData->yshift = (double*)new1DArray(projData->projImgHeight,UNIT_FLOAT64);
  projData->projVolume = (float*)new1DArray( projData->projImgWidth*
					     projData->projImgHeight*projData->usedProjNumber, projData->datatype );
  projData->mask_projVolume = (unsigned char*)new1DArray( projData->projImgWidth*
                                              projData->projImgHeight*projData->usedProjNumber, UNIT_UINT8 );
  float* projection = (float*)new1DArray(NumberOfChannels*AirScn_Sinogram_Height,projData->datatype); // 2017.12.21
    
  float* airprojection;
  float* sum_airprojection;
  if(strcmp(projData->InPrClassString,"PR") == 0)
    {
      fprintf(stderr,"Pre-treatment CT \n");
      airprojection = (float*)new1DArray(NumberOfChannels*AirScn_Sinogram_Height,projData->datatype); // 2017.12.21
      sum_airprojection = (float*)new1DArray(NumberOfChannels*AirScn_Sinogram_Height,projData->datatype); // 2017.12.21
    }
  else
    {
      fprintf(stderr,"In-treatment CT \n");
      airprojection = (float*)new1DArray(NumberOfChannels*AirScn_Sinogram_Height*projData->usedProjNumber,projData->datatype); // 2017.12.21
      sum_airprojection = (float*)new1DArray(NumberOfChannels*AirScn_Sinogram_Height*projData->usedProjNumber,projData->datatype); // 2017.12.21
    }
  
  for(int i = 0; i < projData->projImgHeight; i++)
    {
      angle = (double)(360.0/projData->projImgHeight*i);
      projData->anglesRad[i] = (projData->StartViewAngle + angle) * PI / 180.0; // gantry angle
      projData->xshift[i] = x;  // = 0 for TomoTherapy
      projData->yshift[i] = y;  // = 0 for TomoTherapy
    }
  if( (fp=fopen(projData->projImgFileName,"rb")) == NULL )
    {
      fprintf(stderr,"data file: %s not found\n",projData->projImgFileName);
      return -1;
    }
  if( (fp1=fopen(projData->AirScanImgFileName,"rb")) == NULL )
    {
      fprintf(stderr,"data file: %s not found\n",projData->AirScanImgFileName);
      return -1;
    }

  ////////   Air Scan INPUT    ////////
  //int length = projData->projImgWidth * projData->projImgHeight * projData->datasize; //// 2017.11.11
  int length = NumberOfChannels * AirScn_Sinogram_Height * projData->datasize; //// 2017.12.21
  if(strcmp(projData->InPrClassString,"PR") == 0)
    {
      readSize = fread( airprojection, 1, length, fp1 );
      for(int i = 0; i < (AirScn_Sinogram_Number); i++)
        //for(int i = 0; i < 2; i++)
        {
	  readSize = fread( airprojection, 1, length, fp1 );
	  for(int k = 0; k < NumberOfChannels * AirScn_Sinogram_Height; k++)             //// 2017.12.21
	    {
	      sum_airprojection[k] += airprojection[k]/(AirScn_Sinogram_Number);
	    }
        }
    }
  else
    {
      //length = projData->projImgWidth * projData->projImgHeight * projData->usedProjNumber * projData->datasize;
      //fread( sum_airprojection, 1, length, fp1 );
      readSize = fread( airprojection, 1, NumberOfChannels * AirScn_Sinogram_Height * projData->datasize, fp1 );//// 2017.11.11
      //fread( airprojection, 1, projData->projImgWidth * AirScn_Sinogram_Height * projData->datasize, fp1 );
      for(int i = 0; i < (AirScn_Sinogram_Number); i++)
        //for(int i = 0; i < 2; i++)
        {
	  readSize = fread( airprojection, 1, NumberOfChannels * AirScn_Sinogram_Height * projData->datasize, fp1 );
	  int k = 0;
	  for(int k1 = 0; k1 < projData->projImgHeight; k1++)                                //// 2017.11.11
	    for(int k2 = NumberOfSkips; k2 < (NumberOfSkips + projData->projImgWidth); k2++) //// 2017.11.11
	      {
		int kk = k1*NumberOfChannels  + k2;                                          //// 2017.11.11
		sum_airprojection[k] += airprojection[kk]/(AirScn_Sinogram_Number);
		k++;
		//sum_airprojection[k] += airprojection[k]/(2.0);
		//fprintf(stderr, "%f \n", sum_airprojection[k] );
	      }
	  /*            for(int k1 = 0; k1 < projData->projImgWidth; k1++)
			{
			for(int k = 0; k < projData->projImgHeight; k++)
			{
			sum_airprojection[k1] += airprojection[k * projData->projImgWidth + k1]
			                                 /(AirScn_Sinogram_Number*projData->projImgHeight);
			//sum_airprojection[k] += airprojection[k]/(2.0);
			}
			fprintf(stderr, "%f \n", sum_airprojection[k1] );
			}
	  */
        }
    }
  
  double projection_max = 0.0;
  fclose(fp1);
  free(airprojection);
  ////////   Sinogram INPUT    ////////  
  if(strcmp(projData->InPrClassString,"PR") == 0)
    {
      int kkk=0;//// 2017.11.11
      for(tgtCount = 0; tgtCount < projData->totalProjNumber; tgtCount++)
        {
	  
	  int readSize;
	  readSize = fread( projection, 1, length, fp );

	  for(int k1 = 0; k1 < AirScn_Sinogram_Height; k1++)                                //// 2017.12.21
	    for(int k2 = NumberOfSkips; k2 < (NumberOfSkips + projData->projImgWidth); k2++) //// 2017.11.11
	      {
		int kk = k1*NumberOfChannels  + k2;                                          //// 2017.11.11		
                float dproj = projection[kk]/(sum_airprojection[kk]);                         //// 2017.12.21
                //if (dproj < 0.01) {
                //    dproj = 1;
                //}
		if(k1/nskip*nskip == k1)
		  {
		    projData->projVolume[kkk] = dproj; //// 2017.11.11
		    projData->mask_projVolume[kkk] = 1;//// 2017.11.11
		    kkk++;//// 2017.11.11
		  }
	      }
	  double zz = projData->StartCouchPosition+tgtCount*projData->CouchSpeed;
	  fprintf(stderr,"Z scan step: %lf \n", zz);
        }
    }
  else
    {
      length = NumberOfChannels * projData->datasize;//// 2017.11.11
      fprintf(stderr,"%d %d %d %d\n", projData->totalProjNumber, projData->projImgHeight, projData->projImgWidth, length);
      for(tgtCount = 0; tgtCount < projData->totalProjNumber; tgtCount++)
        {
	  for(int j=0;j<projData->projImgHeight;j++)
            {
	      int k=0;//// 2017.11.11
	      readSize = fread( projection, 1, length, fp );
	      for(int i = NumberOfSkips; i < (NumberOfSkips + projData->projImgWidth); i++) //// 2017.11.11		
                {
		  int kk = j*NumberOfChannels  + i;                                          //// 2017.11.11		
		  float dproj = projection[kk]/(sum_airprojection[k]*2);
		  if (projection_max < dproj && i >= Start_Channel && i < End_Channel)
                    {
		      projection_max = dproj;
		      if(projection_max > 1.0)
                        {
                          fprintf(stderr, "Maximum projection: %lf > 1.0... STOP \n", projection_max);
			  exit (0);
                        }
		      //fprintf(stderr, "Maximum projection: %lf %d % d\n", projection_max, i, j);
                    }
		  projData->mask_projVolume[(tgtCount*projData->projImgHeight+j)*projData->projImgWidth + k] = 1;//// 2017.11.11
		  if (dproj < 0.1)
                    {
		      dproj = 0.0; // MLC block region
		      projData->mask_projVolume[(tgtCount*projData->projImgHeight+j)*projData->projImgWidth + k] = 0;//// 2017.11.11
                    }
		  
		  projData->projVolume[(tgtCount*projData->projImgHeight+j)*projData->projImgWidth + k] = dproj;//// 2017.11.11
		  //		  kji++;
		  //fprintf(stderr, "%d %f %f\n", tgtCount, projection[kji], dproj );
		  k++;//// 2017.11.11
                }
	      if(projData->DataSamplingRateOrig >= 200)
                {
		  readSize = fread( projection, 1, length, fp );
		  readSize = fread( projection, 1, length, fp );
		  readSize = fread( projection, 1, length, fp );
		  readSize = fread( projection, 1, length, fp );
                }
	      
            }
        }
    }
  fprintf(stderr, "Maximum projection: %lf \n", projection_max);
  free(projection);
  free(sum_airprojection);
  fclose(fp);
  
  // projection volume write
  char	filename[256];
  sprintf(filename,"projections.raw");
  FILE*	fp10=fopen(filename,"wb");
  fwrite( projData->projVolume,
	  getByteSizeOfUnit(projData->datatype),
	  projData->projImgWidth*projData->projImgHeight*projData->usedProjNumber,
	  fp10);
  fclose(fp10);
  
  return tgtCount;
}


// load pixels
int
loadPixelData( char* filename, int offset, int length, unsigned short* pixels )
{
  FILE*	fp;
  int	readSize;
  //fprintf(stderr,"loading %s\n",filename);
  if( (fp=fopen(filename,"rb")) == NULL )
    {
      fprintf(stderr,"data file: %s not found\n",filename);
      return -1;
    }
  
  //fprintf(stderr,"loading %s\n",filename);
  
  fseek(fp,offset,SEEK_SET);
  readSize = fread( (void*)pixels, 1, length, fp );
  fclose(fp);
  
  if( readSize != length )
    {
      fprintf(stderr,"WARNING: pixel data size unmatch\n");
    }
  
  return  readSize;
}






