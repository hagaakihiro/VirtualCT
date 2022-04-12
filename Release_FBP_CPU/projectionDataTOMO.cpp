#include <stdio.h>
#include <string.h>
#include <stdlib.h> // for MAC
//#include <malloc.h>

#include "physParamsTOMO.h"
#include "mallocMD.h"
#include "projectionDataTOMO.h"


#define LINEBUF_SIZE	512

static char	lineBuf[LINEBUF_SIZE];

//////

static void
skipLines(FILE* fp,int nLines)
{
	for(int i=0;i<nLines;i++)	fgets(lineBuf,LINEBUF_SIZE,fp);
}

static void
loadInteger(FILE* fp,int* value)
{
	skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%d",value);
}

static void
loadfloat(FILE* fp,float* value)
{
	skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%f",value);
}

static void
loaddouble(FILE* fp,double* value)
{
    skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%lf",value);
}

static void
loadIntegerHex(FILE* fp,int* value)
{
	skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%X",value);
}

static void
loadString(FILE* fp,char* string)
{
	skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%s",string);
}

//////

int
loadData(char* infoFilename,PROJECTION_DATA* projData)
{
	FILE*	fp;

	if( (fp=fopen(infoFilename,"r")) == NULL )
	{
		fprintf(stderr,"info file: $s not found\n",infoFilename);
		return 0;
	}

	// skip comments (6lines)
	skipLines(fp,6);
	
	// load each element
	loadString( fp, projData->InPrClassString );		// In-treatment/Pre-treatment label
    loadInteger( fp, &projData->ReconstructionFilter );	// Reconstruction Filter; Ram-Lak for 1, Shepp-Logan for 2
	loadInteger( fp, &projData->projImgWidth );         // sinogram image width
	loadInteger( fp, &projData->projImgHeight );		// sinogram image height
    loadString( fp, projData->projImgDataType );        // Data type
	loadInteger( fp, &projData->projImgFileOffset );	// projection image offset
    loadString( fp, projData->projImgFileName );		// projection image file name
    loadString( fp, projData->AirScanImgFileName );		// Air scan image file name
	loaddouble( fp, &projData->GantryRotationPeriod );	// Gantry rotation period [s]
    loaddouble( fp, &projData->DataSamplingRate );       // Data Sampling Rate [Hz]
    loaddouble( fp, &projData->StartViewAngle );         // Starting View Angle [degree]
    loaddouble( fp, &projData->CouchSpeed );             // Couch Speed [cm/rotation]
    loaddouble( fp, &projData->StartCouchPosition );     // Start couch position [cm]
	loadInteger( fp, &projData->reconSize );			// reconstruction size (square only)
	loaddouble( fp, &projData->reconPixSize );			// reconstruction pixel size in cm
	loadString( fp, projData->reconTgtClassString );	// reconstruction target class label
	loaddouble( fp, &projData->reconOffset );         	// reconstruction start z position in cm
	loaddouble( fp, &projData->reconStep );				// reconstruction z step in cm
	loadInteger( fp, &projData->reconSlices );          // number of reconstructed slices
	loadInteger( fp, &projData->xfactor );		        // factor for CT value
	loadString( fp, projData->reconDataFilename );		// reconstruction file name

	fclose(fp);

    projData->totalProjNumber = (int)(projData->projImgHeight/
                                      (projData->GantryRotationPeriod * projData->DataSamplingRate));
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
    
	delete2DArray( (void**)projData->projVolume,
				projData->projImgWidth, projData->projImgHeight*projData->usedProjNumber,
				projData->datatype);

	free(projData);
}



// load angle, shift, pixels
int
loadProjectionData(PROJECTION_DATA* projData)
{
	FILE*	fp;
    FILE*   fp1;
	int	tgtCount=0;
	double	angle, x=0, y=0;

    projData->usedProjNumber = projData->totalProjNumber;

    
    projData->anglesRad = (double*)new1DArray(projData->projImgHeight,UNIT_FLOAT64);
	projData->xshift = (double*)new1DArray(projData->projImgHeight,UNIT_FLOAT64);
	projData->yshift = (double*)new1DArray(projData->projImgHeight,UNIT_FLOAT64);
	projData->projVolume = (float**)new2DArray( projData->projImgWidth,
                    projData->projImgHeight*projData->usedProjNumber, projData->datatype );
    float* projection = (float*)new1DArray(projData->projImgWidth*projData->projImgHeight,projData->datatype);

    /*
    float* airprojection;
    float* sum_airprojection;
    if(strcmp(projData->InPrClassString,"PR") == 0)
    {
        fprintf(stderr,"Pre-treatment CT \n");
        airprojection = (float*)new1DArray(projData->projImgWidth*projData->projImgHeight,projData->datatype);
        sum_airprojection = (float*)new1DArray(projData->projImgWidth*(projData->projImgHeight),projData->datatype);
    }
    else
    {
        fprintf(stderr,"In-treatment CT \n");
        airprojection = (float*)new1DArray(projData->projImgWidth*projData->projImgHeight*projData->usedProjNumber,projData->datatype);
        sum_airprojection = (float*)new1DArray(projData->projImgWidth*projData->projImgHeight*projData->usedProjNumber,projData->datatype);
    }
    */
    for(int i = 0; i < projData->projImgHeight; i++)
    {
        angle = (double)(360.0/projData->projImgHeight*i);
        projData->anglesRad[i] = (projData->StartViewAngle + angle) * PI / 180.0;
        projData->xshift[i] = x;
        projData->yshift[i] = y;
        //fprintf(stderr,"ERROR: %d %d \n", projData->projImgWidth, projData->projImgHeight);
    }
    if( (fp=fopen(projData->projImgFileName,"rb")) == NULL )
    {
        fprintf(stderr,"data file: %s not found\n",projData->projImgFileName);
        return -1;
    }
    //if( (fp1=fopen(projData->AirScanImgFileName,"rb")) == NULL )
    //{
    //    fprintf(stderr,"airscan data file: %s not found\n",projData->AirScanImgFileName);
    //    return -1;
    //}
    
    //projData->projImgFileOffset=3*sizeof(float);
    //fseek(fp1,projData->projImgFileOffset,SEEK_SET);
    int length = projData->projImgWidth * projData->projImgHeight * projData->datasize;
    /*
    if(strcmp(projData->InPrClassString,"PR") == 0)
    {
        fread( airprojection, 1, length, fp1 );
        //for(int i = 0; i < (AirScn_Sinogram_Number); i++)
        for(int i = 0; i < 2; i++)
        {
            fread( airprojection, 1, length, fp1 );
            for(int k = 0; k < projData->projImgWidth*(projData->projImgHeight); k++)
            {
                //sum_airprojection[k] += airprojection[k]/(AirScn_Sinogram_Number);
                sum_airprojection[k] += airprojection[k]/(2.0);
            }
        }
    }
    else
    {
        length = projData->projImgWidth * projData->projImgHeight * projData->usedProjNumber * projData->datasize;
        fread( sum_airprojection, 1, length, fp1 );
    }
 
    fclose(fp1);
    free(airprojection);
    */
    //projData->projImgFileOffset=2*sizeof(float);
    //fseek(fp,projData->projImgFileOffset,SEEK_SET);

    if(strcmp(projData->InPrClassString,"PR") == 0)
    {
        for(tgtCount = 0; tgtCount < projData->totalProjNumber; tgtCount++)
        {

            int readSize;
            readSize = fread( projection, 1, length, fp );

            int kji=0;
            //fprintf(stderr,"ERROR: %d %d\n", projData->projImgHeight,projData->projImgWidth);
            for(int j=0;j<projData->projImgHeight;j++)
            for(int i=0;i<projData->projImgWidth;i++)
            {
	      float dproj = projection[kji]; //sum_airprojection[kji];
                //if (dproj < 0.01) {
                //    dproj = 1;
                //}
                projData->projVolume[tgtCount*projData->projImgHeight+j][i] = dproj;
                kji++;
            }
            double zz = projData->StartCouchPosition+tgtCount*projData->CouchSpeed;
            fprintf(stderr,"Z scan step: %lf \n", zz);
        }
    }
    else
    {
        int readSize;
        readSize = fread( projection, 1, length, fp );
        int kji=0;
        for(tgtCount = 0; tgtCount < projData->totalProjNumber; tgtCount++)
        for(int j=0;j<projData->projImgHeight;j++)
        for(int i=0;i<projData->projImgWidth;i++)
            {
	      float dproj = projection[kji];//sum_airprojection[kji];
                //if (dproj < 0.01 || dproj > 1)
                if (projection[kji] < 50000.0 )
                {
                    dproj = 1;
                }
                projData->projVolume[tgtCount*projData->projImgHeight+j][i] = dproj;
                kji++;
            }
    }
 
    free(projection);
    //free(sum_airprojection);
    fclose(fp);
    /*
	// projection volume
	char	filename[256];
	sprintf(filename,"projections.raw");
	FILE*	fp10=fopen(filename,"wb");
	fwrite( get1DArrayOf2DArray((void**)projData->projVolume,projData->datatype),
		getByteSizeOfUnit(projData->datatype),
		projData->projImgWidth*projData->projImgHeight*projData->usedProjNumber,
		fp10);
	fclose(fp10);
    */
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






