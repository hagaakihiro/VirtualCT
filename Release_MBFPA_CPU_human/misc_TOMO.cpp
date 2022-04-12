#include <stdio.h>
#include <math.h>
#include <time.h>
#include <malloc.h> // for Ubuntu

#include "mallocMD.h"
#include "physParamsTOMO.h"


char*
prepareMaskData(int reconSize)
{
	char*	maskPixData = (char*)new1DArray( reconSize*reconSize, UNIT_UINT8 );

	if( maskPixData == NULL )
	{
		fprintf(stderr,"ERROR: not enough memory for renderMaskData\n");
		return NULL;
	}

	double maskRadius = (double)(reconSize-1)/2.0;
	double	sqMaskRad = maskRadius*maskRadius;

	for(int j=0;j<reconSize;j++)
	for(int i=0;i<reconSize;i++)
	{
	  int ji = j*reconSize + i;
		if( ((double)j-maskRadius)*((double)j-maskRadius) +
				((double)i-maskRadius)*((double)i-maskRadius) <= sqMaskRad )
		{
			maskPixData[ji] = 1;
		}
		else
		{
			maskPixData[ji] = 0;
		}
	}
	return maskPixData;
}

void
clearReconstructedImage(int recon2Dsize, int recon3Dsize, short* reconVoxData2D, short* reconVoxData3D)
{
	if( reconVoxData2D == NULL || reconVoxData3D == NULL )
	{
		fprintf(stderr,"ERROR: reconVoxDatas are NULL\n");
		return;
	}

	for(int i=0;i<recon2Dsize;i++) *reconVoxData2D++ = 0.0;
       	for(int i=0;i<recon3Dsize;i++) *reconVoxData3D++ = 0.0;

}

void
clear_npf(int recon2Dsize, float* npf)
{
  if( npf == NULL )
    {
      fprintf(stderr,"ERROR: npf is NULL\n");
      return;
    }
  
  for(int i=0;i<recon2Dsize;i++) *npf++ = 0.0;
}

void
clearSinogramImage(int projsize, float* sinogram)
{
  if( sinogram == NULL )
    {
      fprintf(stderr,"ERROR: Sinogram is NULL\n");
      return;
    }
  
  for(int i=0;i<projsize;i++) *sinogram++ = 0.0;

}


int
prepareTrigonometricTable(int numAngles,double* anglesRad,double** cosTable,double** sinTable)
{
	*cosTable = (double*)new1DArray( numAngles, UNIT_FLOAT64 );
	*sinTable = (double*)new1DArray( numAngles, UNIT_FLOAT64 );

	if( *cosTable == NULL || *sinTable == NULL )
	{
		fprintf(stderr,"ERROR: not enough memory for renderMaskData\n");
		return 0;
	}
	
	for(int i=0;i<numAngles;i++)
	{
		(*cosTable)[i] = cos(anglesRad[i]);
		(*sinTable)[i] = sin(anglesRad[i]);
	}

	return numAngles;
}

static clock_t start, stop;


void
startTimer()
{
	start = clock();
}

void
stopTimer()
{
	stop = clock();
}

void
showLapTime(char* label)
{
	fprintf(stderr,"%s : %3.2lf sec.\n",label,((double)stop-(double)start)/(double)CLOCKS_PER_SEC);
}























