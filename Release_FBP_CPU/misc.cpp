#include <stdio.h>
#include <math.h>
#include <time.h>

#include "mallocMD.h"
#include "physParamsTOMO.h"


char**
prepareMaskData(int reconSize)
{
	char**	maskPixData = (char**)new2DArray( reconSize, reconSize, UNIT_UINT8 );

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
		if( ((double)j-maskRadius)*((double)j-maskRadius) +
				((double)i-maskRadius)*((double)i-maskRadius) <= sqMaskRad )
		{
			maskPixData[j][i] = 1;
		}
		else
		{
			maskPixData[j][i] = 0;
		}
	}
	return maskPixData;
}

void
clearReconstructedImage(int reconSize, int zslices, float*** reconVoxData, short*** reconVoxDataShort)
{
	int length = reconSize*reconSize*zslices;
	float* pixPtr = (float*)get1DArrayOf3DArray((void***)reconVoxData,UNIT_FLOAT32);
	int length2 = reconSize*reconSize*zslices;
       	short* pixPtr2 = (short*)get1DArrayOf3DArray((void***)reconVoxDataShort,UNIT_SINT16);

	if( pixPtr == NULL )
	{
		fprintf(stderr,"ERROR: reconVoxDatas are NULL\n");
		return;
	}

	for(int i=0;i<length;i++) *pixPtr++ = 0.0;
       	for(int i=0;i<length2;i++) *pixPtr2++ = 0;

}

void
clearSinogramImage(int projWidth, int projHeight_zslicesHeight, float** fltdSinogram)
{
	int length = projWidth*projHeight_zslicesHeight;
	float* pixPtr = (float*)get1DArrayOf2DArray((void**)fltdSinogram,UNIT_FLOAT32);

	if( pixPtr == NULL )
	{
		fprintf(stderr,"ERROR: fltdSinogram is NULL\n");
		return;
	}

	for(int i=0;i<length;i++) *pixPtr++ = 0.0;

}


void
shiftAngles(int numAngles,double* angles,double dAngle)
{
	if( angles == NULL )
	{
		fprintf(stderr,"WARNING: angles is NULL\n");
		return;
	}
	for(int i=0;i<numAngles;i++)	*angles++ += dAngle;
}

void
prepareAngles(int numAngles,double* angles,int reconSize,int** intStartAngle,int** intEndAngle)
{
	if( intStartAngle == NULL || intEndAngle == NULL )
	{
		fprintf(stderr,"ERROR: not enough memory for inteangle\n");
		return;
	}
	double x,y,aangle,tangle,gapv;
	int gapangle0,gapangle1,angle180,angle360;
	for(int j=0;j<reconSize;j++)
	for(int i=0;i<reconSize;i++)
	     {
	       y=reconSize/2-j-0.5;
	       x=-reconSize/2+i+0.5;
	       aangle = atan(fabs(x/y));
	       tangle = atan(fabs(y/x));
	       gapv = 5.0*PI/180.0;
	       if(tangle < gapv ) aangle = 0.5*PI - 5.0*PI/180.0;
	       //if(x>0 && y>0)fprintf(stderr,"%lf, %lf, %lf\n",x,y,aangle*180/PI);
	       //	       intStartAngle[j][i] = 0;
	       for(int ii=1;ii<numAngles;ii++)
	       {
		 if( x > 0.0 && y > 0.0 && PI-aangle > angles[ii-1] && PI-aangle < angles[ii] ) 
		   {
		     intStartAngle[j][i] = ii-1;
		     //if(j==100) fprintf(stderr,"%d, %d, %lf, S%d\n",i,j,aangle*180/PI, ii);
		   }
		 else if( x > 0.0 && y > 0.0 && -aangle > angles[ii-1] && -aangle < angles[ii] ) 
		   {
		     intEndAngle[j][i] = ii;
		     //if(j==100) fprintf(stderr,"%d, %d, %lf, E%d\n",i,j,aangle*180/PI, ii);
		   }
		 else if( x < 0.0 && y > 0.0 && aangle-PI > angles[ii-1] && aangle-PI < angles[ii] ) 
		   {
		     intStartAngle[j][i] = ii-1;
		     //fprintf(stderr,"%d, %d, %lf, S%d\n",i,j,aangle*180/PI, ii);
		   }
		 else if( x < 0.0 && y > 0.0 && aangle > angles[ii-1] && aangle < angles[ii] ) 
		   {
		     intEndAngle[j][i] = ii;
		     //fprintf(stderr,"%d, %d, %lf, E%d\n",i,j,aangle*180/PI, ii);
		   }
		 else if( x < 0.0 && y < 0.0 && -aangle > angles[ii-1] && -aangle < angles[ii] ) 
		 //else if( x < 0.0 && y < 0.0 && -PI+aangle > angles[ii-1] && -PI+aangle < angles[ii] ) 
		   {
		     intStartAngle[j][i] = ii-1;
		     //fprintf(stderr,"%d, %d, %lf, %d\n",i,j,aangle*180/PI, ii);
		   }
		 else if( x < 0.0 && y < 0.0 && PI-aangle > angles[ii-1] && PI-aangle < angles[ii] ) 
		 //else if( x < 0.0 && y < 0.0 && aangle > angles[ii-1] && aangle < angles[ii] ) 
		   {
		     intEndAngle[j][i] = ii;
		     //fprintf(stderr,"%d, %d, %lf, %d\n",i,j,aangle*180/PI, ii);
		   }
		 else if( x > 0.0 && y < 0.0 && aangle > angles[ii-1] && aangle < angles[ii] )
		   {
		     intStartAngle[j][i] = ii-1;
		     //fprintf(stderr,"%d, %d, %lf, S%d\n",i,j,aangle*180/PI, ii);
		   }
		 else if( x > 0.0 && y < 0.0 && aangle-PI > angles[ii-1] && aangle-PI < angles[ii] )
		   {
		     intEndAngle[j][i] = ii;
		     //fprintf(stderr,"%d, %d, %lf, E%d\n",i,j,aangle*180/PI, ii);
		   }
		 if(fabs(angles[ii-1]-angles[ii]) > PI/2) gapangle0 = ii;
		 if(angles[ii] > 0.0 && angles[ii-1] < 0.0 ) gapangle1 = ii;
	       }
	       if(fabs(x) < 0.2*fabs(y) && y > 0 ) 
		   {
		     intStartAngle[j][i] = gapangle0-1;
		     intEndAngle[j][i] = gapangle1;
		   }
	       if(fabs(x) < 0.2*fabs(y) && y < 0 ) 
		   {
		     intStartAngle[j][i] = gapangle1-1;
		     intEndAngle[j][i] = gapangle0;
		   }
	       //intStartAngle[j][i] = gapangle0-1;
	       //intEndAngle[j][i] = gapangle1;
	       /*if(x*x+y*y < 25 ) 
		   {
		     intStartAngle[j][i] = 0;//gapangle0-1;
		     intEndAngle[j][i] = numAngles;//gapangle1;
		   }*/
            }

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























