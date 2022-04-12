#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h> // for Ubuntu

#include "physParamsTOMO.h"
#include "mallocMD.h"
#include "filterSinogramTOMO.h"

double**
prepareConstMap(int projWidth, int ReconstructionFilter )
{
	double**	constMap = (double**)new2DArray( End_Channel-Start_Channel, End_Channel-Start_Channel, UNIT_FLOAT64 );

	if( constMap == NULL )
	{
		fprintf(stderr,"ERROR: not enough memory for prepareConstMap\n");
		return NULL;
	}
    /*
    // WEBB_84
	double theConst = 1.0/(8.0*PI*PI*DETECTOR_PITCH_MM_AT_ISO_CENT * 1024 / projWidth);
	double gamma_x_cutoffFreq = GAMMA * CUTOFF_FREQ;
	double PI_x_cutoffFreq = PI * CUTOFF_FREQ;
	double one_minus_alpha = 1.0 - ALPHA;
	double omega = PI/GAMMA;
	double cos_omega = cos(omega);
	double sin_omega = sin(omega);
    */
    
    double arc_spacing = (Corrected_Delta_Equal_angle_spacing);
    
	for(int m = 0; m < End_Channel-Start_Channel; m++)
	{
		for(int n = 0; n < End_Channel-Start_Channel; n++)
		{
            int MN = (m-n);

			if( m == n )
			{
                if(ReconstructionFilter == 1)
                {
                    //Ram-Lak (Ramachandran and Lakshminarayanan)
                    constMap[m][n] = 0.125/(arc_spacing*arc_spacing);
                }
                else if(ReconstructionFilter == 2)
                {
                    //Shepp and Logan
                    double oddterm = PI*PI*(arc_spacing);
                    constMap[m][n] = 0.5/oddterm*(1.0/tan((0.5)*arc_spacing)-1.0/tan((-0.5)*arc_spacing));
                }
                //fprintf(stderr,"m = %d, n = %d %e\n", m, n, constMap[m][n]);
			}
			else if ( MN % 2 != 0)
			{
					
                if(ReconstructionFilter == 1)
                {
                    //Ram-Lak (Ramachandran and Lakshminarayanan)
                    double oddterm = PI*sin(MN*arc_spacing);
                    constMap[m][n] = -0.5/oddterm/oddterm;
                }
                else if(ReconstructionFilter == 2)
                {
                    //Shepp and Logan
                    double oddterm = PI*PI*arc_spacing;
                    constMap[m][n] = 0.5/oddterm*(1.0/tan((MN+0.5)*arc_spacing)-1.0/tan((MN-0.5)*arc_spacing));
                }
                /*
                // WEBB_84
                double cos_theta = cos(theta);
                double sin_theta = sin(theta);
                double theta = PI_x_cutoffFreq * MN;
				double delta = gamma_x_cutoffFreq * MN;
				double eta = 1.0/(1.0-1.0/delta);
				double xi = 1.0/(1.0+1.0/delta);

                constMap[m][n] = theConst /MN/MN *
							( 2.0*ALPHA*theta*sin_theta
					 			+ 2.0*ALPHA*(cos_theta-1.0)
								+ xi*one_minus_alpha*theta*(sin_theta*cos_omega+cos_theta*sin_omega)
					 			+ eta*one_minus_alpha*theta*(sin_theta*cos_omega-cos_theta*sin_omega)
								+ xi*xi*one_minus_alpha*(cos_theta*cos_omega-sin_theta*sin_omega -1.0)
								+ eta*eta*one_minus_alpha*(cos_theta*cos_omega+sin_theta*sin_omega -1.0)
								);
                 */
			}
            else if ( MN % 2 == 0)
            {
                if(ReconstructionFilter == 1)
                {
                    //Ram-Lak (Ramachandran and Lakshminarayanan)
                    constMap[m][n] = 0.0;
                }
                else if(ReconstructionFilter == 2)
                {
                    //Shepp and Logan
                    double oddterm = PI*PI*arc_spacing;
                    constMap[m][n] = 0.5/oddterm*(1.0/tan((MN+0.5)*arc_spacing)-1.0/tan((MN-0.5)*arc_spacing));
                }
            }
		}
	}

	return constMap;
}

void
prepareCorrectedDetectorInterpolation(int projWidth, double* Corrected_Detector_Weight0,
                                      double* Corrected_Detector_Weight1, double*Corrected_Detector_Weight2)
{
    double dalph = (DETECTOR_PITCH_CM)/DETECTOR_RADIUS_CURVATURE_CM;
    double alph_max = (isocenter_channel-Start_Channel)*DETECTOR_PITCH_CM/(DETECTOR_RADIUS_CURVATURE_CM);

    for(int k = 0;k < projWidth ;k++)
    {
        double alph0 = alph_max-dalph*(k-1);
        double alph1 = alph_max-dalph*(k);
        double alph2 = alph_max-dalph*(k+1);
        double xx0 = sqrt((DETECTOR_RADIUS_CURVATURE_CM*DETECTOR_RADIUS_CURVATURE_CM)
                      +(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM)
                      *(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM)
                      +2.0*DETECTOR_RADIUS_CURVATURE_CM*(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM)*cos(alph0));
        double xx1 = sqrt((DETECTOR_RADIUS_CURVATURE_CM*DETECTOR_RADIUS_CURVATURE_CM)
                      +(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM)
                      *(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM)
                      +2.0*DETECTOR_RADIUS_CURVATURE_CM*(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM)*cos(alph1));
        double xx2 = sqrt((DETECTOR_RADIUS_CURVATURE_CM*DETECTOR_RADIUS_CURVATURE_CM)
                      +(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM)
                      *(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM)
                      +2.0*DETECTOR_RADIUS_CURVATURE_CM*(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM)*cos(alph2));
        double  beta0 = acos((DETECTOR_RADIUS_CURVATURE_CM*cos(alph0)+(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM))/xx0);
        double  beta1 = acos((DETECTOR_RADIUS_CURVATURE_CM*cos(alph1)+(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM))/xx1);
        double  beta2 = acos((DETECTOR_RADIUS_CURVATURE_CM*cos(alph2)+(DIST_BTWN_SRC_AND_DETECTOR_CM-DETECTOR_RADIUS_CURVATURE_CM))/xx2);
        if(alph0 < 0) beta0 = -beta0;
        if(alph1 < 0) beta1 = -beta1;
        if(alph2 < 0) beta2 = -beta2;
        double arc_length = DIST_BTWN_SRC_AND_DETECTOR_CM*(Max_beta-Corrected_Delta_Equal_angle_spacing*k);
        double arc_length0 = DIST_BTWN_SRC_AND_DETECTOR_CM*(beta0);
        double arc_length1 = DIST_BTWN_SRC_AND_DETECTOR_CM*(beta1);
        double arc_length2 = DIST_BTWN_SRC_AND_DETECTOR_CM*(beta2);
        if( arc_length < arc_length1)
        {
            Corrected_Detector_Weight0[k] = (arc_length1 - arc_length)/fabs(arc_length1-arc_length0);
            Corrected_Detector_Weight1[k] = Corrected_Detector_Weight0[k];
            Corrected_Detector_Weight2[k] = 0.0;
        }
        else
        {
            Corrected_Detector_Weight0[k] = (arc_length - arc_length1)/fabs(arc_length2-arc_length1);
            Corrected_Detector_Weight1[k] = 0.0;
            Corrected_Detector_Weight2[k] = Corrected_Detector_Weight0[k];
        }
        //fprintf(stderr,"k = %d %lf %lf %lf \n", k,
        //                 Corrected_Detector_Weight0[k] );
    }
}

void
filterSinogram( int projHeight, int zsliceHeight_min, int zsliceHeight_max,
			float** projVolume, float** fltdVoxData,
			double** constMap,
            double* Corrected_Detector_Weight0, double* Corrected_Detector_Weight1, double* Corrected_Detector_Weight2,
            double* xshift, double* yshift,
            char* InPrClassString
            )
{
	if( fltdVoxData == NULL )
	{
		fprintf(stderr,"ERROR: fltdVoxData is NULL\n");
		return;
	}
    //WEB84
	//double omega = PI/GAMMA;
	//double coeff = CUTOFF_FREQ*CUTOFF_FREQ/4.0/
	//  (DETECTOR_PITCH_MM_AT_ISO_CENT * 1024 / projWidth) *(ALPHA*0.5+(1.0-ALPHA)*(1.0/omega*sin(omega)+1.0/omega/omega*(cos(omega)-1.0)));
    
    // Ideal angle destance for "Equal angle spacing algorithm"
    for(int th = zsliceHeight_min; th <= zsliceHeight_max; th++)
	{
        int centWidth = (End_Channel-Start_Channel)/2;
        for(int n = 0; n < projHeight; n++)
        {
          for(int m = 0; m < End_Channel-Start_Channel; m++)
          {
              //Equal Angle Spacing
              for(int k = 0;k < End_Channel-Start_Channel;k++) // p-integral (gamma-integral)
              {
                  double  coeff = DIST_BTWN_SRC_AND_ISOCENT_CM * cos((k-centWidth) * (Corrected_Delta_Equal_angle_spacing));
                  double  projdat = (1.0 - Corrected_Detector_Weight0[k])*projVolume[th*projHeight+n][k+Start_Channel]
                                + Corrected_Detector_Weight1[k]*projVolume[th*projHeight+n][k-1+Start_Channel]
                                + Corrected_Detector_Weight2[k]*projVolume[th*projHeight+n][k+1+Start_Channel];
                  fltdVoxData[th*projHeight+n][m] += constMap[m][k] * (-log(projdat) * coeff);
                  //if(m==1) fprintf(stderr,"m = %d, k = %d %lf %lf %lf %lf\n", m, k,
                  //              rr, rr1, Corrected_Detector_Weight0[k], Corrected_Detector_Weight1[k] );
              }
              if((strcmp(InPrClassString,"IN") == 0) &&
                 (fltdVoxData[th*projHeight+n][m] < 1000.0 || fltdVoxData[th*projHeight+n][m] > 300000.0))
              {
                  fltdVoxData[th*projHeight+n][m] = 0.0;
              }
          }
		//fprintf(stderr,"th = %d, n = %d \n", th, n);

      }
    }
    int zmid = (zsliceHeight_min+zsliceHeight_max)/2;
	char	filename[256];
    sprintf(filename,"fltsinogram.float32.raw");
    FILE*	fp=fopen(filename,"wb");
    for(int n = 0; n < projHeight; n++)
        for(int m = 0; m < End_Channel-Start_Channel; m++)
        {
            float flt=fltdVoxData[zmid*projHeight+n][m];
                fwrite( &flt, getByteSizeOfUnit(UNIT_FLOAT32),1,fp);
        }
    fclose(fp);
    fprintf(stderr,"Write fltdVoxData \n");
	
}


