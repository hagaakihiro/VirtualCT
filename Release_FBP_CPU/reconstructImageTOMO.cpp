#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "physParamsTOMO.h"
#include "mallocMD.h"
#include "reconstructImageTOMO.h"


void
reconstructImage( int reconSize, float*** reconVoxData, char** maskData,
			int sngrmWidth, int sngrmHeight, float** sngrmVoxData,
			double* cosTable, double* sinTable, 
			double reconOffset, double reconStep, int reconSlices,
			double StartCouchPosition, double CouchSpeed, int usedProjNumber,
			double xfactor, double scale,
            char* InPrClassString
            )
{
    double centerPos = (double)reconSize/2.0 - 0.5;
    double detector_Z_pitch = (double)(CouchSpeed/sngrmHeight);

	if( reconVoxData == NULL )
	{
		fprintf(stderr,"ERROR: reconVoxData is NULL\n");
		return;
	}

    int usesinogram = int(sngrmHeight/(UsedReconAngle)*2);
	
	for(int z = 0;z < reconSlices;z++)
	{
		double zz = reconOffset + z * reconStep; // reconstruction z coordinate//スライス中心のz座標
		//int z_image = (int)((zz-StartCouchPosition)/detector_Z_pitch-sngrmHeight/(UsedReconAngle)+0.00001);
		int z_image = 0;
        //full sinogram height//最初のサイノ番号
		//int z_image = (int)((zz-StartCouchPosition)/detector_Z_pitch - reconStep/(2.0*detector_Z_pitch));//Limited angle
		//int iz_image = (int)(z_image/sngrmHeight)*sngrmHeight; //Num. of rotation
		int iz_image = 0;
        int thdb0 = z_image - iz_image;//ガントリー回転目のなかでのサイノグラム番号

        fprintf(stderr,"zz-c_position = %lf, z_image = %d %d %d \n",zz-StartCouchPosition,z_image,
                iz_image, thdb0);
        
        for(int y = 0;y < reconSize;y++)
        {
            for(int x = 0;x < reconSize;x++)
            {
                if( maskData[y][x] != 0 )
                {
                    double xx = -(x - centerPos) * scale; // phase; CW rotation!
                    double yy = -(y - centerPos) * scale;
				
                    double count = 0.0;
                    //for(int th = 0; th < sngrmHeight; th++)
                    for(int th = 0; th < usesinogram; th++)
                    //for(int th = 0; th < (int)(reconStep/detector_Z_pitch); th++)//Limited angle
                    {
                        double dbeta = 2.0*PI/(sngrmHeight);
                        int th1 = th + thdb0;
                        if(th1 >= sngrmHeight) th1 = th1 - sngrmHeight;
                        else if (th1 < 0) th1 = th1 + sngrmHeight;
                        double	cos_thdb1 = cosTable[th1];
                        double	sin_thdb1 = sinTable[th1];
                        
                        double X =  xx*cos_thdb1 + yy*sin_thdb1;
                        double Y = -xx*sin_thdb1 + yy*cos_thdb1;
                        // inversed by Y-coodinate; The image domain usually uses it.
                        double L2 = (-Y+DIST_BTWN_SRC_AND_ISOCENT_CM)*(-Y+DIST_BTWN_SRC_AND_ISOCENT_CM)+X*X;
                        double fGAM = atan(X/(DIST_BTWN_SRC_AND_ISOCENT_CM - Y))/(Corrected_Delta_Equal_angle_spacing)
                                    + (End_Channel-Start_Channel)/2;
                        int GAM = (int)(fGAM);
                        double w_GAM = fGAM - GAM;
                        //if(x==100 && y==100) fprintf(stderr,"%d %lf %lf %d %lf\n", th1, sin_thdb1, cos_thdb1, GAM, fGAM);
                        if( GAM >= 0 && GAM < End_Channel-Start_Channel -1 )
                        {
                            int iphase = 0;
                            if(th < usesinogram/2) iphase = 1;
                            else if(th > usesinogram/2) iphase = -1;
                            if((z_image+iphase*sngrmHeight)+th < sngrmHeight) iphase = 0;
                            if((z_image+iphase*sngrmHeight)+th > (usedProjNumber-1)*sngrmHeight) iphase = 0;

                            double sum1 = (1-w_GAM)*sngrmVoxData[z_image+th][GAM]
                                            + w_GAM*sngrmVoxData[z_image+th][GAM+1];
                            //double sum2 = (1-w_GAM)*sngrmVoxData[(z_image+iphase*sngrmHeight)+th][GAM]
                            //                + w_GAM*sngrmVoxData[(z_image+iphase*sngrmHeight)+th][GAM+1];
                            //double z_weight = 0.5*fabs(th-(sngrmHeight/2.0))/(sngrmHeight/2.0);
                            double sum2 = 0;
                            double z_weight = 0.0;
                        
                            double sum = (1.0-z_weight)*sum1 + z_weight*sum2;
                            reconVoxData[z][y][x] += (float)( 1.0/L2 * sum * dbeta * xfactor);//使用するアングルで係数を変える
                            //printf("L2=%lf, sum=%lf, dbeta=%lf, xfactor=%lf\n", L2, sum, dbeta, xfactor);
                        }
                        else
                        {
                            //fprintf(stderr, "GAM = %d not within the detector array \n", GAM);
                        }
                    } // th
                } // mask
            } // x
        } // y
    } // z
}
