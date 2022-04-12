#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h> // for Ubuntu
#ifdef _OPENMP
#include <omp.h>
#endif

#include "physParamsTOMO.h"
#include "mallocMD.h"
#include "lengthproj_TOMO.h"


void
lengthproj_TOMO( int ite, int reconSize, float* npf_float, float* projVolume,
		 double reconScale, 
		 int projImgWidth, int projImgHeight,
		 double* cosTable, double* sinTable,
		 double zstart, 
		 double CouchSpeed,
		 char* MaskData,
		 float* systemMatrixfloat, float* xraytracedetector,
		 int thinout)
{
  /*int OSEM_rate = 1;
    int ite_minus = ite - 1;
    int OSEM_initial_offset = ite_minus % OSEM_rate;
    if(ite == 0){
    printf("OSEM_rate = %d, ", OSEM_rate);
    OSEM_rate = 1;
    OSEM_initial_offset = 0;
    }
    printf("OSEM_initial_offset = %d, ", OSEM_initial_offset);
  */

  double rr; // X-ray path length in Voxels

  int sngrmHeight=projImgHeight;
  double detector_Z_pitch = (double)(CouchSpeed/sngrmHeight);

  //z-coodinate of reconstruction slice 
  double rz = zstart; // zstart = z_position
            
  //start sinogram number: 
  int start_angle_num = (int)(rz/detector_Z_pitch-sngrmHeight/(UsedReconAngle)+0.00000001); 
  
  //end sinogram number: 
  int end_angle_num = start_angle_num + sngrmHeight/(UsedReconAngle)*2;

  // number of rotation in start_angle
  int num_rotation = (int)(start_angle_num/sngrmHeight);

#ifdef _OPENMP
#pragma omp parallel for private(rr)
#endif
  for(int jy = 0;jy < reconSize;jy++){
    //fprintf(stderr,"used threads = %d \n",omp_get_num_threads());
    double ry = reconSize*reconScale*0.5 - (jy+0.5)*reconScale;
    
    for(int jx = 0;jx < reconSize;jx++){
      int  mas = jy*reconSize + jx;
      if( MaskData[mas] != 0 ){
	double rx = -reconSize*reconScale*0.5 + (jx+0.5)*reconScale;
	
	
	float npf = 0.0;
	int j_one = 0;
      int sinogram_thinout = end_angle_num - start_angle_num;
      for(int j_thin=0;j_thin<sinogram_thinout;j_thin++) // sinogram range for one slice;
          {
int jjj = start_angle_num + j_thin;
        {
	    int j_onerotation = jjj - (num_rotation * sngrmHeight);
	    if(j_onerotation >= sngrmHeight) j_onerotation = j_onerotation - sngrmHeight;
	    int masth = jy*reconSize*sngrmHeight + jx*sngrmHeight +j_onerotation;
	    
	    // Source location in cartesian coodinate
	    double sx = DIST_BTWN_SRC_AND_ISOCENT_CM*sinTable[j_onerotation];
	    double sy = DIST_BTWN_SRC_AND_ISOCENT_CM*cosTable[j_onerotation];
	    double sz = rz;
	    //if(jx == 200 && jy == 200) fprintf(stderr, "%d %d %lf \n", jjj, j_onerotation, sinTable[j_onerotation]);
	    
	    double zslope = fabs(rz-sz);
	    double yslope = fabs(ry-sy);
	    double xslope = fabs(rx-sx);
	    
	    
	    if(xslope > yslope){
	      double xx0 = rx - 0.5*reconScale;
	      double yy0 = (xx0-sx)/(rx-sx)*(ry-sy)+sy;
	      double zz0 = (xx0-sx)/(rx-sx)*(rz-sz)+sz;
	      double xx1 = rx + 0.5*reconScale;
	      double yy1 = (xx1-sx)/(rx-sx)*(ry-sy)+sy;
	      double zz1 = (xx1-sx)/(rx-sx)*(rz-sz)+sz;
	      rr = sqrt((xx0-xx1)*(xx0-xx1)+(yy0-yy1)*(yy0-yy1)+(zz0-zz1)*(zz0-zz1));
	    }
	    else{
	      double yy0 = ry - 0.5*reconScale;
	      double xx0 = (yy0-sy)/(ry-sy)*(rx-sx)+sx;
	      double zz0 = (yy0-sy)/(ry-sy)*(rz-sz)+sz;
	      double yy1 = ry + 0.5*reconScale;
	      double xx1 = (yy0-sy)/(ry-sy)*(rx-sx)+sx;
	      double zz1 = (yy0-sy)/(ry-sy)*(rz-sz)+sz;
	      rr = sqrt((xx0-xx1)*(xx0-xx1)+(yy0-yy1)*(yy0-yy1)+(zz0-zz1)*(zz0-zz1));
	    }
	    systemMatrixfloat[masth] = rr;
	    
	    // Detector Source location in cartesian coodinate
	    double dx = DIST_BTWN_DETERCTORSRC_AND_ISOCENT_CM*sinTable[j_onerotation];
	    double dy = DIST_BTWN_DETERCTORSRC_AND_ISOCENT_CM*cosTable[j_onerotation];
	    double dz = rz;
	    
	    double aa = (sy-ry)/(sx-rx);
	    double bb = -sx*(sy-ry)/(sx-rx)+sy;
	    
	    double A = 1.0+aa*aa;
	    double B = 2.0*aa*bb-2.0*aa*dy-2.0*dx;
	    double C = dx*dx+bb*bb + dy*dy -2.0*bb*dy-DETECTOR_RADIUS_CURVATURE_CM*DETECTOR_RADIUS_CURVATURE_CM;
	    
	    double px1 = (-B + sqrt(B*B - 4.0*A*C) ) / (2.0*A);
	    double py1 = aa*px1+bb;
	    
	    double px2 = (-B - sqrt(B*B - 4.0*A*C) ) / (2.0*A);
	    double py2 = aa*px2+bb;
	    
	    double px, py;//detector location
	    if(px1*px1+py1*py1 > px2*px2+py2*py2){
	      px = px2;
	      py = py2;
	    }
	    else{
	      px = px1;
	      py = py1;
	    }
	    
	    //detector center
	    double ox = -DIST_BTWN_ISO_AND_DETECTOR_CM*sinTable[j_onerotation];
	    double oy = -DIST_BTWN_ISO_AND_DETECTOR_CM*cosTable[j_onerotation];
	    
	    double dis_p_to_o = (px-ox)*(px-ox)+(py-oy)*(py-oy);
	    double dis_iso_to_detector = DETECTOR_RADIUS_CURVATURE_CM;
	    double detector_angle = acos(1.0 - dis_p_to_o / (2.0*dis_iso_to_detector*dis_iso_to_detector));
	    
	    
	    double rotated_px = px*cosTable[j_onerotation] - py*sinTable[j_onerotation];
	    double alpha;
	    double delta_angle = Delta_Equal_angle_spacing_alpha;
	    
	    double dj;
	    if(rotated_px > 0.0) dj = isocenter_channel - detector_angle/delta_angle;
	    else dj = isocenter_channel + detector_angle/delta_angle;
	    int int_dj = (int)dj;
	    double delta_dj = dj - (double)int_dj;
	    int int_dj2 = int_dj + 1;
	    //if(int_dj >= End_Channel)int_dj2 = int_dj;
	    xraytracedetector[masth] = (float)dj;

	    if(Start_Channel <= dj && dj < End_Channel - 1){
	      
	      double np_para = ((1.0-delta_dj)*projVolume[j_one*projImgWidth+int_dj] 
				+ delta_dj*projVolume[j_one*projImgWidth+int_dj2]);
	      //printf("projImgWidth:%d, int_dj:%d, Start_Channel:%d\n", projImgWidth, int_dj, Start_Channel);
	      npf += (float)(np_para * rr);
	    }
	    else if(Start_Channel > dj){
	      double np_para = projVolume[j_one*projImgWidth+Start_Channel];
	      npf += (float)(np_para * rr);
	    }
	    else if(dj >= End_Channel - 1){
	      double np_para = projVolume[j_one*projImgWidth+End_Channel-1];
	      npf += (float)(np_para * rr);
	    }
	    //if(ite!=0)printf("start:%d, end:%d, jjj= %d\n",start_angle_num,end_angle_num, jjj);
	    j_one++;
	  } // Thin out
	} // projection angle
	
	int vnum = jx + jy*reconSize;
	npf_float[vnum] = npf;
      }//mask
    }//jx=0
  }//jy=0
  
  return;
}


void
lengthproj_TOMO_2( int ite, int reconSize, float* npf_float, float* projVolume,
		   double reconScale, 
		   int projImgWidth, int projImgHeight,
		   double zstart,
		   double CouchSpeed,
		   char* MaskData,
		   float* systemMatrixfloat, float* xraytracedetector,
		   int thinout)
{
  int thinstep = int((ite-1)/thinout);
  int thinnumb = (ite-1)-thinstep*thinout;
  double rr;

  int sngrmHeight=projImgHeight;
  double detector_Z_pitch = (double)(CouchSpeed/sngrmHeight);

  //z-coodinate of reconstruction slice 
  double rz = zstart;
            
  //start sinogram number: 
  int start_angle_num = (int)(rz/detector_Z_pitch-sngrmHeight/(UsedReconAngle)+0.00000001); 
      
  //end sinogram number: 
  int end_angle_num = start_angle_num + sngrmHeight/(UsedReconAngle)*2;
      
      
  // number of ratation in start_angle
  int num_rotation = (int)(start_angle_num/sngrmHeight);

#ifdef _OPENMP
#pragma omp parallel for private(rr)
#endif
  for(int jy=0;jy<reconSize;jy++){
      
    for(int jx=0;jx<reconSize;jx++){
      int  mas = jy*reconSize + jx;
      if( MaskData[mas] != 0 ){
	      
	float npf = 0.0;
	int j_one = thinnumb;

	int sinogram_thinout = int((end_angle_num - start_angle_num)/thinout);
	for(int j_thin=0;j_thin<sinogram_thinout;j_thin++) // sinogram range for one slice;
	  {
	  int jjj = start_angle_num + thinnumb + j_thin * thinout;
	    int j_onerotation = jjj - num_rotation*sngrmHeight;
	    if(j_onerotation >= sngrmHeight){
	      j_onerotation = j_onerotation - sngrmHeight;
	    }
	    int masth = jy*reconSize*sngrmHeight + jx*sngrmHeight +j_onerotation;
	    
	    rr = systemMatrixfloat[masth];
		      
	    float dj = xraytracedetector[masth];

	    int int_dj = (int)dj;
	    double delta_dj = dj - (double)int_dj;
	    int int_dj2 = int_dj + 1;
	    //if(int_dj >= End_Channel)int_dj2 = int_dj;
		      
	    if(Start_Channel <= dj && dj < End_Channel - 1){
	      double np_para = ((1.0-delta_dj)*projVolume[j_one*projImgWidth+int_dj] 
				+ delta_dj*projVolume[j_one*projImgWidth+int_dj2]);
	      //printf("projImgWidth:%d, int_dj:%d, Start_Channel:%d\n", projImgWidth, int_dj, Start_Channel);
	      npf += (float)(np_para * rr);
	    }
	    else if(Start_Channel > dj){
	      double np_para = projVolume[j_one*projImgWidth+Start_Channel];
	      npf += (float)(np_para * rr);
	    }
	    else if(dj >= End_Channel - 1){
	      double np_para = projVolume[j_one*projImgWidth+End_Channel-1];
	      npf += (float)(np_para * rr);
	    }
	    //if(ite!=0)printf("start:%d, end:%d, jjj= %d\n",start_angle_num,end_angle_num, jjj);
	    j_one = j_one + thinout;
	    //} // Thin out
	} // projection angle
	      
	int vnum = jx + jy*reconSize;
	npf_float[vnum] = npf;

      }//mask
    }//jx=0
  }//jy=0
  
  return;
}


