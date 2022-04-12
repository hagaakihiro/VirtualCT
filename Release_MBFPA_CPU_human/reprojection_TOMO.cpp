#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h> // for Ubuntu
#ifdef _OPENMP
#include <omp.h>
#endif

#include "physParamsTOMO.h"
#include "mallocMD.h"
#include "reprojection_TOMO.h"

void
reprojection_TOMO( int ite, int reconSize, float* reconImagefloat, float* ProjImagefloat,
		   double reconScale,
		   int projImgWidth, int projImgHeight,
		   double* cosTable, double* sinTable,
		   double zstart,
		   double CouchSpeed,
		   double xfactor, double* xalpha, char* MaskData,
		   float* reprojection_systemMatrixfloat,
		   int* reprojection_vnum, int* reprojection_sum,
		   int thinout, float* Attenuation, float* pdf_pe, int num_mat, int NE, int ke)
{
  int nxm = reconSize*0.5;
  int nym = reconSize*0.5;
  int pri_size = reconSize*reconSize;
  int sngrmHeight=projImgHeight;
  double detector_Z_pitch = (double)(CouchSpeed/sngrmHeight);
	
  //start sinogram number: OK
  //int start_angle_num = (int)(zstart/detector_Z_pitch-sngrmHeight/(UsedReconAngle)+0.00000001); 
  int start_angle_num = 0;
  //end sinogram number: OK
  int end_angle_num = start_angle_num + sngrmHeight/(UsedReconAngle)*2; 
  
  double rmax = reconSize*0.5 * 1.0;
  int sinogram_thinout = int((end_angle_num - start_angle_num)/thinout);//20171023
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int j_thin=0;j_thin<sinogram_thinout;j_thin++) // sinogram range for one slice;
  {
    int jjj = start_angle_num + j_thin * thinout;
      {
      //fprintf(stderr,"used threads = %d \n",omp_get_num_threads());
      
      int num_rotation = (int)(jjj/sngrmHeight);
      int j_onerotation = jjj-num_rotation*sngrmHeight;// OK
      int j_one = jjj - start_angle_num;
      
      // Source location in cartesian coodinate
      double sx = DIST_BTWN_SRC_AND_ISOCENT_CM*sinTable[j_onerotation]; //OK
      double sy = DIST_BTWN_SRC_AND_ISOCENT_CM*cosTable[j_onerotation]; //OK
      double sz = 0.0;
      
      double pz = sz; 
      
      for(int jj=Start_Channel;jj<End_Channel;jj++){	    // Width
	double delta_angle = Delta_Equal_angle_spacing_alpha;
	double theta = (isocenter_channel-jj)*Delta_Equal_angle_spacing_alpha;
	
	int jjjj = j_one*projImgWidth+jj;
	
	ProjImagefloat[jjjj]=0.0;
	
	double ox0 = DETECTOR_RADIUS_CURVATURE_CM*sin(theta);
	double oy0 = -DETECTOR_RADIUS_CURVATURE_CM*cos(theta);
	
	double px0 = ox0;
	double py0 = oy0 + DIST_BTWN_DETERCTORSRC_AND_ISOCENT_CM;	
	
	double px = +px0*cosTable[j_onerotation] + py0*sinTable[j_onerotation];
	double py = -px0*sinTable[j_onerotation] + py0*cosTable[j_onerotation];
	
	double xx[4], yy[4], zz[4];
	int m = 0;
	yy[m] = double(nym)*reconScale;
	if(fabs((yy[m]-sy)/(py-sy)*(px-sx)+sx) < double(nxm)*reconScale){
	  xx[m] = (yy[m]-sy)/(py-sy)*(px-sx)+sx;
	  zz[m] = sz;
	  m++;
	}
	yy[m] = -double(nym)*reconScale;
	if(fabs((yy[m]-sy)/(py-sy)*(px-sx)+sx) < double(nxm)*reconScale){
	  xx[m] = (yy[m]-sy)/(py-sy)*(px-sx)+sx;
	  zz[m] = sz;
	  m++;
	}
	xx[m] = double(nxm)*reconScale;
	if(fabs((xx[m]-sx)/(px-sx)*(py-sy)+sy) < double(nym)*reconScale){
	  yy[m] = (xx[m]-sx)/(px-sx)*(py-sy)+sy;
	  zz[m]= sz;
	  m++;
	}
	xx[m] = -double(nxm)*reconScale;
	if(fabs((xx[m]-sx)/(px-sx)*(py-sy)+sy) < double(nym)*reconScale){
	  yy[m] = (xx[m]-sx)/(px-sx)*(py-sy)+sy;
	  zz[m]=sz;
	  m++;
	}
	double xi, yi, zi, xf, yf, zf; // Initial and Final Voxels, in which X-ray passes through
	if((xx[0]-sx)*(xx[0]-sx)+(yy[0]-sy)*(yy[0]-sy)+(zz[0]-sz)*(zz[0]-sz) < (xx[1]-sx)*(xx[1]-sx)+(yy[1]-sy)*(yy[1]-sy)+(zz[1]-sz)*(zz[1]-sz)){
	  xi=xx[0]; yi=yy[0]; zi=zz[0]; xf=xx[1]; yf=yy[1]; zf=zz[1];
	}
	else{
	  xi=xx[1]; yi=yy[1]; zi=zz[1]; xf=xx[0]; yf=yy[0]; zf=zz[0];
	}
	//fprintf(stderr,"%lf %lf %lf %lf %lf %lf\n", xi, yi, zi, xf, yf, zf);
	double xb0 = xi, yb0 = yi, zb0 = zi, xb1, yb1, zb1;
	double xii, yii, zii, xiii, yiii, ziii, xiiii, yiiii, ziiii, yj, xj, rr, transp, total_length;
	int xma, yma, zma, xmb, ymb, zmb, xmc, ymc, zmc, icheck, vnum;
	if(fabs(xi-xf) > fabs(yi-yf)){ // dx is larger than dy
	  int iks = 0;
	  double pha = (xf-xi)/fabs(xi-xf);
	  int xin, ik=0;
	  double chy, chz;
	  if((xi < 0 && (xf-xi) < 0) || (xi > 0 && (xf-xi) > 0)) xin = int(xi) + (xf-xi)/fabs(xf-xi);
	  else xin = int(xi);
	  xii = xin;
	  yii = (xii-sx)/(px-sx)*(py-sy)+sy;
	  //zii = (xii-sx)/(px-sx)*(pz-sz)+sz;
	  transp = 0.0;
	  total_length = 0.0;
	  while( xii <= double(nxm)*reconScale && xii >= -double(nxm)*reconScale &&  yii <= double(nym)*reconScale && yii >= -double(nym)*reconScale){
	    xii = xin + pha * ik*reconScale;
	    ymb=yma;
	    yj = yii;
	    yii = (xii-sx)/(px-sx)*(py-sy)+sy;
	    yma=(int)(yii/reconScale);
	    icheck = 0;
	    chy = 0;
	    chz = 0;
	    if((ik != 0 && yma != ymb) || (ik != 0 && (yii/fabs(yii) != yj/fabs(yj)))){
	      if((yii < 0 && (yf-yi) < 0) || (yii > 0 && (yf-yi) > 0)) yiii = (double)(yma*reconScale);
	      else yiii = (double)(yma*reconScale) - (yf-yi)/fabs(yf-yi)*reconScale;
	      xiii = (yiii-sy)/(py-sy)*(px-sx)+sx;
	      chy = (xii-xiii)*(xii-xiii)+(yii-yiii)*(yii-yiii);
	      icheck++;
	    }
	    if( icheck == 1 && chy != 0){
	      iks++;
	      xb1=xiii;
	      yb1=yiii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1));
	      xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	      if(xmc >= 0 && ymc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 ){
		vnum = xmc + ymc*reconSize;
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		if(reconImagefloat[vnum] < 0 || r >= rmax) reconImagefloat[vnum]=0;
		for(int mat = 0; mat < num_mat; mat++)
		{
		  double ph = 1.0;
		  if(reconImagefloat[mat*pri_size + vnum] < 0 || r >= rmax) ph = 0.0;
		  transp += reconImagefloat[mat*pri_size + vnum] * Attenuation[mat*NE + ke] * rr * ph;
		}
		//transp += reconImagefloat[vnum] * rr;
		total_length += rr;
		int masth = iks*projImgWidth*projImgHeight + j_onerotation*projImgWidth+jj;
		if(masth > projImgWidth*projImgHeight*reconSize*3) fprintf(stderr,"*** %d\n", iks);
		reprojection_systemMatrixfloat[masth] = rr;
		reprojection_vnum[masth] = vnum;
	      }
	      xb0=xb1;
	      yb0=yb1;
	    }
	    iks++;
	    xb1=xii;
	    yb1=yii;
	    rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1));
	    xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	    ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	    if(xmc >= 0 && ymc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 ){
	      vnum = xmc + ymc*reconSize;
	      double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
	      if(reconImagefloat[vnum] < 0 || r >= rmax) reconImagefloat[vnum] = 0;
	      for(int mat = 0; mat < num_mat; mat++)
		{
		  double ph = 1.0;
		  if(reconImagefloat[mat*pri_size + vnum] < 0 || r >= rmax) ph = 0.0;
		  transp += reconImagefloat[mat*pri_size + vnum] * Attenuation[mat*NE + ke] * rr * ph;
		}
	      //transp += reconImagefloat[vnum] * rr;
	      total_length += rr;
	      int masth = iks*projImgWidth*projImgHeight + j_onerotation*projImgWidth+jj;
	      reprojection_systemMatrixfloat[masth] = rr;
	      reprojection_vnum[masth] = vnum;
	    }
	    xb0=xb1;
	    yb0=yb1;
	    ik++;
	  } // while
	  int jone = j_onerotation*projImgWidth+jj;
	  reprojection_sum[jone] = iks;
	  double intens = pdf_pe[ke] * exp(-transp/xfactor); // * photon count n_i/n_0;
	  ProjImagefloat[jjjj] = (float)(intens);
	  //fprintf(stderr,"%d %d \n", jjjj, reprojection_sum[jjjj]);
	  
	} // dx dy
	else{ // dy is larger than dx
	  int iks = 0;
	  double pha = (yf-yi)/fabs(yi-yf);
	  int yin, ik=0;
	  double chy, chz;
	  if((yi < 0 && (yf-yi) < 0) || (yi > 0 && (yf-yi) > 0)) yin = int(yi) + (yf-yi)/fabs(yf-yi);
	  else yin = int(yi);
	  yii = yin;
	  xii = (yii-sy)/(py-sy)*(px-sx)+sx;
	  transp = 0.0;
	  total_length = 0.0;
	  while( xii <= double(nxm)*reconScale && xii >= -double(nxm)*reconScale &&  yii <= double(nym)*reconScale && yii >= -double(nym)*reconScale){
	    yii = yin + pha * ik*reconScale;
	    xmb=xma;
	    xj = xii;
	    xii = (yii-sy)/(py-sy)*(px-sx)+sx;
	    xma=(int)(xii/reconScale);
	    icheck = 0;
	    chy = 0;
	    chz = 0;
	    if((ik != 0 && xma != xmb) || (ik != 0 && (xii/fabs(xii) != xj/fabs(xj)))){
	      if((xii < 0 && (xf-xi) < 0) || (xii > 0 && (xf-xi) > 0)) xiii = (double)(xma*reconScale);
	      else xiii = (double)(xma*reconScale) - (xf-xi)/fabs(xf-xi)*reconScale;
	      yiii = (xiii-sx)/(px-sx)*(py-sy)+sy;
	      chy = (xii-xiii)*(xii-xiii)+(yii-yiii)*(yii-yiii);
	      icheck++;
	    }
	    if( icheck == 1 && chy != 0){
	      iks++;
	      xb1=xiii;
	      yb1=yiii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1));
	      xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	      if(xmc >= 0 && ymc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 ){
		vnum = xmc + ymc*reconSize;
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		if(reconImagefloat[vnum] < 0 || r >= rmax) reconImagefloat[vnum]=0;
		for(int mat = 0; mat < num_mat; mat++)
		{
		  double ph = 1.0;
		  if(reconImagefloat[mat*pri_size + vnum] < 0 || r >= rmax) ph = 0.0;
		  transp += reconImagefloat[mat*pri_size + vnum] * Attenuation[mat*NE + ke] * rr * ph;
		}
		//transp += reconImagefloat[vnum] * rr;
		total_length += rr;
		int masth = iks*projImgWidth*projImgHeight + j_onerotation*projImgWidth+jj;
		reprojection_systemMatrixfloat[masth] = rr;
		reprojection_vnum[masth] = vnum;
	      }
	      xb0=xb1;
	      yb0=yb1;
	    }
	    iks++;
	    xb1=xii;
	    yb1=yii;
	    rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1));
	    xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	    ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	    if(xmc >= 0 && ymc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 ){
	      vnum = xmc + ymc*reconSize;
	      double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
	      if(reconImagefloat[vnum] < 0 || r >= rmax) reconImagefloat[vnum]=0;
	      for(int mat = 0; mat < num_mat; mat++)
		{
		  double ph = 1.0;
		  if(reconImagefloat[mat*pri_size + vnum] < 0 || r >= rmax) ph = 0.0;
		  transp += reconImagefloat[mat*pri_size + vnum] * Attenuation[mat*NE + ke] * rr * ph;
		}
	      //transp += reconImagefloat[vnum] * rr;
	      int masth = iks*projImgWidth*projImgHeight + j_onerotation*projImgWidth+jj;
	      reprojection_systemMatrixfloat[masth] = rr;
	      reprojection_vnum[masth] = vnum;						
	    }
	    xb0=xb1;
	    yb0=yb1;
	    ik++;
	  } //while
	  int jone = j_onerotation*projImgWidth+jj;
	  reprojection_sum[jone] = iks;
	  double intens = pdf_pe[ke] * exp(-transp/xfactor); // * photon count n_i/n_0;
	  ProjImagefloat[jjjj] = (float)(intens);              
	} //dx dy	  
      }
    } // Thinout
  } // jjj 
  
  return;
}

void
reprojection_TOMO_2( int ite, int reconSize, float* reconImagefloat, float* ProjImagefloat,
		     double reconScale,
		     int projImgWidth, int projImgHeight,
		     double* cosTable, double* sinTable,
		     double zstart,
		     double CouchSpeed,
		     double xfactor, double* xalpha, char* MaskData,
		     float* reprojection_systemMatrixfloat,
		     int* reprojection_vnum, int* reprojection_sum,
		     int thinout, float* Attenuation, float* pdf_pe, int num_mat, int NE, int ke)
{
  int pri_size = reconSize*reconSize;
  int thinstep = int((ite-1)/thinout);
  int thinnumb = (ite-1)-thinstep*thinout;
  int sngrmHeight=projImgHeight;
  double detector_Z_pitch = (double)(CouchSpeed/sngrmHeight);

  //start sinogram number: OK
  //int start_angle_num = (int)(zstart/detector_Z_pitch-sngrmHeight/(UsedReconAngle)+0.00000001); 
  int start_angle_num = 0;
  //end sinogram number: OK
  int end_angle_num = (int)(start_angle_num + sngrmHeight/(UsedReconAngle)*2); 
  int sinogram_thinout = int((end_angle_num - start_angle_num)/thinout);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int j_thin=0;j_thin<sinogram_thinout;j_thin++) // sinogram range for one slice;
    {
      int jjj = start_angle_num + thinnumb + j_thin * thinout;
	{
	  int num_rotation = (int)(jjj/sngrmHeight);
	  int j_onerotation = jjj-num_rotation*sngrmHeight;// OK
	  int j_one = jjj - start_angle_num;
	  
	  for(int jj=Start_Channel;jj<End_Channel;jj++) // Width
	    {	  
	      int jjjj = j_one*projImgWidth+jj;
	      int jone = j_onerotation*projImgWidth+jj;
	      //ProjImagefloat[jjjj]=0.0;
	      double transp = 0.0;
	      for( int iks = 1; iks < reprojection_sum[jone]; iks++)
		{
		  int masth = iks*projImgWidth*projImgHeight + j_onerotation*projImgWidth+jj;
		  float rr = reprojection_systemMatrixfloat[masth];
		  int vnum = reprojection_vnum[masth];
		  for(int mat = 0; mat < num_mat; mat++)
		    {
		      double ph = 1.0;
		      //if(reconImagefloat[mat*pri_size + vnum] < 0 || r >= rmax) ph = 0.0;
		      transp += reconImagefloat[mat*pri_size + vnum] * Attenuation[mat*NE + ke] * rr * ph;
		    }
		  //transp += reconImagefloat[vnum] * rr;
		}
	      double intens = pdf_pe[ke] * exp(-transp/xfactor); // * photon count n_i/n_0;
	      ProjImagefloat[jjjj] = (float)(intens);//ProjImagefloat[jjjj] = 50.0;			
	    } // jj width
	} // Thinout
    } // jjj height

    return;
}

