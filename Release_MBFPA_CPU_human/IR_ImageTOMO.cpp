#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h> // for Ubuntu

#include "physParamsTOMO.h"
#include "mallocMD.h"
#include "IR_ImageTOMO.h"

double double_TV(float a, float a1, float a2, float a3, float a4, float a7, float a8)
{
    double dist1 = sqrt(double((a-a2)*(a-a2)+(a-a4)*(a-a4)))+0.00000001;
    double dist2 = sqrt(double((a-a1)*(a-a1)+(a1-a7)*(a1-a7)))+0.00000001;
    double dist3 = sqrt(double((a3-a8)*(a3-a8)+(a-a3)*(a-a3)))+0.00000001;
    return ((a-a2)+(a-a4))/dist1+(a-a1)/dist2+(a-a3)/dist3;
    //    else return 0.0;
}

void
IR_ImageTOMO( const int ite, int reconSize, float* reconImagefloat_before, float* reconImagefloat,
              float* Prior_image_2D, float* npf_float, float* npf_float_re, double wbefore,
              double weight_tv, double weight_prior, char* MaskData, const int thinout)
{
  for(int kky = 0; kky < reconSize; kky++){
    for(int kkx = 0; kkx < reconSize; kkx++){
      long kk = kky*reconSize + kkx;
      if( MaskData[kk] != 0 ){
	double cost_tv = 0.0, cost_prior = 0.0;
	//TotalVariation
	if(kky < reconSize-1 && kkx < reconSize-1 && kkx > 0 && kky > 0){
	  int kk = kky*reconSize + kkx;
	  int k1 = kky*reconSize + kkx-1;
	  int k2 = kky*reconSize + kkx+1;
	  int k3 = (kky-1)*reconSize + kkx;
	  int k4 = (kky+1)*reconSize + kkx;
	  int k7 = (kky+1)*reconSize + kkx-1;
	  int k8 = (kky-1)*reconSize + kkx+1;
	  cost_tv = - double_TV(reconImagefloat_before[kk],
				reconImagefloat_before[k1],
				reconImagefloat_before[k2],
				reconImagefloat_before[k3],
				reconImagefloat_before[k4],
				reconImagefloat_before[k7],
				reconImagefloat_before[k8]);
	  cost_prior = - double_TV(reconImagefloat_before[kk]-Prior_image_2D[kk],
				   reconImagefloat_before[k1]-Prior_image_2D[k1],
				   reconImagefloat_before[k2]-Prior_image_2D[k2],
				   reconImagefloat_before[k3]-Prior_image_2D[k3],
				   reconImagefloat_before[k4]-Prior_image_2D[k4],
				   reconImagefloat_before[k7]-Prior_image_2D[k7],
				   reconImagefloat_before[k8]-Prior_image_2D[k8]);
	}
	
	double lambda = weight_tv * cost_tv + weight_prior * cost_prior + 0.000001;//for penalty term
	//Direct method
	//reconImagefloat[kk] = reconImagefloat_before[kk]*(
	//				   wbefore+(1-wbefore)*npf_float_re[kk]/(npf_float[kk]+lambda));
	//Gradient Decent
//	double del_F = npf_float_re[kk] - npf_float[kk] + lambda;//20171023
	    double del_F = npf_float_re[kk] - npf_float[kk]/thinout + lambda;//20171023 add
	    if(ite > 300) 
	    wbefore = reconImagefloat_before[kk]/(npf_float[kk]/thinout - lambda);
	//double stepsize = reconImagefloat_before[kk]/(npf_float[kk]-lambda);
	//double stepsize = 0.01;
	reconImagefloat[kk] = (float)(reconImagefloat_before[kk] + wbefore * del_F);
	if(reconImagefloat[kk] < 0.0) reconImagefloat[kk] = 0.0;
	//reconImagefloat[kk] = reconImagefloat_before[kk]*(1.0+npf_float_re[kk]/(npf_float[kk]+lambda))*0.5;
      }
    }
  }

}
