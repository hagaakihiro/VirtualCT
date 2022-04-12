
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
		   int thinout, float* Attenuation, float* pdf_pe, int num_mat, int NE, int ke);

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
		     int thinout, float* Attenuation, float* pdf_pe, int num_mat, int NE, int ke);




