void
lengthproj_TOMO( int ite, int reconSize, float* npf_float, float* projVolume,
		 double reconScale,
		 int projImgWidth, int projImgHeight,
		 double* cosTable, double* sinTable,
		 double zstart,
		 double CouchSpeed,
		 char* MaskData,
		 float* systemMatrixfloat, float* xraytracedetector,
		 int thinout);

void
lengthproj_TOMO_2( int ite, int reconSize, float* npf_float, float* projVolume,
		   double reconScale,
		   int projImgWidth, int projImgHeight,
		   double zstart,
		   double CouchSpeed,
		   char* MaskData,
		   float* systemMatrixfloat, float* xraytracedetector,
		   int thinout);

