
double**	prepareConstMap(int projWidth, int ReconstructionFilter);

void    prepareCorrectedDetectorInterpolation(int projWidth, double* Corrected_Detector_Weight0,
                                              double* Corrected_Detector_Weight1, double*Corrected_Detector_Weight2);


void	filterSinogram( int projHeight, int zsliceHeight_min, int zsliceHeight_max,
				float** projVolume, float** fltdVoxData,
				double** constMap,
                double* Corrected_Detector_Weight0, double* Corrected_Detector_Weight1, double* Corrected_Detector_Weight2,
                double* xshift, double* yshift,
                char* InPrClassString
                       );
