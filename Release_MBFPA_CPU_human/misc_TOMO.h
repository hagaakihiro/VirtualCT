
char*	prepareMaskData(int reconSize);

void	clearReconstructedImage(int reconSize, int zslices, double* reconVoxData, short* reconVoxDataShort);
void    clear_npf(int recon2Dsize, float* npf);
void    clearSinogramImage(int projsize, float* sinogram);

//void	shiftAngles(int numAngles,double* angles,double dAngle);
//void	prepareAngles(int numAngles,double* angles,int reconSize,int** intStartAngle,int** intEndAngle);
int     prepareTrigonometricTable(int numAngles,double* anglesRad,double** cosTable,double** sinTable);

void	startTimer();
void	stopTimer();
void	showLapTime(char* label);



