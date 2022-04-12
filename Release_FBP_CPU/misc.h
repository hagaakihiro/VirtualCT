
char**	prepareMaskData(int reconSize);

void	clearReconstructedImage(int reconSize, int zslices, float*** reconVoxData, short*** reconVoxDataShort);
void    clearSinogramImage(int projWidth, int projHeight_zslicesHeight, float** fltdSinogram);

void	shiftAngles(int numAngles,double* angles,double dAngle);
void	prepareAngles(int numAngles,double* angles,int reconSize,int** intStartAngle,int** intEndAngle);
int     prepareTrigonometricTable(int numAngles,double* anglesRad,double** cosTable,double** sinTable);

void	startTimer();
void	stopTimer();
void	showLapTime(char* label);



