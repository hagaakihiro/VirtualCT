
typedef struct projectionData
{
	int     usedProjNumber;

    int     totalProjNumber;
	int     projImgWidth;
	int     projImgHeight;

    char	projImgDataType[16];
    
	char	projImgFileName[256];
    
    char    AirScanImgFileName[256];

	int     projImgFileOffset;

	int     reconSize;
    double   reconPixSize;
	char	reconTgtClassString[16];
	double   reconOffset;
	double   reconStep;
	int     reconSlices;
	int     xfactor;
	char	InPrClassString[16];
    int     datasize;
    int     datatype;
    
    
    double GantryRotationPeriod;     // Gantry rotation period [s]
    double DataSamplingRate;         // Data Sampling Rate [Hz]
    double StartViewAngle;           // Starting View Angle [degree]
    double CouchSpeed;               // Couch Speed [mm/rotation]
    double StartCouchPosition;       // Start couch position [mm]
    
    int     ReconstructionFilter;
    
	char	reconDataFilename[256];

	double*	anglesRad;
	double*	xshift;
	double*	yshift;
    int*    inteangle;

	float**	projVolume;

} PROJECTION_DATA;


int	loadData(char* infoFilename,PROJECTION_DATA* projData);

PROJECTION_DATA*	newProjectionData();
void	deleteProjectionData(PROJECTION_DATA* projData);

int	loadProjectionData(PROJECTION_DATA* projData);
int	loadPixelData( char* filename, int offset, int length, unsigned short* pixels );

