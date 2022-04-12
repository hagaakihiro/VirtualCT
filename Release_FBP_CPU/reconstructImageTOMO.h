void
reconstructImage( int reconSize, float*** reconVoxData, char** maskData,
                 int sngrmWidth, int sngrmHeight, float** sngrmVoxData,
                 double* cosTable, double* sinTable, 
				 double reconOffset, double reconStep, int reconSlices,
				 double StartCouchPosition, double CouchSpeed, int usedProjNumber,
				 double xfactor, double scale,
                 char* InPrClassString
                );
