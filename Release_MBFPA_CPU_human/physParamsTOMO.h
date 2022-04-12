#define PI	3.141592653589793238462643383279

//Elekta XVI system
//#define DETECTOR_PITCH_CM             0.04
//#define DETECTOR_PITCH_MM_AT_ISO_CENT	0.26
//#define DIST_BTWN_SRC_AND_ISOCENT_MM	1000.0
//#define MAX_PIXVALUE                  65536.0
//#define AIR_PIXVALUE                  10000.0
//#define MAX_AIR_PIXDIFF	(MAX_PIXVALUE-AIR_PIXVALUE)

//#define MAX_PIXVALUE                    65536.0
//#define AIR_PIXVALUE                    30000.0
//#define MAX_AIR_PIXDIFF	(MAX_PIXVALUE-AIR_PIXVALUE)

//WEB84 reconstruction parameters
//#define CUTOFF_FREQ	1.60
//#define GAMMA		1.00
//#define ALPHA		0.50
//CUDA

//#define DimBlk_x 10
//#define DimBlk_y 10
// Higuchi Geometry
#define DETECTOR_PITCH_CM               		0.15
#define DIST_BTWN_SRC_AND_ISOCENT_CM			85.0 //57.0   // 85 in Tomo
#define DIST_BTWN_SRC_AND_DETECTOR_CM			142.6 //104.0  // 142.6 in Tomo
#define DETECTOR_RADIUS_CURVATURE_CM    		142.6 //104.0  // 99.8 in Tomo
#define DIST_BTWN_DETERCTORSRC_AND_ISOCENT_CM   	85.0 //57.0    // 42.2 in Tomo
#define DIST_BTWN_ISO_AND_DETECTOR_CM			57.6 //57.0   // 57.6 in Tomo = 142.6-85

#define Delta_Equal_angle_spacing       		DETECTOR_PITCH_CM/DIST_BTWN_SRC_AND_DETECTOR_CM
#define NumberOfChannels                		609
#define NumberOfSkips                		        27 // 2017.11.11
//Minnesota
//#define Start_Channel                   		65
//#define End_Channel                     		592
//#define isocenter_channel               		328.01

//Tokyo
//#define Start_Channel                   		27
//#define End_Channel                     		554
//#define isocenter_channel               		290.15
#define Start_Channel                   		0      // 2017.11.11
//#define End_Channel                     		527    // 2017.11.11
//#define isocenter_channel               		263.15  // 2017.11.11
#define End_Channel                     		609  // 2017.11.11
#define isocenter_channel               		304.5  // 2017.11.11

//Madison
//#define Start_Channel                   		27
//#define End_Channel                     		554
//#define isocenter_channel               		290.15

//#define Max_alpha					0.33003507 // Detector pitch = 0.125
//#define Max_alpha					0.3301143//0.231140351//0.3301143  // Detector pitch = 0.12503
#define Max_alpha					0.326086956 // Detector pitch = 0.15 and Channel = 620
#define Delta_Equal_angle_spacing_alpha			Max_alpha/((End_Channel-Start_Channel-1)/2)
//#define Max_beta					0.231177
//#define Max_beta					0.23219637
#define Max_beta					0.326086956 // Detector pitch = 0.15 and Channel = 620
//#define Max_beta                                        0.231140351 
#define Corrected_Delta_Equal_angle_spacing 	        Max_beta/((End_Channel-Start_Channel-1)/2)

#define AirScn_Sinogram_Height                          800
#define AirScn_Sinogram_Number                          5

#define UsedReconAngle					2 // if = 2 => 360 degree reconstruction, if = 4 => 180 degree reconstruction

#define SHOW_OK	fprintf(stderr,"OK\n");

// Treatment Planning CT
#define pCTsize					512
#define pCTw					0.2   // cm
#define pCTpixelsize                            0.107421875   // cm

// pre MV CT
#define pMVCTsize				512
#define pMVCTw					0.1   // cm
#define pMVCTpixelsize                          0.1   // cm



// For gaussian noise
#define MAX_PIXVALUE		131071.0
#define AIR_PIXVALUE		0.0
#define MAX_AIR_PIXDIFF	(MAX_PIXVALUE-AIR_PIXVALUE)

#define Rho_e_H2O   3.343
#define Avogadro    6.0221
