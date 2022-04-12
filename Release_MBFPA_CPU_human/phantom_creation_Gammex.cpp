#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main()
{
  ///////////////  Prior material image production ///////////////////
  FILE*	fp;
  FILE* fp3;
  char filename[128];
  int mnum=0, jmax;
  float dens[76], hdens[76], cdens[76], ndens[76], odens[76], pdens[76], cadens[76];

  // H, C, N, O, P, Ca, Mg, Si
  int num_material = 8;
  int reconsize = 512;
  float yscale = 1.0;
  float xscale = 1.0;
  int recon2Dsize = reconsize*reconsize*num_material;
  int reconsize2 = reconsize*reconsize;
  short*   Pri_reconImageshort =  (short*)malloc(recon2Dsize*sizeof(short) );
  short*   Pri_reconImageshort_dens =  (short*)malloc(reconsize2*sizeof(short) );
  /*
  if((fp = fopen("material_density_mod.txt","r")) == NULL ) {
                printf("Can not open\n");
                return;
        }
  while(fscanf(fp, "%d %f %f %f %f %f %f %f\n", &jmax, &dens[mnum], &hdens[mnum], &cdens[mnum], 
	       &ndens[mnum], &odens[mnum], &pdens[mnum], &cadens[mnum]) != EOF)
    {
      printf("%d %f %f \n", mnum, dens[mnum], hdens[mnum]);
      mnum++;
    }
  fclose(fp);
  
  int mat0[3] = {1,2,3};
  int mat11[6] = {4,9,14,26,39,49};
  int mat12[6] = {5,10,18,30,42,50};
  int mat13[6] = {6,11,22,35,45,51};
  int mat21[3] = {52,59,73};
  int mat22[3] = {53,63,75};
  int mat23[3] = {54,68,76};

  for( int i = 0; i < 3; i++)
    for(int j = 0; j < 6; j++)
      for(int k = 0; k < 3; k++)
	{
	  sprintf(filename,"ph3D300SL16bit.dat");///////////////////////////////////////////////////////////////////////////
	  int ict1 = 1310;
	  int ict2 = 1965;
	  if( (fp=fopen(filename,"rb")) == NULL )
	    {
	      fprintf(stderr,"ERROR: No output CT file \n");
	      return;
	    }
	  int a = fread( Pri_reconImageshort_0, 1, pri_size * sizeof(short), fp );
	  fclose(fp);
	  
	  int material1 = mat0[i]-1;
	  int material2 = mat11[j]-1;
	  if(2 == mat0[i]) material2 = mat12[j]-1;
	  if(3 == mat0[i]) material2 = mat13[j]-1;
	  int material3 = mat21[k]-1;
	  if(2 == mat0[i]) material3 = mat22[k]-1;
	  if(3 == mat0[i]) material3 = mat23[k]-1;
	  sprintf(filename,"ph3D300SL16bit_%d_%d_%d.dat", material1, material2, material3);///////////////////////////////////////////////////////////////////////////
	  fprintf(stderr,"%s %f %f %f \n", filename, dens[material1],dens[material2],dens[material3]);
	  fp=fopen(filename,"wb");
	  for(int ii = 0; ii < pri_size; ii++) 
	    {
	      int i0 = ii/reconSize/reconSize;
	      int i1 = (ii-i0*reconSize*reconSize)/reconSize;
	      int i2 = ii-i1*reconSize-i0*reconSize*reconSize;
	      if(Pri_reconImageshort_0[ii] < 0) Pri_reconImageshort_0[ii] = 0;
	      if(Pri_reconImageshort_0[ii] == ict1) Pri_reconImageshort_0[ii] = (short)(dens[material2]*1000);
	      else if(Pri_reconImageshort_0[ii] == ict2) Pri_reconImageshort_0[ii] = (short)(dens[material3]*1000);
	      else if(Pri_reconImageshort_0[ii] == 0 && (i1 <= 210 && i1 >= 80) && (i2 <= 195 && i2 >= 75) && (i0 <= 194 && i0 >= 108))  
		Pri_reconImageshort_0[ii] = (short)(dens[material1]*1000);
	      else Pri_reconImageshort_0[ii] = 0;
	      
	    }
	  fwrite( Pri_reconImageshort_0,
		  getByteSizeOfUnit(UNIT_SINT16),pri_size,fp);
	  fclose(fp);
	  
	  float rho1 = dens[material1]*1000;
	  float rho2 = dens[material2]*1000;
	  float rho3 = dens[material3]*1000;
	  float wei1[6], wei2[6], wei3[6];
	  wei1[0] = hdens[material1]; wei1[1] = cdens[material1]; wei1[2] = ndens[material1]; wei1[3] = odens[material1]; wei1[4] = pdens[material1], wei1[5] = cadens[material1];  
	  wei2[0] = hdens[material2]; wei2[1] = cdens[material2]; wei2[2] = ndens[material2]; wei2[3] = odens[material2]; wei2[4] = pdens[material2], wei2[5] = cadens[material2];  
	  wei3[0] = hdens[material3]; wei3[1] = cdens[material3]; wei3[2] = ndens[material3]; wei3[3] = odens[material3]; wei3[4] = pdens[material3], wei3[5] = cadens[material3];  
	  
	  //exit(0);
	  //
	  */
  float mat[14][8], dens_gammex[14];
  for(int i = 0; i<14; i++)
    for(int j = 0; j<8; j++)
      {
	mat[i][j] = 0.0;
	dens_gammex[i] = 0.0;
      }
  // Air
  dens_gammex[0] = 0.002;
  mat[0][0] = 0.0; 
  mat[0][1] = 0.0;
  mat[0][2] = 0.75;
  mat[0][3] = 0.232;
  mat[0][4] = 0.0;
  mat[0][5] = 0.018;
  // LN300
  dens_gammex[1] = 0.28;
  mat[1][0] = 8.46; 
  mat[1][1] = 59.38;
  mat[1][2] = 1.96;
  mat[1][3] = 18.14;
  mat[1][4] = 0.0;
  mat[1][5] = 0.0;
  mat[1][6] = 11.19;
  mat[1][7] = 0.78;
  // LN450
  dens_gammex[2] = 0.40;
  mat[2][0] = 8.47; 
  mat[2][1] = 59.57;
  mat[2][2] = 1.97;
  mat[2][3] = 18.11;
  mat[2][4] = 0.0;
  mat[2][5] = 0.0;
  mat[2][6] = 11.21;
  mat[2][7] = 0.58;
  // AP6
  dens_gammex[3] = 0.942;
  mat[3][0] = 9.06; 
  mat[3][1] = 72.30;
  mat[3][2] = 2.25;
  mat[3][3] = 16.27;
  mat[3][4] = 0.0;
  mat[3][5] = 0.0;
  // BR12
  dens_gammex[4] = 0.977;
  mat[4][0] = 8.59; 
  mat[4][1] = 70.11;
  mat[4][2] = 2.33;
  mat[4][3] = 17.90;
  mat[4][4] = 0.0;
  mat[4][5] = 0.95;
  // Solid water
  dens_gammex[5] = 1.018;
  mat[5][0] = 8.02; 
  mat[5][1] = 67.23;
  mat[5][2] = 2.41;
  mat[5][3] = 19.91;
  mat[5][4] = 0.0;
  mat[5][5] = 2.31;
  // SR2-Brain
  dens_gammex[6] = 1.053;
  mat[6][0] = 10.83; 
  mat[6][1] = 72.54;
  mat[6][2] = 1.69;
  mat[6][3] = 14.86;
  mat[6][4] = 0.0;
  mat[6][5] = 0.0;
  // LV1
  dens_gammex[7] = 1.097;
  mat[7][0] = 8.06; 
  mat[7][1] = 67.01;
  mat[7][2] = 2.47;
  mat[7][3] = 20.01;
  mat[7][4] = 0.0;
  mat[7][5] = 2.31;
  // IB3
  dens_gammex[8] = 1.143;
  mat[8][0] = 6.67; 
  mat[8][1] = 55.64;
  mat[8][2] = 1.96;
  mat[8][3] = 23.52;
  mat[8][4] = 3.23;
  mat[8][5] = 8.86;
  // B200
  dens_gammex[9] = 1.154;
  mat[9][0] = 6.65; 
  mat[9][1] = 55.52;
  mat[9][2] = 1.98;
  mat[9][3] = 23.64;
  mat[9][4] = 3.24;
  mat[9][5] = 8.87;
  // CB2-30%
  dens_gammex[10] = 1.335;
  mat[10][0] = 6.68; 
  mat[10][1] = 53.48;
  mat[10][2] = 2.12;
  mat[10][3] = 25.61;
  mat[10][4] = 0.0;
  mat[10][5] = 12.01;
  // CB2-50%
  dens_gammex[11] = 1.56;
  mat[11][0] = 4.77; 
  mat[11][1] = 41.63;
  mat[11][2] = 1.52;
  mat[11][3] = 32.00;
  mat[11][4] = 0.0;
  mat[11][5] = 20.02;
  // SB3-Cortical
  dens_gammex[12] = 1.825;
  mat[12][0] = 3.41; 
  mat[12][1] = 31.41;
  mat[12][2] = 1.84;
  mat[12][3] = 36.50;
  mat[12][4] = 0.0;
  mat[12][5] = 26.81;
  // Water
  dens_gammex[13] = 1.0;
  mat[13][0] = 1.00798*2/(1.00798*2+15.9994); 
  mat[13][1] = 0.0;
  mat[13][2] = 0.0;
  mat[13][3] = 15.9994/(1.00798*2+15.9994);
  mat[13][4] = 0.0;
  mat[13][5] = 0.0;

  for(int i = 1; i<13; i++)
    {
      for(int j = 0; j<8; j++)
	{
	  mat[i][j] = mat[i][j]/100;
	}
    }
  for(int i = 0; i<14; i++)
    {
      dens_gammex[i] = dens_gammex[i]*1000;
    }
  for(int i = 0; i<14; i++)
    {
      float sum_mat = 0.0;
      for(int j = 0; j<8; j++)
	{
	  sum_mat += mat[i][j];
	}
      printf("%f %f\n",sum_mat, sum_mat*dens_gammex[i]);
    }
  //exit(0);
  // Elements; H C N O P Ca Mg Si
  // Materials; Air, LN300, LN450,,,
  ///////// Material density (True) //////////
  float rs = 14.0;
  int in = 0, inn;
  float r0 = 55;
  float r00 = 105;
    for(int iphase = 0; iphase < num_material; iphase++)
      {
	inn = 0;
	for(int iy = 0; iy < reconsize; iy++)
	  {
	    float y = yscale*(reconsize/2 - iy);
	    for(int ix = 0; ix < reconsize; ix++) 
	      {
		float x = xscale*(ix - reconsize/2);
		float rr = sqrt(x*x + y*y);
		Pri_reconImageshort[in] = (int)(2*mat[0][iphase]);
		if(rr <= 330/2) Pri_reconImageshort[in] = (int)(dens_gammex[5]*mat[5][iphase]);  // Solid water
		if(rr <= 330/2) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[5]);  // Solid water
		//float r1 = sqrt(x*x + (y-r0)*(y-r0));
		//float r2 = sqrt((x-r0/sqrt(2.0))*(x-r0/sqrt(2.0)) + (y-r0/sqrt(2.0))*(y-r0/sqrt(2.0)));
		//float r3 = sqrt((x-r0)*(x-r0) + (y)*(y));
		//float r4 = sqrt((x-r0/sqrt(2.0))*(x-r0/sqrt(2.0)) + (y+r0/sqrt(2.0))*(y+r0/sqrt(2.0)));
		//float r5 = sqrt(x*x + (y+r0)*(y+r0));
		//float r6 = sqrt((x+r0/sqrt(2.0))*(x+r0/sqrt(2.0)) + (y+r0/sqrt(2.0))*(y+r0/sqrt(2.0)));
		//float r7 = sqrt((x+r0)*(x+r0) + (y)*(y));
		//float r8 = sqrt((x+r0/sqrt(2.0))*(x+r0/sqrt(2.0)) + (y-r0/sqrt(2.0))*(y-r0/sqrt(2.0)));

		float dtheta = 45*3.14159265/180;
		for(int ith = 0; ith < 8; ith++)
		  {
		    float theta = ith*dtheta + 22.5*3.14159265/180;
		    float x1 = r00*sin(theta);
		    float y1 = r00*cos(theta);
		    float r01 = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1));
		    if(r01 <= rs && ith == 0) Pri_reconImageshort[in] = (int)(dens_gammex[1]*mat[1][iphase]);
		    else if(r01 <= rs && ith == 1) Pri_reconImageshort[in] = (int)(dens_gammex[5]*mat[5][iphase]);
		    else if(r01 <= rs && ith == 2) Pri_reconImageshort[in] = (int)(dens_gammex[8]*mat[8][iphase]);
		    else if(r01 <= rs && ith == 3) Pri_reconImageshort[in] = (int)(dens_gammex[5]*mat[5][iphase]);
		    else if(r01 <= rs && ith == 4) Pri_reconImageshort[in] = (int)(dens_gammex[7]*mat[7][iphase]);
		    else if(r01 <= rs && ith == 5) Pri_reconImageshort[in] = (int)(dens_gammex[9]*mat[9][iphase]);
		    else if(r01 <= rs && ith == 6) Pri_reconImageshort[in] = (int)(dens_gammex[2]*mat[2][iphase]);
		    else if(r01 <= rs && ith == 7) Pri_reconImageshort[in] = (int)(dens_gammex[5]*mat[5][iphase]);
		    if(r01 <= rs && ith == 0) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[1]);
		    else if(r01 <= rs && ith == 1) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[5]);
		    else if(r01 <= rs && ith == 2) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[8]);
		    else if(r01 <= rs && ith == 3) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[5]);
		    else if(r01 <= rs && ith == 4) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[7]);
		    else if(r01 <= rs && ith == 5) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[9]);
		    else if(r01 <= rs && ith == 6) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[2]);
		    else if(r01 <= rs && ith == 7) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[5]);
		  }
		float r1 = sqrt(x*x + (y-r0)*(y-r0));
		float r2 = sqrt((x-r0/sqrt(2.0))*(x-r0/sqrt(2.0)) + (y-r0/sqrt(2.0))*(y-r0/sqrt(2.0)));
		float r3 = sqrt((x-r0)*(x-r0) + (y)*(y));
		float r4 = sqrt((x-r0/sqrt(2.0))*(x-r0/sqrt(2.0)) + (y+r0/sqrt(2.0))*(y+r0/sqrt(2.0)));
		float r5 = sqrt(x*x + (y+r0)*(y+r0));
		float r6 = sqrt((x+r0/sqrt(2.0))*(x+r0/sqrt(2.0)) + (y+r0/sqrt(2.0))*(y+r0/sqrt(2.0)));
		float r7 = sqrt((x+r0)*(x+r0) + (y)*(y));
		float r8 = sqrt((x+r0/sqrt(2.0))*(x+r0/sqrt(2.0)) + (y-r0/sqrt(2.0))*(y-r0/sqrt(2.0)));
		if(r1 <= rs) Pri_reconImageshort[in] = (int)(dens_gammex[12]*mat[12][iphase]);
		else if(r2 <= rs ) Pri_reconImageshort[in] = (int)(dens_gammex[6]*mat[6][iphase]);
		else if(r3 <= rs ) Pri_reconImageshort[in] = (int)(dens_gammex[3]*mat[3][iphase]);
		else if(r4 <= rs ) Pri_reconImageshort[in] = (int)(dens_gammex[11]*mat[11][iphase]);
		else if(r5 <= rs ) Pri_reconImageshort[in] = (int)(dens_gammex[5]*mat[5][iphase]);
		else if(r6 <= rs ) Pri_reconImageshort[in] = (int)(dens_gammex[10]*mat[10][iphase]);
		else if(r7 <= rs ) Pri_reconImageshort[in] = (int)(dens_gammex[13]*mat[13][iphase]);
		else if(r8 <= rs ) Pri_reconImageshort[in] = (int)(dens_gammex[4]*mat[4][iphase]);
		if(r1 <= rs) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[12]);
		else if(r2 <= rs ) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[6]);
		else if(r3 <= rs ) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[3]);
		else if(r4 <= rs ) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[11]);
		else if(r5 <= rs ) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[5]);
		else if(r6 <= rs ) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[10]);
		else if(r7 <= rs ) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[13]);
		else if(r8 <= rs ) Pri_reconImageshort_dens[inn] = (int)(dens_gammex[4]);
		in++;
		inn++;
		      //Pri_reconImageshort[i+pri_size*iphase] = 0;
		      //if(Pri_reconImageshort_0[i] == (short)(dens[material1]*1000)) Pri_reconImageshort[i+pri_size*iphase] = (short)(wei1[iphase]*1000);
		      //if(Pri_reconImageshort_0[i] == (short)(dens[material2]*1000)) Pri_reconImageshort[i+pri_size*iphase] = (short)(wei2[iphase]*1000);
		      //if(Pri_reconImageshort_0[i] == (short)(dens[material3]*1000)) Pri_reconImageshort[i+pri_size*iphase] = (short)(wei3[iphase]*1000);
	      }
	  }
      }
  
  sprintf(filename,"Weight_input_gammex.raw");
  fp3=fopen(filename,"wb");	      
  fwrite( Pri_reconImageshort,
	  sizeof(short),reconsize*reconsize*num_material,fp3);
  fclose(fp3);
  free(Pri_reconImageshort);
  
  sprintf(filename,"dens_input_gammex.raw");
  fp3=fopen(filename,"wb");	      
  fwrite( Pri_reconImageshort_dens,
	  sizeof(short),reconsize*reconsize,fp3);
  fclose(fp3);
  free(Pri_reconImageshort_dens);  
	    /*
	  fwrite( req_reconImageshort,
		  getByteSizeOfUnit(UNIT_SINT16),
		  reconSize*reconSize*num_material,fp3);
	  fclose(fp3);
	    */
	  
	  //for(int i = 0; i < projImgWidth*projImgHeight; i++)
	  //  projData->projVolume[i] = ProjImagefloat[i];
 	  //////////////////////// Virtual projection image production -- END -- ////////////////////////////	  
	    //	}
  return 0;
}

