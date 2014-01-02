#include "contour.h"

// uchar_data : volume data which has unsigned char values as an array
// dim        : integer array with 3 elements. dim[0]:x,dim[1]:y,dim[2]:z
// array_size : will be set as the array size of isoval,area,...,gradient
//              usually 256
// isoval     : array with regularly increasing from min to max values
// area,min_vol,max_vol,gradient : array of quantitative properties
//            : needs to be normalized
void getContourSpectrum(unsigned char* uchar_data, int* dim, int&
 array_size, float* isoval , float* area, float* min_vol, float* max_vol,
 float* gradient)
{
 	int i;
	array_size=256;
 	ConDataset* the_data;

 	the_data = newDatasetReg(CONTOUR_UCHAR, CONTOUR_REG_3D, 1, 1, dim,uchar_data);
	
	Signature	*sig;
	sig=getSignatureFunctions(the_data, 0,0);

	for (i=0;i<array_size;i++) {
		isoval[i]=sig[0].fx[i];
		area[i]=sig[0].fy[i];
		min_vol[i]=sig[1].fy[i];
		max_vol[i]=sig[2].fy[i];
		gradient[i]=sig[3].fy[i];
	}

	delete the_data;
	delete sig;
}
