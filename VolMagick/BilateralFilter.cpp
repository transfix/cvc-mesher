/* $Id: BilateralFilter.cpp,v 1.2 2008/02/01 20:12:12 transfix Exp $ */

#include <math.h>
#include "VolMagick.h"

namespace VolMagick
{
  Voxels& Voxels::bilateralFilter(double radiometricSigma, double spatialSigma, unsigned int filterRadius)
  {
    int i,j,k,x,y,z,c,index;
    int filterDiameter = filterRadius*2+1;
    double fsample, weight, normalizedDiff, factor;
    double sum, denom;
    bool bool1, bool2;
    double radiometricTable[256];
    double *spatialMask = new double[filterDiameter*filterDiameter*filterDiameter];

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::BilateralFilter, ZDim());

    //compute the radiometric table
    for(c=0; c<256; c++)
      {
	factor = c * ((max() - min()) / 255.0);
	radiometricTable[c] = exp((double)(factor*factor)/(-radiometricSigma*radiometricSigma*2.0));
      }

    //compute the spatial mask
    index = 0;
    for (k=-int(filterRadius); k<=int(filterRadius); k++)
      for (j=-int(filterRadius); j<=int(filterRadius); j++)
	for (i=-int(filterRadius); i<=int(filterRadius); i++)
	  spatialMask[index++] = exp((double)(k*k+j*j+i*i)/(-spatialSigma*spatialSigma*2.0));
    
    for(k=0; k<int(ZDim()); k++)
      {
	for(j=0; j<int(YDim()); j++)
	  for(i=0; i<int(XDim()); i++)
	    {
	      fsample = (*this)(i,j,k);
	      sum = 0; denom = 0;
	      
	      for(z=0; z<filterDiameter; z++)
		{
		  bool1 = k+z>=int(filterRadius) && k+z<int(ZDim())+int(filterRadius);
		  for(y=0; y<filterDiameter; y++)
		    {
		      bool2 = bool1 && (j+y>=int(filterRadius) && j+y<int(YDim())+int(filterRadius));
		      for(x=0; x<filterDiameter; x++)
			{
			  if(i+x>=int(filterRadius) && i+x<int(XDim())+int(filterRadius) && bool2)
			    {
			      normalizedDiff = fabs(fsample - (*this)(i+x-filterRadius,j+y-filterRadius,k+z-filterRadius));
			      normalizedDiff /= (max()-min());
			      normalizedDiff *= 255.0;
			      weight = radiometricTable[int(normalizedDiff)]*
				spatialMask[z*filterDiameter*filterDiameter+y*filterDiameter+x];
			      denom += weight;
			      sum += weight * (*this)(i+x-filterRadius,j+y-filterRadius,k+z-filterRadius);
			    }
			}
		    }
		}
	      
	      (*this)(i,j,k, sum/denom);
	    }

	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::BilateralFilter, k);
      }

    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::BilateralFilter);

    delete [] spatialMask;
    return *this;
  }
};
