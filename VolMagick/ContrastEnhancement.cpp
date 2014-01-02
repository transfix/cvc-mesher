/* $Id: ContrastEnhancement.cpp,v 1.2 2008/02/01 20:12:12 transfix Exp $ */

#include <math.h>
#include "VolMagick.h"

namespace VolMagick
{
  static inline void contrastEnhancementSlice(double resistor,
					      Voxels *paramin,
					      Voxels *paramax,
					      Voxels *imgavg,
					      unsigned int k)
  {
    int i, j;
    int xdim = paramin->XDim(), ydim = paramin->YDim();

    Voxels tmpmin(Dimension(xdim,ydim,1),paramin->voxelType());
    Voxels tmpmax(Dimension(xdim,ydim,1),paramin->voxelType());
    Voxels lcmin(Dimension(xdim,ydim,1),paramin->voxelType());
    Voxels lcmax(Dimension(xdim,ydim,1),paramin->voxelType());

    for(j=0; j<ydim; j++)
      for(i=0; i<xdim; i++)
	{
	  lcmin(i,j,0,  (*paramin)(i,j,k));
	  lcmax(i,j,0,  (*paramax)(i,j,k));
	  tmpmin(i,j,0, (*paramin)(i,j,k)); 
	  tmpmax(i,j,0, (*paramax)(i,j,k));
	}

    /* Bottom-up */
    for(i=1; i<xdim; i++)
      {
	(*imgavg)(i,0,k, (*imgavg)(i,0,k) + resistor*((*imgavg)(i-1,0,k)-(*imgavg)(i,0,k)));
	if(tmpmin(i-1) < tmpmin(i))
	  tmpmin(i, tmpmin(i) + resistor*(tmpmin(i-1)-tmpmin(i)));
	if(tmpmax(i-1) > tmpmax(i))
	  tmpmax(i, tmpmax(i) + resistor*(tmpmax(i-1)-tmpmax(i)));
      }
   
    for(i=xdim-2; i>=0; i--)
      {
	(*imgavg)(i,0,k, (*imgavg)(i,0,k) + resistor*((*imgavg)(i+1,0,k)-(*imgavg)(i,0,k)));
	if(tmpmin(i+1) < tmpmin(i))
	  tmpmin(i, tmpmin(i) + resistor*(tmpmin(i+1)-tmpmin(i)));
	if(tmpmax(i+1) > tmpmax(i))
	  tmpmax(i, tmpmax(i) + resistor*(tmpmax(i+1)-tmpmax(i)));
      }

    for(j=1; j<ydim; j++)
      {
	for(i=0; i<xdim; i++)
	  {
	    (*imgavg)(i,j,k, (*imgavg)(i,j,k) + resistor*((*imgavg)(i,j-1,k)-(*imgavg)(i,j,k)));
	    if(tmpmin(i,j-1,0) < tmpmin(i,j,0))
	      tmpmin(i,j,0, tmpmin(i,j,0) + resistor*(tmpmin(i,j-1,0)-tmpmin(i,j,0)));
	    if(tmpmax(i,j-1,0) > tmpmax(i,j,0))
	      tmpmax(i,j,0, tmpmax(i,j,0) + resistor*(tmpmax(i,j-1,0)-tmpmax(i,j,0)));
	  }
	
	for(i=1; i<xdim; i++)
	  {
	    (*imgavg)(i,j,k, (*imgavg)(i,j,k) + resistor*((*imgavg)(i-1,j,k)-(*imgavg)(i,j,k)));
	    if(tmpmin(i-1,j,0) < tmpmin(i,j,0))
	      tmpmin(i,j,0, tmpmin(i,j,0) + resistor*(tmpmin(i-1,j,0)-tmpmin(i,j,0)));
	    if(tmpmax(i-1,j,0) > tmpmax(i,j,0))
	      tmpmax(i,j,0, tmpmax(i,j,0) + resistor*(tmpmax(i-1,j,0)-tmpmax(i,j,0)));
	  }
	
	for(i=xdim-2; i>=0; i--)
	  {
	    (*imgavg)(i,j,k, (*imgavg)(i,j,k) + resistor*((*imgavg)(i+1,j,k)-(*imgavg)(i,j,k)));
	    if(tmpmin(i+1,j,0) < tmpmin(i,j,0))
	      tmpmin(i,j,0, tmpmin(i,j,0) + resistor*(tmpmin(i+1,j,0)-tmpmin(i,j,0)));
	    if(tmpmax(i+1,j,0) > tmpmax(i,j,0))
	      tmpmax(i,j,0, tmpmax(i,j,0) + resistor*(tmpmax(i+1,j,0)-tmpmax(i,j,0)));
	  }
      }

    /* Top-down */
    j=ydim-1;
    for(i=1; i<xdim; i++)
      {
	(*imgavg)(i,j,k, (*imgavg)(i,j,k) + resistor*((*imgavg)(i-1,j,k)-(*imgavg)(i,j,k)));
	if(lcmin(i-1,j,0) < lcmin(i,j,0))
	  lcmin(i,j,0, lcmin(i,j,0) + resistor*(lcmin(i-1,j,0)-lcmin(i,j,0)));
	if(lcmax(i-1,j,0) > lcmax(i,j,0))
	  lcmax(i,j,0, lcmax(i,j,0) + resistor*(lcmax(i-1,j,0)-lcmax(i,j,0)));
      }
   
    for(i=xdim-2; i>=0; i--)
      {
	(*imgavg)(i,j,k, (*imgavg)(i,j,k) + resistor*((*imgavg)(i+1,j,k)-(*imgavg)(i,j,k)));
	if(lcmin(i+1,j,0) < lcmin(i,j,0))
	  lcmin(i,j,0, lcmin(i,j,0) + resistor*(lcmin(i+1,j,0)-lcmin(i,j,0)));
	if(lcmax(i+1,j,0) > lcmax(i,j,0))
	  lcmax(i,j,0, lcmax(i,j,0) + resistor*(lcmax(i+1,j,0)-lcmax(i,j,0)));
      }
  
    for(j=ydim-2; j>=0; j--)
      {
	for(i=0; i<xdim; i++)
	  {
	    (*imgavg)(i,j,k, (*imgavg)(i,j,k) + resistor*((*imgavg)(i,j+1,k)-(*imgavg)(i,j,k)));
	    if(lcmin(i,j+1,0) < lcmin(i,j,0))
	      lcmin(i,j,0, lcmin(i,j,0) + resistor*(lcmin(i,j+1,0)-lcmin(i,j,0)));
	    if(lcmax(i,j+1,0) > lcmax(i,j,0))
	      lcmax(i,j,0, lcmax(i,j,0) + resistor*(lcmax(i,j+1,0)-lcmax(i,j,0)));
	  }

	for(i=1; i<xdim; i++)
	  {
	    (*imgavg)(i,j,k, (*imgavg)(i,j,k) + resistor*((*imgavg)(i-1,j,k)-(*imgavg)(i,j,k)));
	    if(lcmin(i-1,j,0) < lcmin(i,j,0))
	      lcmin(i,j,0, lcmin(i,j,0) + resistor*(lcmin(i-1,j,0)-lcmin(i,j,0)));
	    if(lcmax(i-1,j,0) > lcmax(i,j,0))
	      lcmax(i,j,0, lcmax(i,j,0) + resistor*(lcmax(i-1,j,0)-lcmax(i,j,0)));
	  }
    
	for(i=xdim-2; i>=0; i--)
	  {
	    (*imgavg)(i,j,k, (*imgavg)(i,j,k) + resistor*((*imgavg)(i+1,j,k)-(*imgavg)(i,j,k)));
	    if(lcmin(i+1,j,0) < lcmin(i,j,0))
	      lcmin(i,j,0, lcmin(i,j,0) + resistor*(lcmin(i+1,j,0)-lcmin(i,j,0)));
	    if(lcmax(i+1,j,0) > lcmax(i,j,0))
	      lcmax(i,j,0, lcmax(i,j,0) + resistor*(lcmax(i+1,j,0)-lcmax(i,j,0)));
	  }
      }
    
    for(j=0; j<ydim; j++)
      for(i=0; i<xdim; i++)
	{
	  (*paramin)(i,j,k, MIN(lcmin(i,j,0),tmpmin(i,j,0)));
	  (*paramax)(i,j,k, MAX(lcmax(i,j,0),tmpmax(i,j,0)));
	}
  }

  Voxels& Voxels::contrastEnhancement(double resistor)
  {
    int i,j,k, curstep=0;
    double origmin, origmax, lmin, lmax, img, avg;
    double window, a, b, c, alpha;

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::ContrastEnhancement, ZDim()*3);

    //clamp the resistor value between [0.0,1.0]
    resistor = MIN(1.0,MAX(0.0,resistor));
    
    origmin = min();
    origmax = max();

    map(0.0,255.0); //not sure if this is necessary

    Voxels upmin(*this), upmax(*this);
    Voxels downmin(*this), downmax(*this), imgavg(*this);

    /* Bottom-up propagation */
    contrastEnhancementSlice(resistor,&upmin,&upmax,&imgavg,0);
    if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::ContrastEnhancement, curstep++);

    for(k=1; k<int(ZDim()); k++)
      {
	/* propagation from lower slice */
	for(j=0; j<int(YDim()); j++)
	  for(i=0; i<int(XDim()); i++)
	    {
	      imgavg(i,j,k, imgavg(i,j,k) + resistor*(imgavg(i,j,k-1)-imgavg(i,j,k)));
	      if(upmin(i,j,k-1) < upmin(i,j,k))
		upmin(i,j,k, upmin(i,j,k) + resistor*(upmin(i,j,k-1)-upmin(i,j,k)));
	      if(upmax(i,j,k-1) > upmax(i,j,k))
		upmax(i,j,k, upmax(i,j,k) + resistor*(upmax(i,j,k-1)-upmax(i,j,k)));
	    }
	contrastEnhancementSlice(resistor,&upmin,&upmax,&imgavg,k);
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::ContrastEnhancement, curstep++);
      }

    /* Top-down propagation */
    contrastEnhancementSlice(resistor,&downmin,&downmax,&imgavg,ZDim()-1);
    if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::ContrastEnhancement, curstep++);

    for(k=ZDim()-2; k>=0; k--)
      {
	/* propagation from upper slice */
	for(j=0; j<int(YDim()); j++)
	  for(i=0; i<int(XDim()); i++) 
	    {
	      imgavg(i,j,k, imgavg(i,j,k) + resistor*(imgavg(i,j,k+1)-imgavg(i,j,k)));
	      if(downmin(i,j,k+1) < downmin(i,j,k))
		downmin(i,j,k, downmin(i,j,k) + resistor*(downmin(i,j,k+1)-downmin(i,j,k)));
	      if(downmax(i,j,k+1) > downmax(i,j,k))
		downmax(i,j,k, downmax(i,j,k) + resistor*(downmax(i,j,k+1)-downmax(i,j,k)));
	    }
	contrastEnhancementSlice(resistor,&downmin,&downmax,&imgavg,k);
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::ContrastEnhancement, curstep++);
      }

    /* stretching */
    for(k=0; k<int(ZDim()); k++)
      {
	for(j=0; j<int(YDim()); j++)
	  for(i=0; i<int(XDim()); i++)
	    {
	      lmin = MIN(upmin(i,j,k),downmin(i,j,k));
	      lmax = MAX(upmax(i,j,k),downmax(i,j,k));
	      img = (*this)(i,j,k);
	      avg = imgavg(i,j,k);
	      window = lmax - lmin;
	      window = sqrt(window*(510-window));

	      if(lmin != lmax)
		{
		  img = window*(img-lmin)/(lmax-lmin);
		  avg = window*(avg-lmin)/(lmax-lmin);
		}

	      alpha = (avg-img)/(181.019 * window);
	      if(alpha != 0) 
		{
		  a = 0.707*alpha;
		  b = 1.414*alpha*(img - window) - 1;
		  c = 0.707*alpha*img*(img-2*window) + img;
		  imgavg(i,j,k, lmin+(-b-sqrt(b*b-4*a*c))/(2*a));
		}
	      else 
		imgavg(i,j,k, img+lmin);
	    }
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::ContrastEnhancement, curstep++);
      }

    imgavg.unsetMinMax(); //we need to recalculate min and max for imgavg
    copy(imgavg);

    map(origmin,origmax); //restore the original min/max and relative values

    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::ContrastEnhancement);
    
    return *this;
  }
};
