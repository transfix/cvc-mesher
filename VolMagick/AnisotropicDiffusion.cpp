/* $Id: AnisotropicDiffusion.cpp,v 1.2 2008/02/01 20:12:12 transfix Exp $ */

#include <math.h>
#include "VolMagick.h"

namespace VolMagick
{
  Voxels& Voxels::anisotropicDiffusion(unsigned int iterations)
  {
    double cn, cs, ce, cw, cu, cd;
    double delta_n, delta_s, delta_e, delta_w, delta_u, delta_d;
    double K_para = 3;
    double Lamda_para = 0.16f;
    unsigned int m;
    VolMagick::uint64 i,j,k, stepnum = 0;
    Voxels tempt(*this);

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::AnisotropicDiffusion, ZDim()*iterations);

    for (m=0; m<iterations; m++)
      {
        for (k=0; k<ZDim(); k++)
	  {
	    for (j=0; j<YDim(); j++) 
	      for (i=0; i<XDim(); i++)
		{
		if (j < YDim()-1)
		  delta_s = (*this)(i,j+1,k) - (*this)(i,j,k);
		else
		  delta_s = 0.0;
		if (j > 0)
		  delta_n = (*this)(i,j-1,k) - (*this)(i,j,k);
		else 
		  delta_n = 0.0;
		if (i < XDim()-1)
		  delta_e = (*this)(i+1,j,k) - (*this)(i,j,k);
		else 
		  delta_e = 0.0;
		if (i > 0)
		  delta_w = (*this)(i-1,j,k) - (*this)(i,j,k);
		else
		  delta_w = 0.0;
		if (k < ZDim()-1)
		  delta_u = (*this)(i,j,k+1) - (*this)(i,j,k);
		else 
		  delta_u = 0.0;
		if (k > 0)
		  delta_d = (*this)(i,j,k-1) - (*this)(i,j,k);
		else
		  delta_d = 0.0;
		
		cn = 1.0f / (1.0f + ((delta_n * delta_n) / (K_para * K_para)));
		cs = 1.0f / (1.0f + ((delta_s * delta_s) / (K_para * K_para)));
		ce = 1.0f / (1.0f + ((delta_e * delta_e) / (K_para * K_para)));
		cw = 1.0f / (1.0f + ((delta_w * delta_w) / (K_para * K_para)));
		cu = 1.0f / (1.0f + ((delta_u * delta_u) / (K_para * K_para)));
		cd = 1.0f / (1.0f + ((delta_d * delta_d) / (K_para * K_para)));
	  
		tempt(i,j,k, 
		      (*this)(i,j,k) + 
		      Lamda_para * (cn * delta_n + cs * delta_s + ce * delta_e + 
				    cw * delta_w + cu * delta_u + cd * delta_d));
		}
	    
	    if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::AnisotropicDiffusion, stepnum++);
	  }
    
	(*this) = tempt;
      }

    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::AnisotropicDiffusion);

    return *this;
  }
}
