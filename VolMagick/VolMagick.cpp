/* $Id: VolMagick.cpp,v 1.4 2008/08/15 21:53:04 transfix Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <new>

#include "VolMagick.h"
#include "endians.h"

#ifdef _MSC_VER
#define SNPRINTF _snprintf
#else
#define SNPRINTF snprintf
#endif

//trilinear interpolation function (see VolMagick::Voxels::resize())
static inline double getTriVal(double val[8], double x, double y, double z,
			       double resX, double resY, double resZ)
{
  double x_ratio, y_ratio, z_ratio;
  double temp1,temp2,temp3,temp4,temp5,temp6;
  
  x_ratio=x/resX;
  y_ratio=y/resY;
  z_ratio=z/resZ;
  
  if( x_ratio == 1 ) x_ratio = 0;
  if( y_ratio == 1 ) y_ratio = 0;
  if( z_ratio == 1 ) z_ratio = 0;
	
  temp1 = val[0] + (val[1]-val[0])*x_ratio;
  temp2 = val[4] + (val[5]-val[4])*x_ratio;
  temp3 = val[2] + (val[3]-val[2])*x_ratio;
  temp4 = val[6] + (val[7]-val[6])*x_ratio;
  temp5 = temp1  + (temp3-temp1)*y_ratio;
  temp6 = temp2  + (temp4-temp2)*y_ratio;
  
  return temp5  + (temp6-temp5)*z_ratio;
}

static inline void geterrstr(int errnum, char *strerrbuf, size_t buflen)
{
#ifdef HAVE_STRERROR_R
  strerror_r(errnum,strerrbuf,buflen);
#else
  SNPRINTF(strerrbuf,buflen,"%s",strerror(errnum)); /* hopefully this is thread-safe on the target system! */
#endif
}

namespace VolMagick
{
  static const VoxelOperationStatusMessenger* vosmDefault = NULL;
  void setDefaultMessenger(const VoxelOperationStatusMessenger* vosm) { vosmDefault = vosm; }

  Voxels::Voxels(const Dimension& d, VoxelType vt) 
    : _dimension(d), _voxelType(vt), _minIsSet(false), _maxIsSet(false), _vosm(vosmDefault)
  {
    try
      {
	_voxels.reset(new unsigned char[XDim()*YDim()*ZDim()*voxelSize()]);
	memset(_voxels.get(),0,XDim()*YDim()*ZDim()*voxelSize());
      }
    catch(std::bad_alloc& e)
      {
	throw MemoryAllocationError("Could not allocate memory for voxels!");
      }
  }

  Voxels::Voxels(const void *v, const Dimension& d, VoxelType vt)
    : _dimension(d), _voxelType(vt), _minIsSet(false), _maxIsSet(false), _vosm(vosmDefault)
  {
    try
      {
	_voxels.reset(new unsigned char[XDim()*YDim()*ZDim()*voxelSize()]);
	memcpy(_voxels.get(),v,XDim()*YDim()*ZDim()*VoxelTypeSizes[_voxelType]);
      }
    catch(std::bad_alloc& e)
      {
	throw MemoryAllocationError("Could not allocate memory for voxels!");
      }
  }

  Voxels::Voxels(const Voxels& v)
    : _dimension(v.dimension()), 
      _voxelType(v.voxelType()), _minIsSet(false), _maxIsSet(false), _vosm(vosmDefault)
  {
    _voxels = v._voxels;
    if(v.minIsSet() && v.maxIsSet())
      {
	min(v.min());
	max(v.max());
      }
  }

  Voxels::~Voxels() { }

  void Voxels::dimension(const Dimension& d)
  {
    if(d.isNull()) throw NullDimension("Null volume dimension.");

    if(dimension() == d) return;

    Voxels bak(*this); //backup voxels into bak

    //allocate for the new dimension
    try
      {
	//in case this throws...
	boost::shared_array<unsigned char> tmp(new unsigned char[d.xdim*d.ydim*d.zdim*voxelSize()]);
	_voxels = tmp;
      }
    catch(std::bad_alloc& e)
      {
	throw MemoryAllocationError("Could not allocate memory for voxels!");
      }
    
    _dimension = d;
    memset(_voxels.get(),0,XDim()*YDim()*ZDim()*voxelSize());

    //copy the voxels back
    for(uint64 k = 0; k < ZDim() && k < bak.ZDim(); k++)
      for(uint64 j = 0; j < YDim() && j < bak.YDim(); j++)
	for(uint64 i = 0; i < XDim() && i < bak.XDim(); i++)
	  (*this)(i,j,k, bak(i,j,k));
  }
  
  void Voxels::voxelType(VoxelType vt)
  {
    if(voxelType() == vt) return;

    Voxels bak(*this); // backup voxels into bak

    //allocate for the new voxel size
    try
      {
	//in case this throws...
	boost::shared_array<unsigned char> tmp(new unsigned char[XDim()*YDim()*ZDim()*VoxelTypeSizes[vt]]);
	_voxels = tmp;
      }
    catch(std::bad_alloc& e)
      {
	throw MemoryAllocationError("Could not allocate memory for voxels!");
      }

    _voxelType = vt;
    memset(_voxels.get(),0,XDim()*YDim()*ZDim()*voxelSize());

    //copy the voxels back
    uint64 len = XDim()*YDim()*ZDim();
    for(uint64 i = 0; i<len; i++)
      (*this)(i,bak(i));
  }

  void Voxels::calcMinMax() const
  {
    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::CalculatingMinMax, ZDim());
    double val;
    size_t len = XDim()*YDim()*ZDim(), i, count=0;
    if(len == 0) return;
    val = (*this)(0);
    _min = _max = val;
    for(i=0; i<len; i++)
      {
	val = (*this)(i);
	if(val < _min) _min = val;
	if(val > _max) _max = val;

	if(_vosm && (i % (XDim()*YDim())) == 0)
	  _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMinMax, count++);
      }
    _minIsSet = _maxIsSet = true;
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::CalculatingMinMax);
  }

  double Voxels::min(uint64 off_x, uint64 off_y, uint64 off_z,
		     const Dimension& dim) const
  {
    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::CalculatingMin, dim[2]);

    double val;
    uint64 i,j,k;
    val = (*this)(0,0,0);
    for(k=0; k<dim[2]; k++)
      {
	for(j=0; j<dim[1]; j++)
	  for(i=0; i<dim[0]; i++)
	    if(val > (*this)(i+off_x,j+off_y,k+off_z))
	      val = (*this)(i+off_x,j+off_y,k+off_z);
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMin, k);
      }
    
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::CalculatingMin);
    return val;
  }

  double Voxels::max(uint64 off_x, uint64 off_y, uint64 off_z,
		     const Dimension& dim) const
  {
    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::CalculatingMax, dim[2]);

    double val;
    uint64 i,j,k;
    val = (*this)(0,0,0);
    for(k=0; k<dim[2]; k++)
      {
	for(j=0; j<dim[1]; j++)
	  for(i=0; i<dim[0]; i++)
	    if(val < (*this)(i+off_x,j+off_y,k+off_z))
	      val = (*this)(i+off_x,j+off_y,k+off_z);
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMax, k);
      }
    
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::CalculatingMax);
    return val;
  }
  
  Voxels& Voxels::copy(const Voxels& vox)
  {
    //voxelType(vox.voxelType());
    //dimension(vox.dimension());
    //memcpy(_voxels,*vox,XDim()*YDim()*ZDim()*VoxelTypeSizes[voxelType()]);
    _voxelType = vox._voxelType;
    _dimension = vox._dimension;
    _voxels = vox._voxels;
    if(vox.minIsSet() && vox.maxIsSet())
      {
	min(vox.min());
	max(vox.max());
      }
    else
      unsetMinMax();
    return *this;
  }

  Voxels& Voxels::sub(uint64 off_x, uint64 off_y, uint64 off_z,
		      const Dimension& subvoldim)
  {
    if(off_x+subvoldim[0]-1 >= dimension()[0] || 
       off_y+subvoldim[1]-1 >= dimension()[1] || 
       off_z+subvoldim[2]-1 >= dimension()[2])
      throw IndexOutOfBounds("Subvolume offset and/or dimension is out of bounds");

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::SubvolumeExtraction, subvoldim[2]);

    Voxels tmp(*this); // back up this object into tmp

    dimension(subvoldim); // change this object's dimension to the subvolume dimension

    //copy the subvolume voxels
    for(uint64 k=0; k<dimension()[2]; k++)
      {
	for(uint64 j=0; j<dimension()[1]; j++)
	  for(uint64 i=0; i<dimension()[0]; i++)
	    (*this)(i,j,k,tmp(i+off_x,j+off_y,k+off_z));
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::SubvolumeExtraction, k);
      }

    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::SubvolumeExtraction);
    return *this;
  }

  Voxels& Voxels::fill(double val)
  {
    return fillsub(0,0,0,dimension(),val);
  }

  Voxels& Voxels::fillsub(uint64 off_x, uint64 off_y, uint64 off_z,
			  const Dimension& subvoldim, double val)
  {
    if(off_x+subvoldim[0]-1 >= dimension()[0] || 
       off_y+subvoldim[1]-1 >= dimension()[1] || 
       off_z+subvoldim[2]-1 >= dimension()[2])
      throw IndexOutOfBounds("Subvolume offset and/or dimension is out of bounds");

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::Fill, subvoldim[2]);

    for(uint64 k=0; k<subvoldim[2]; k++)
      {
	for(uint64 j=0; j<subvoldim[1]; j++)
	  for(uint64 i=0; i<subvoldim[0]; i++)
	    (*this)(i+off_x,j+off_y,k+off_z,val);
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::Fill, k);
      }

    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::Fill);
    return *this;
  }

  Voxels& Voxels::map(double min_, double max_)
  {
    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::Map, ZDim());
    uint64 len = XDim()*YDim()*ZDim(), count=0;
    for(uint64 i=0; i<len; i++)
      {
	(*this)(i,min_ + (((*this)(i) - min())/(max() - min()))*(max_ - min_));
	if(_vosm && (i % (XDim()*YDim())) == 0)
	  _vosm->step(this, VoxelOperationStatusMessenger::Map, count++);
      }
    min(min_); max(max_); // set the new min and max
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::Map);
    return *this;
  }

  Voxels& Voxels::resize(const Dimension& newdim)
  {
    double inSpaceX, inSpaceY, inSpaceZ;
    double val[8];
    uint64 resXIndex = 0, resYIndex = 0, resZIndex = 0;
    uint64 ValIndex[8];
    double xPosition = 0, yPosition = 0, zPosition = 0;
    double xRes = 0, yRes = 0, zRes = 0;
    uint64 i,j,k;
    double x,y,z;

    if(newdim.isNull()) throw NullDimension("Null voxels dimension.");

    if(dimension() == newdim) return *this; //nothing needs to be done

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::Resize, newdim[2]);

    Voxels newvox(newdim,voxelType());

    //we require a dimension of at least 2^3 
    if(newdim < Dimension(2,2,2)) 
      {
	//resize this object as if it was 2x2x2
	resize(Dimension(2,2,2));

	//copy it into newvox
	newvox.copy(*this);

	//change this object's dimension to the real dimension (destroying voxel values, hence the backup)
	dimension(newdim);

	for(k=0; k<ZDim(); k++)
	  for(j=0; j<YDim(); j++)
	    for(i=0; i<XDim(); i++)
	      (*this)(i,j,k,newvox(i,j,k));

	return *this;
      }

    // inSpace calculation
    inSpaceX = (double)(dimension()[0]-1)/(newdim[0]-1);
    inSpaceY = (double)(dimension()[1]-1)/(newdim[1]-1);
    inSpaceZ = (double)(dimension()[2]-1)/(newdim[2]-1);

    for(k = 0; k < newvox.ZDim(); k++)
      {
	z = double(k)*inSpaceZ;
	resZIndex = uint64(z);
	zPosition = z - uint64(z);
	zRes = 1;
	
	for(j = 0; j < newvox.YDim(); j++)
	  {
	    y = double(j)*inSpaceY;
	    resYIndex = uint64(y);
	    yPosition = y - uint64(y);
	    yRes =  1;

	    for(i = 0; i < newvox.XDim(); i++)
	      {
		x = double(i)*inSpaceX;
		resXIndex = uint64(x);
		xPosition = x - uint64(x);
		xRes = 1;

		// find index to get eight voxel values
		ValIndex[0] = resZIndex*dimension()[0]*dimension()[1] + resYIndex*dimension()[0] + resXIndex;
		ValIndex[1] = ValIndex[0] + 1;
		ValIndex[2] = resZIndex*dimension()[0]*dimension()[1] + (resYIndex+1)*dimension()[0] + resXIndex;
		ValIndex[3] = ValIndex[2] + 1;
		ValIndex[4] = (resZIndex+1)*dimension()[0]*dimension()[1] + resYIndex*dimension()[0] + resXIndex;
		ValIndex[5] = ValIndex[4] + 1;
		ValIndex[6] = (resZIndex+1)*dimension()[0]*dimension()[1] + (resYIndex+1)*dimension()[0] + resXIndex;
		ValIndex[7] = ValIndex[6] + 1;

		if(resXIndex>=dimension()[0]-1)
		  {
		    ValIndex[1] = ValIndex[0];
		    ValIndex[3] = ValIndex[2];
		    ValIndex[5] = ValIndex[4];
		    ValIndex[7] = ValIndex[6];
		  }
		if(resYIndex>=dimension()[1]-1)
		  {
		    ValIndex[2] = ValIndex[0];
		    ValIndex[3] = ValIndex[1];
		    ValIndex[6] = ValIndex[4];
		    ValIndex[7] = ValIndex[5];
		  }
		if(resZIndex>=dimension()[2]-1) 
		  {
		    ValIndex[4] = ValIndex[0];
		    ValIndex[5] = ValIndex[1];
		    ValIndex[6] = ValIndex[2];
		    ValIndex[7] = ValIndex[3];
		  }

		for(int Index = 0; Index < 8; Index++) 
		  val[Index] = (*this)(ValIndex[Index]);
		  
		newvox(i,j,k,
		       getTriVal(val, xPosition, yPosition, zPosition, xRes, yRes, zRes));
	      }
	  }

	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::Resize, k);
      }

    copy(newvox); //make this into a copy of the interpolated voxels

    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::Resize);

    return *this;
  }

  Voxels& Voxels::composite(const Voxels& compVox, int64 off_x, int64 off_y, int64 off_z, const CompositeFunction& func)
  {
    uint64 i,j,k;

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::Composite, compVox.ZDim());

    for(k=0; k<compVox.ZDim(); k++)
      {
	for(j=0; j<compVox.YDim(); j++)
	  for(i=0; i<compVox.XDim(); i++)
	    if((int64(i)+off_x >= 0) && (int64(i)+off_x < int64(XDim())) &&
	       (int64(j)+off_y >= 0) && (int64(j)+off_y < int64(YDim())) &&
	       (int64(k)+off_z >= 0) && (int64(k)+off_z < int64(ZDim())))
	      (*this)(int64(i) + off_x, int64(j) + off_y, int64(k) + off_z,
		      func(compVox,i,j,k,
			   *this, int64(i) + off_x, int64(j) + off_y, int64(k) + off_z));
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::Composite, k);
      }

    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::Composite);
    
    return *this;
  }

  Volume& Volume::copy(const Volume& vol)
  {
    Voxels::copy(vol);
    boundingBox(vol.boundingBox());
    desc(vol.desc());
    return *this;
  }

  Volume& Volume::sub(uint64 off_x, uint64 off_y, uint64 off_z,
		      const Dimension& subvoldim
#ifdef _MSC_VER 
			, int brain_damage //avoiding VC++ error C2555
#endif
			  )
  {
    Voxels::sub(off_x,off_y,off_z,subvoldim);
    boundingBox().setMin(XMin()+XSpan()*off_x,
			 YMin()+YSpan()*off_y,
			 ZMin()+ZSpan()*off_z);
    boundingBox().setMax(XMin()+XSpan()*(off_x+XDim()-1),
			 YMin()+YSpan()*(off_y+YDim()-1),
			 ZMin()+YSpan()*(off_z+ZDim()-1));
    return *this;
  }

#if 0
  Volume& Volume::compositeObj(const Volume& compVol, double off_x, double off_y, double off_z, const CompositeFunction& func)
  {
    return *this; // finish me!
  }
#endif

  Volume& Volume::sub(const BoundingBox& subvolbox)
  {
    //keep the span of the subvolume as close as possible to the original volume
    return *this; //finish me!
  }

  Volume& Volume::sub(const BoundingBox& subvolbox, const Dimension& subvoldim)
  {
    if(!subvolbox.isWithin(boundingBox()))
      throw SubVolumeOutOfBounds("Subvolume bounding box must be within the bounding box of the original volume.");

    Volume subvol(subvoldim,
		  voxelType(),
		  subvolbox);

    subvol.desc(desc());

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::SubvolumeExtraction, subvoldim[2]);

    uint64 i,j,k;
    for(k=0; k<subvol.ZDim(); k++)
      {
	for(j=0; j<subvol.YDim(); j++)
	  for(i=0; i<subvol.XDim(); i++)
	    subvol(i,j,k, this->interpolate(subvolbox.minx + i*subvol.XSpan(),
					    subvolbox.miny + j*subvol.YSpan(),
					    subvolbox.minz + k*subvol.ZSpan()));
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::SubvolumeExtraction, k);
      }

    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::SubvolumeExtraction);

    copy(subvol);
    return *this;
  }

  double Volume::interpolate(double obj_x, double obj_y, double obj_z) const
  {
    double inSpaceX, inSpaceY, inSpaceZ;
    double val[8];
    uint64 resXIndex = 0, resYIndex = 0, resZIndex = 0;
    uint64 ValIndex[8];
    double xPosition = 0, yPosition = 0, zPosition = 0;
    double xRes = 1, yRes = 1, zRes = 1;
    uint64 i,j,k;
    double x,y,z;

    if(!boundingBox().contains(obj_x,obj_y,obj_z)) 
      throw IndexOutOfBounds("Coordinates are outside of bounding box");

    resXIndex = uint64((obj_x - XMin())/XSpan());
    resYIndex = uint64((obj_y - YMin())/YSpan());
    resZIndex = uint64((obj_z - ZMin())/ZSpan());

    // find index to get eight voxel values
    ValIndex[0] = resZIndex*dimension()[0]*dimension()[1] + resYIndex*dimension()[0] + resXIndex;
    ValIndex[1] = ValIndex[0] + 1;
    ValIndex[2] = resZIndex*dimension()[0]*dimension()[1] + (resYIndex+1)*dimension()[0] + resXIndex;
    ValIndex[3] = ValIndex[2] + 1;
    ValIndex[4] = (resZIndex+1)*dimension()[0]*dimension()[1] + resYIndex*dimension()[0] + resXIndex;
    ValIndex[5] = ValIndex[4] + 1;
    ValIndex[6] = (resZIndex+1)*dimension()[0]*dimension()[1] + (resYIndex+1)*dimension()[0] + resXIndex;
    ValIndex[7] = ValIndex[6] + 1;
    
    if(resXIndex>=dimension()[0]-1)
      {
	ValIndex[1] = ValIndex[0];
	ValIndex[3] = ValIndex[2];
	ValIndex[5] = ValIndex[4];
	ValIndex[7] = ValIndex[6];
      }
    if(resYIndex>=dimension()[1]-1)
      {
	ValIndex[2] = ValIndex[0];
	ValIndex[3] = ValIndex[1];
	ValIndex[6] = ValIndex[4];
	ValIndex[7] = ValIndex[5];
      }
    if(resZIndex>=dimension()[2]-1) 
      {
	ValIndex[4] = ValIndex[0];
	ValIndex[5] = ValIndex[1];
	ValIndex[6] = ValIndex[2];
	ValIndex[7] = ValIndex[3];
      }
    
    for(int Index = 0; Index < 8; Index++) 
      val[Index] = (*this)(ValIndex[Index]);

    xPosition = obj_x - (double(resXIndex)*XSpan() + XMin());
    yPosition = obj_y - (double(resYIndex)*YSpan() + YMin());
    zPosition = obj_z - (double(resZIndex)*ZSpan() + ZMin());

    xRes = XSpan();
    yRes = YSpan();
    zRes = ZSpan();

    return getTriVal(val, xPosition, yPosition, zPosition, xRes, yRes, zRes);
  }

  Volume& Volume::combineWith(const Volume& vol, const Dimension& dim)
  {
    BoundingBox combbox = boundingBox() + vol.boundingBox();
    Volume combvol(dim,voxelType(),combbox);

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::CombineWith, dim[2]);

    for(uint64 k = 0; k < combvol.ZDim(); k++)
      {
	for(uint64 j = 0; j < combvol.YDim(); j++)
	  for(uint64 i = 0; i < combvol.XDim(); i++)
	    {
	      double x = i*combvol.XSpan() + XMin();
	      double y = j*combvol.YSpan() + YMin();
	      double z = k*combvol.ZSpan() + ZMin();
	      
	      //TODO: consider using composite algoritms or average these or something...
	      if(vol.boundingBox().contains(x,y,z))
		combvol(i,j,k, vol.interpolate(x,y,z));
	      else if(boundingBox().contains(x,y,z))
		combvol(i,j,k, interpolate(x,y,z));
	    }

	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CombineWith, k);
      }

    (*this) = combvol;

    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::CombineWith);

    return (*this);
  }

  Volume& Volume::combineWith(const Volume& vol)
  {
    return combineWith(vol,dimension());
  }

  void VolumeFileInfo::read(const std::string& file)
  {
    if(file.rfind(".rawiv") == file.size() - std::string(".rawiv").size()) // check extension
      readRawIV(file);
    else if(file.rfind(".rawv") == file.size() - std::string(".rawv").size()) // check extension
      readRawV(file);
    else if(file.rfind(".mrc") == file.size() - std::string(".mrc").size()) // check extension
      readMRC(file);
    else if(file.rfind(".inr") == file.size() - std::string(".inr").size()) // check extension
      readINR(file);
    else if(file.rfind(".vol") == file.size() - strlen(".vol") ||
	    file.rfind(".xmp") == file.size() - strlen(".xmp") ||
	    file.rfind(".spi") == file.size() - strlen(".spi"))
      readSpider(file);
    else
      throw UnsupportedVolumeFileType("VolMagick::VolumeFileInfo::read(): Cannot read " + file);

    _filename = file;
  }

  void VolumeFileInfo::calcMinMax(unsigned int var, unsigned int time) const
  {
    Volume vol;
    readVolumeFile(vol,filename(),var,time);
    _min[var][time] = vol.min();
    _max[var][time] = vol.max();
    _minIsSet[var][time] = true;
    _maxIsSet[var][time] = true;
  }

  /*
    I/O functions!
  */

  void readRawIV(Volume& vol,
		 const std::string& filename, 
		 unsigned int var, unsigned int time,
		 uint64 off_x, uint64 off_y, uint64 off_z,
		 const Dimension& subvoldim);

  void createRawIV(const std::string& filename,
		   const BoundingBox& boundingBox,
		   const Dimension& dimension,
		   const std::vector<VoxelType>& voxelTypes,
		   unsigned int numVariables, unsigned int numTimesteps,
		   double min_time, double max_time);

  void writeRawIV(const Volume& wvol, 
		  const std::string& filename,
		  unsigned int var, unsigned int time,
		  uint64 off_x, uint64 off_y, uint64 off_z);

  void readRawV(Volume& vol,
		const std::string& filename, 
		unsigned int var, unsigned int time,
		uint64 off_x, uint64 off_y, uint64 off_z,
		const Dimension& subvoldim);

  void createRawV(const std::string& filename,
		  const BoundingBox& boundingBox,
		  const Dimension& dimension,
		  const std::vector<VoxelType>& voxelTypes,
		  unsigned int numVariables, unsigned int numTimesteps,
		  double min_time, double max_time);

  void writeRawV(const Volume& wvol, 
		 const std::string& filename,
		 unsigned int var, unsigned int time,
		 uint64 off_x, uint64 off_y, uint64 off_z);

  void readMRC(Volume& vol,
	       const std::string& filename, 
	       unsigned int var, unsigned int time,
	       uint64 off_x, uint64 off_y, uint64 off_z,
	       const Dimension& subvoldim);

  void createMRC(const std::string& filename,
		 const BoundingBox& boundingBox,
		 const Dimension& dimension,
		 const std::vector<VoxelType>& voxelTypes,
		 unsigned int numVariables, unsigned int numTimesteps,
		 double min_time, double max_time);
  
  void writeMRC(const Volume& wvol, 
		const std::string& filename,
		unsigned int var, unsigned int time,
		uint64 off_x, uint64 off_y, uint64 off_z);

  void readINR(Volume& vol,
	       const std::string& filename, 
	       unsigned int var, unsigned int time,
	       uint64 off_x, uint64 off_y, uint64 off_z,
	       const Dimension& subvoldim);

  void createINR(const std::string& filename,
		 const BoundingBox& boundingBox,
		 const Dimension& dimension,
		 const std::vector<VoxelType>& voxelTypes,
		 unsigned int numVariables, unsigned int numTimesteps,
		 double min_time, double max_time);
  
  void writeINR(const Volume& wvol, 
		const std::string& filename,
		unsigned int var, unsigned int time,
		uint64 off_x, uint64 off_y, uint64 off_z);
  
  void readSpider(Volume& vol,
	       const std::string& filename, 
	       unsigned int var, unsigned int time,
	       uint64 off_x, uint64 off_y, uint64 off_z,
	       const Dimension& subvoldim);

  void createSpider(const std::string& filename,
		 const BoundingBox& boundingBox,
		 const Dimension& dimension,
		 const std::vector<VoxelType>& voxelTypes,
		 unsigned int numVariables, unsigned int numTimesteps,
		 double min_time, double max_time);
  
  void writeSpider(const Volume& wvol, 
		const std::string& filename,
		unsigned int var, unsigned int time,
		uint64 off_x, uint64 off_y, uint64 off_z);

  void readVolumeFile(Volume& vol, 
		      const std::string& filename,
		      unsigned int var, unsigned int time)
  {
    VolumeFileInfo volinfo(filename);
    readVolumeFile(vol,filename,var,time,0,0,0,volinfo.dimension());
  }

  void readVolumeFile(Volume& vol,
		      const std::string& filename, 
		      unsigned int var, unsigned int time,
		      uint64 off_x, uint64 off_y, uint64 off_z,
		      const Dimension& subvoldim)
  {
    if(filename.rfind(".rawiv") == filename.size() - strlen(".rawiv")/*std::string(".rawiv").size()*/) // check extension
      readRawIV(vol,filename,var,time,off_x,off_y,off_z,subvoldim);
    else if(filename.rfind(".rawv") == filename.size() - strlen(".rawv")/*std::string(".rawv").size()*/) // check extension
      readRawV(vol,filename,var,time,off_x,off_y,off_z,subvoldim);
    else if(filename.rfind(".mrc") == filename.size() - strlen(".mrc")/*std::string(".mrc").size()*/)
      readMRC(vol,filename,var,time,off_x,off_y,off_z,subvoldim);
    else if(filename.rfind(".inr") == filename.size() - strlen(".inr")/*std::string(".inr").size()*/)
      readINR(vol,filename,var,time,off_x,off_y,off_z,subvoldim);
    else if(filename.rfind(".vol") == filename.size() - strlen(".vol") ||
	    filename.rfind(".xmp") == filename.size() - strlen(".xmp") ||
	    filename.rfind(".spi") == filename.size() - strlen(".spi"))
      readSpider(vol,filename,var,time,off_x,off_y,off_z,subvoldim);
    else
      throw UnsupportedVolumeFileType("VolMagick::readVolumeFile(): Cannot read " + filename);
  }

  void readVolumeFile(std::vector<Volume>& vols,
		      const std::string& filename)
  {
    VolumeFileInfo volinfo(filename);
    Volume vol;
    vols.clear();
    for(unsigned int var=0; var<volinfo.numVariables(); var++)
      for(unsigned int time=0; time<volinfo.numTimesteps(); time++)
	{
	  readVolumeFile(vol,filename,var,time);
	  vols.push_back(vol);
	}
  }

  void writeVolumeFile(const Volume& vol, 
		       const std::string& filename,
		       unsigned int var, unsigned int time,
		       uint64 off_x, uint64 off_y, uint64 off_z)
  {
    if(filename.rfind(".rawiv") == filename.size() - strlen(".rawiv")/*std::string(".rawiv").size()*/) // check extension
      writeRawIV(vol,filename,var,time,off_x,off_y,off_z);
    else if(filename.rfind(".rawv") == filename.size() - strlen(".rawv")/*std::string(".rawv").size()*/) // check extension
      writeRawV(vol,filename,var,time,off_x,off_y,off_z);
    else if(filename.rfind(".mrc") == filename.size() - strlen(".mrc")/*std::string(".mrc").size()*/)
      writeMRC(vol,filename,var,time,off_x,off_y,off_z);
    else if(filename.rfind(".inr") == filename.size() - strlen(".inr")/*std::string(".inr").size()*/)
      writeINR(vol,filename,var,time,off_x,off_y,off_z);
    else if(filename.rfind(".vol") == filename.size() - strlen(".vol") ||
	    filename.rfind(".xmp") == filename.size() - strlen(".xmp") ||
	    filename.rfind(".spi") == filename.size() - strlen(".spi"))
      writeSpider(vol,filename,var,time,off_x,off_y,off_z);
    else
      throw UnsupportedVolumeFileType("VolMagick::writeVolumeFile(): Cannot write " + filename);
  }

  void writeVolumeFile(const std::vector<Volume>& vols,
		       const std::string& filename)
  {
    uint64 i;

    if(vols.size() == 0) return;
    
    //create types vector
    std::vector<VoxelType> voxelTypes;
    for(i=0; i<vols.size(); i++) voxelTypes.push_back(vols[i].voxelType());

    //create the file and write the volume info
    createVolumeFile(filename,vols[0].boundingBox(),vols[0].dimension(),voxelTypes,vols.size());
    for(i=0; i<vols.size(); i++)
      writeVolumeFile(vols[i],filename,i);
  }

  void createVolumeFile(const std::string& filename,
			const BoundingBox& boundingBox,
			const Dimension& dimension,
			const std::vector<VoxelType>& voxelTypes,
			unsigned int numVariables, unsigned int numTimesteps,
			double min_time, double max_time)
  {
    if(filename.rfind(".rawiv") == filename.size() - std::string(".rawiv").size())
      createRawIV(filename,boundingBox,dimension,voxelTypes,numVariables,numTimesteps,min_time,max_time);
    else if(filename.rfind(".rawv") == filename.size() - std::string(".rawv").size())
      createRawV(filename,boundingBox,dimension,voxelTypes,numVariables,numTimesteps,min_time,max_time);
    else if(filename.rfind(".mrc") == filename.size() - std::string(".mrc").size())
      createMRC(filename,boundingBox,dimension,voxelTypes,numVariables,numTimesteps,min_time,max_time);
    else if(filename.rfind(".inr") == filename.size() - std::string(".inr").size())
      createINR(filename,boundingBox,dimension,voxelTypes,numVariables,numTimesteps,min_time,max_time);
    else if(filename.rfind(".vol") == filename.size() - strlen(".vol") ||
	    filename.rfind(".xmp") == filename.size() - strlen(".xmp") ||
	    filename.rfind(".spi") == filename.size() - strlen(".spi"))
      createSpider(filename,boundingBox,dimension,voxelTypes,numVariables,numTimesteps,min_time,max_time);
    else
      throw UnsupportedVolumeFileType("VolMagick::createVolumeFile(): Cannot create " + filename);
  }

  void calcGradient(std::vector<Volume>& grad, const Volume& vol, VoxelType vt)
  {
    double dx,dy,dz,length;
    int i, j, k;

    if(vol.messenger()) vol.messenger()->start(&vol, VoxelOperationStatusMessenger::CalcGradient, vol.ZDim());

    Volume gradx(vol.dimension(),vt,vol.boundingBox());
    Volume grady(vol.dimension(),vt,vol.boundingBox());
    Volume gradz(vol.dimension(),vt,vol.boundingBox());

    //central differences algorithm
    for(k=0; k<int(vol.ZDim()); k++)
      {
	for(j=0; j<int(vol.YDim()); j++)
	  for(i=0; i<int(vol.XDim()); i++)
	    {
	      dx = (vol(MIN(i+1,int(vol.XDim())-1),j,k) - vol(MAX(i-1,0),j,k))/2.0;
	      dy = (vol(i,MIN(j+1,int(vol.YDim())-1),k) - vol(i,MAX(j-1,0),k))/2.0;
	      dz = (vol(i,j,MIN(k+1,int(vol.ZDim())-1)) - vol(i,j,MAX(k-1,0)))/2.0;
	      length = sqrt(dx*dx+dy*dy+dz*dz);
	      if(length>0.0)
		{
		  dx /= length;
		  dy /= length;
		  dz /= length;
		}
	      
	      switch(vt)
		{
		case UChar:
		  dx = dx*double((~char(0))>>1)+double((~char(0))>>1);
		  dy = dy*double((~char(0))>>1)+double((~char(0))>>1);
		  dz = dz*double((~char(0))>>1)+double((~char(0))>>1);
		  break;
		case UShort:
		  dx = dx*double((~short(0))>>1)+double((~short(0))>>1);
		  dy = dy*double((~short(0))>>1)+double((~short(0))>>1);
		  dz = dz*double((~short(0))>>1)+double((~short(0))>>1);
		  break;
		case UInt:
		  dx = dx*double((~int(0))>>1)+double((~int(0))>>1);
		  dy = dy*double((~int(0))>>1)+double((~int(0))>>1);
		  dz = dz*double((~int(0))>>1)+double((~int(0))>>1);
		  break;
		default: break;
		}
	      
	      gradx(i,j,k, dx);
	      grady(i,j,k, dy);
	      gradz(i,j,k, dz);
	    }

	if(vol.messenger()) vol.messenger()->step(&vol, VoxelOperationStatusMessenger::CalcGradient, k);
      }

    grad.clear();
    grad.push_back(gradx);
    grad.push_back(grady);
    grad.push_back(gradz);

    if(vol.messenger()) vol.messenger()->end(&vol, VoxelOperationStatusMessenger::CalcGradient);
  }

  void sub(Volume& dest, const Volume& vol, 
	   uint64 off_x, uint64 off_y, uint64 off_z,
	   const Dimension& subvoldim)
  {
    if(!(Dimension(off_x+subvoldim[0],off_y+subvoldim[1],off_z+subvoldim[2]) <= vol.dimension()))
      throw IndexOutOfBounds("Subvolume offset and dimension exceeds the boundary of input volume.");
    
    dest.unsetMinMax();
    dest.dimension(subvoldim);
    dest.voxelType(vol.voxelType());
    dest.boundingBox(BoundingBox(vol.XMin()+off_x*vol.XSpan(),
				 vol.YMin()+off_y*vol.YSpan(),
				 vol.ZMin()+off_z*vol.ZSpan(),
				 vol.XMin()+(off_x+subvoldim[0]-1)*vol.XSpan(),
				 vol.YMin()+(off_y+subvoldim[1]-1)*vol.YSpan(),
				 vol.ZMin()+(off_z+subvoldim[2]-1)*vol.ZSpan()));

    for(uint64 k=0; k<subvoldim[2]; k++)
      for(uint64 j=0; j<subvoldim[1]; j++)
	for(uint64 i=0; i<subvoldim[0]; i++)
	  dest(i,j,k, vol(i+off_x,j+off_y,k+off_z));
  }
};
