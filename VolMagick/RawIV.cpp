/* $Id: RawIV.cpp,v 1.4 2008/08/15 21:53:04 transfix Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "VolMagick.h"
#include "endians.h"

#ifdef __WINDOWS__ 
#define SNPRINTF _snprintf
#define FSEEK fseek
#else
#define SNPRINTF snprintf
#define FSEEK fseeko
#endif

static inline void geterrstr(int errnum, char *strerrbuf, size_t buflen)
{
#ifdef HAVE_STRERROR_R
  strerror_r(errnum,strerrbuf,buflen);
#else
  SNPRINTF(strerrbuf,buflen,"%s",strerror(errnum)); /* hopefully this is thread-safe on the target system! */
#endif
}

struct RawIVHeader
{
  float min[3];
  float max[3];
  unsigned int numVerts;
  unsigned int numCells;
  unsigned int dim[3];
  float origin[3];
  float span[3];
} rawivHeader;

namespace VolMagick
{
  void VolumeFileInfo::readRawIV(const std::string& file)
  {
    char buf[256];
    RawIVHeader rawivHeader;

    FILE *input;
    size_t i;
    VoxelType vt;

    memset(buf,0,256);

    if((input = fopen(file.c_str(),"rb")) == NULL)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error opening file '" + file + "': " + buf;
	throw ReadError(errStr);
      }

    if(fread(&rawivHeader, sizeof(rawivHeader), 1, input) != 1)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error reading file '" + file + "': " + buf;
        fclose(input);
        throw ReadError(errStr);
      }

    if(!big_endian())
      {
        for(i=0; i<3; i++) SWAP_32(&(rawivHeader.min[i]));
        for(i=0; i<3; i++) SWAP_32(&(rawivHeader.max[i]));
        SWAP_32(&(rawivHeader.numVerts));
        SWAP_32(&(rawivHeader.numCells));
        for(i=0; i<3; i++) SWAP_32(&(rawivHeader.dim[i]));
        for(i=0; i<3; i++) SWAP_32(&(rawivHeader.origin[i]));
        for(i=0; i<3; i++) SWAP_32(&(rawivHeader.span[i]));
      }

    struct stat s;
    if(stat(file.c_str(), &s)==-1)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error stat-ing file '" + file + "': " + buf;
	fclose(input);
	throw ReadError(errStr);
      }

    /* error checking */
    if(rawivHeader.dim[0] == 0 || rawivHeader.dim[1] == 0 || rawivHeader.dim[2] == 0)
      {
	fclose(input);
	throw InvalidRawIVHeader("Dimension values cannot be zero");
      }
    if(rawivHeader.numVerts != (rawivHeader.dim[0]*rawivHeader.dim[1]*rawivHeader.dim[2]))
      {
	fclose(input);
	throw InvalidRawIVHeader("Number of vertices does not match calculated value");
      }
    if(rawivHeader.numCells != ((rawivHeader.dim[0]-1)*(rawivHeader.dim[1]-1)*(rawivHeader.dim[2]-1)))
      {
	fclose(input);
	throw InvalidRawIVHeader("Number of cells does not match calculated value");
      }
    for(i=0; i<3; i++)
      if(fabs(rawivHeader.span[i] - ((rawivHeader.max[i] - rawivHeader.min[i])/(rawivHeader.dim[i] - 1))) > 0.01f)
	{
	  SNPRINTF(buf,256,"Span value (%f) does not match calculation (%f)",
		   rawivHeader.span[i],((rawivHeader.max[i] - rawivHeader.min[i])/(rawivHeader.dim[i] - 1)));
	  fclose(input);
	  throw InvalidRawIVHeader(buf);
	}

    /* get voxel type */
    uint64 voxelTypeSize = uint64(uint64(s.st_size-68)/(uint64(rawivHeader.dim[0])*uint64(rawivHeader.dim[1])*uint64(rawivHeader.dim[2])));
    switch(voxelTypeSize)
      {
      case 1: vt = UChar; break;
      case 2: vt = UShort; break;
      case 4: vt = Float; break;
      case 8: vt = Double; break;
      default:
	fclose(input);
	SNPRINTF(buf,256,"Cannot determine voxel type (size %lld)",voxelTypeSize);
	throw InvalidRawIVHeader(buf);
      }

    if(uint64(s.st_size-68) != (uint64(rawivHeader.dim[0])*uint64(rawivHeader.dim[1])*uint64(rawivHeader.dim[2])*uint64(VoxelTypeSizes[vt])))
      {
	SNPRINTF(buf,256,"Volume dimensions do not match filesize. (file size: %lld, calculated size: %lld)",
		 uint64(s.st_size),(uint64(rawivHeader.dim[0])*uint64(rawivHeader.dim[1])*uint64(rawivHeader.dim[2])*uint64(VoxelTypeSizes[vt]) + 68));
	fclose(input);
	throw InvalidRawIVHeader(buf);
      }

    /**** at this point, the volume header has no errors, so we may start modifiying this object ****/
    _numVariables = 1;
    _numTimesteps = 1;
    _tmin = _tmax = 0.0;
    
    _dimension = Dimension(rawivHeader.dim);
    _boundingBox = BoundingBox(rawivHeader.min[0],rawivHeader.min[1],rawivHeader.min[2],
			       rawivHeader.max[0],rawivHeader.max[1],rawivHeader.max[2]);

    _voxelTypes.clear();
    _voxelTypes.push_back(vt);
    _names.clear();
    _names.push_back("No Name");

    /* new volume, so min/max is now unset */
    
    _minIsSet.clear();
    _minIsSet.resize(_numVariables); for(i=0; i<_minIsSet.size(); i++) _minIsSet[i].resize(_numTimesteps);
    _min.clear();
    _min.resize(_numVariables); for(i=0; i<_min.size(); i++) _min[i].resize(_numTimesteps);
    _maxIsSet.clear();
    _maxIsSet.resize(_numVariables); for(i=0; i<_maxIsSet.size(); i++) _maxIsSet[i].resize(_numTimesteps);
    _max.clear();
    _max.resize(_numVariables); for(i=0; i<_max.size(); i++) _max[i].resize(_numTimesteps);
    
    /* 
       Check if min/max is available.

       Here we use an extension to the rawiv format if it's available.  Since the origin values are not really
       used, we can use them to encode the min and max density values.  If the origin X value has the bytes
       0xBAADBEEF, then origin y is the minimum density value and origin z is the maximum.
    */
    unsigned int *tmp_int = (unsigned int*)(&(rawivHeader.origin[0])); /* use a pointer so the compiler doesn't try to convert the float value to int */
    if(*tmp_int == 0xBAADBEEF)
      {
	_min[0][0] = rawivHeader.origin[1];
	_max[0][0] = rawivHeader.origin[2];
	_minIsSet[0][0] = _maxIsSet[0][0] = true;
      }

    fclose(input);
  }

  void readRawIV(Volume& vol,
		 const std::string& filename, 
		 unsigned int var, unsigned int time,
		 uint64 off_x, uint64 off_y, uint64 off_z,
		 const Dimension& subvoldim)
  {
    char buf[256];
    RawIVHeader rawivHeader;
    
    FILE *input;
    size_t i,j,k;
    VoxelType vt;

    memset(buf,0,256);

    if(var > 0)
      throw IndexOutOfBounds("Variable index out of bounds.");
    if(time > 0)
      throw IndexOutOfBounds("Timestep index out of bounds.");
    if(subvoldim.isNull())
      throw IndexOutOfBounds("Specified subvolume dimension is null.");

    if((input = fopen(filename.c_str(),"rb")) == NULL)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error opening file '" + filename + "': " + buf;
	throw ReadError(errStr);
      }

    if(fread(&rawivHeader, sizeof(rawivHeader), 1, input) != 1)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error reading file '" + filename + "': " + buf;
        fclose(input);
        throw ReadError(errStr);
      }

    if(!big_endian())
      {
        for(i=0; i<3; i++) SWAP_32(&(rawivHeader.min[i]));
        for(i=0; i<3; i++) SWAP_32(&(rawivHeader.max[i]));
        SWAP_32(&(rawivHeader.numVerts));
        SWAP_32(&(rawivHeader.numCells));
        for(i=0; i<3; i++) SWAP_32(&(rawivHeader.dim[i]));
        for(i=0; i<3; i++) SWAP_32(&(rawivHeader.origin[i]));
        for(i=0; i<3; i++) SWAP_32(&(rawivHeader.span[i]));
      }

    struct stat s;
    stat(filename.c_str(), &s);

    /* error checking */
    if(rawivHeader.dim[0] == 0 || rawivHeader.dim[1] == 0 || rawivHeader.dim[2] == 0)
      {
	fclose(input);
	throw InvalidRawIVHeader("Dimension values cannot be zero");
      }
    if(rawivHeader.numVerts != (rawivHeader.dim[0]*rawivHeader.dim[1]*rawivHeader.dim[2]))
      {
	fclose(input);
	throw InvalidRawIVHeader("Number of vertices does not match calculated value");
      }
    if(rawivHeader.numCells != ((rawivHeader.dim[0]-1)*(rawivHeader.dim[1]-1)*(rawivHeader.dim[2]-1)))
      {
	fclose(input);
	throw InvalidRawIVHeader("Number of cells does not match calculated value");
      }
    for(i=0; i<3; i++)
      if(fabs(rawivHeader.span[i] - ((rawivHeader.max[i] - rawivHeader.min[i])/(rawivHeader.dim[i] - 1))) > 0.01f)
	{
	  SNPRINTF(buf,256,"Span value (%f) does not match calculation (%f)",
		   rawivHeader.span[i],((rawivHeader.max[i] - rawivHeader.min[i])/(rawivHeader.dim[i] - 1)));
	  fclose(input);
	  throw InvalidRawIVHeader(buf);
	}

    /* get voxel type */
    uint64 voxelTypeSize = uint64(uint64(s.st_size-68)/(uint64(rawivHeader.dim[0])*uint64(rawivHeader.dim[1])*uint64(rawivHeader.dim[2])));
    switch(voxelTypeSize)
      {
      case 1: vt = UChar; break;
      case 2: vt = UShort; break;
      case 4: vt = Float; break;
      case 8: vt = Double; break;
      default:
	fclose(input);
	SNPRINTF(buf,256,"Cannot determine voxel type (size %lld)",voxelTypeSize);
	throw InvalidRawIVHeader(buf);
      }

    if(uint64(s.st_size-68) != (uint64(rawivHeader.dim[0])*uint64(rawivHeader.dim[1])*uint64(rawivHeader.dim[2])*uint64(VoxelTypeSizes[vt])))
      {
	SNPRINTF(buf,256,"Volume dimensions do not match filesize. (file size: %lld, calculated size: %lld)",
		 uint64(s.st_size),(uint64(rawivHeader.dim[0])*uint64(rawivHeader.dim[1])*uint64(rawivHeader.dim[2])*uint64(VoxelTypeSizes[vt]) + 68));
	fclose(input);
	throw InvalidRawIVHeader(buf);
      }

    if((off_x + subvoldim[0] - 1 >= rawivHeader.dim[0]) ||
       (off_y + subvoldim[1] - 1 >= rawivHeader.dim[1]) ||
       (off_z + subvoldim[2] - 1 >= rawivHeader.dim[2]))
      {
	fclose(input);
	throw IndexOutOfBounds("Subvolume specified is outside volume dimensions");
      }

    /**** all errors have been checked for, at this point the rawiv file should be valid... we may now write to 'vol' ****/
    try
      {
	vol.voxelType(vt);

	/* 
	   Check if min/max is available.
	   
	   Here we use an extension to the rawiv format if it's available.  Since the origin values are not really
	   used, we can use them to encode the min and max density values.  If the origin X value has the bytes
	   0xBAADBEEF, then origin y is the minimum density value and origin z is the maximum.
	*/
	unsigned int *tmp_int = (unsigned int*)(&(rawivHeader.origin[0])); /* use a pointer so the compiler doesn't try to convert the float value to int */
	if(*tmp_int == 0xBAADBEEF)
	  {
	    vol.min(rawivHeader.origin[1]);
	    vol.max(rawivHeader.origin[2]);
	  }
	
	vol.dimension(subvoldim);
	BoundingBox subvolbox;
	subvolbox.setMin(rawivHeader.min[0]+rawivHeader.span[0]*off_x,
			 rawivHeader.min[1]+rawivHeader.span[1]*off_y,
			 rawivHeader.min[2]+rawivHeader.span[2]*off_z);
	subvolbox.setMax(rawivHeader.min[0]+rawivHeader.span[0]*(off_x+subvoldim[0]-1),
			 rawivHeader.min[1]+rawivHeader.span[1]*(off_y+subvoldim[1]-1),
			 rawivHeader.min[2]+rawivHeader.span[2]*(off_z+subvoldim[2]-1));
	vol.boundingBox(subvolbox);
      }
    catch(MemoryAllocationError& e)
      {
	fclose(input);
	throw e;
      }

    /*
      read the volume data
    */
    off_t file_offx, file_offy, file_offz;
    for(k=off_z; k<=(off_z+subvoldim[2]-1); k++)
      {
	file_offz = 68+k*rawivHeader.dim[0]*rawivHeader.dim[1]*vol.voxelSize();
	for(j=off_y; j<=(off_y+subvoldim[1]-1); j++)
	  {
	    file_offy = j*rawivHeader.dim[0]*vol.voxelSize();
	    file_offx = off_x*vol.voxelSize();
	    //seek and read a scanline at a time
	    if(FSEEK(input,file_offx+file_offy+file_offz,SEEK_SET) == -1)
	      {
		geterrstr(errno,buf,256);
		std::string errStr = "Error reading volume data in file '" + filename + "': " + buf;
		fclose(input);
		throw ReadError(errStr);
	      }
	    if(fread(*vol+
		     (k-off_z)*vol.XDim()*vol.YDim()*vol.voxelSize()+
		     (j-off_y)*vol.XDim()*vol.voxelSize(),
		     vol.voxelSize(),vol.XDim(),input) != vol.XDim())
	      {
		geterrstr(errno,buf,256);
		std::string errStr = "Error reading volume data in file '" + filename + "': " + buf;
		fclose(input);
		throw ReadError(errStr);
	      }
	  }
      }

    /* swap the volume data if on little endian machine */
    if(!big_endian())
      {
	size_t len = vol.XDim()*vol.YDim()*vol.ZDim();
	switch(vol.voxelType())
	  {
	  case UShort: for(i=0;i<len;i++) SWAP_16(*vol+i*vol.voxelSize()); break;
	  case Float:  for(i=0;i<len;i++) SWAP_32(*vol+i*vol.voxelSize()); break;
	  case Double: for(i=0;i<len;i++) SWAP_64(*vol+i*vol.voxelSize()); break;
	  default: break; /* no swapping needed for unsigned char data, and unsigned int is not defined for rawiv */
	  }
      }

    fclose(input);
  }

  void createRawIV(const std::string& filename,
		   const BoundingBox& boundingBox,
		   const Dimension& dimension,
		   const std::vector<VoxelType>& voxelTypes,
		   unsigned int numVariables, unsigned int numTimesteps,
		   double min_time, double max_time)
  {
    char buf[256];
    RawIVHeader rawivHeader;
     
    FILE *output;
    size_t i,j,k;

    memset(buf,0,256);

    if(boundingBox.isNull())
      throw InvalidBoundingBox("Bounding box must not be null");
    if(dimension.isNull())
      throw InvalidBoundingBox("Dimension must not be null");
    if(numVariables > 1)
      {
	char buf[256];
	sprintf(buf,"RawIV format only supports 1 variable (%d requested)",numVariables);
	throw InvalidRawIVHeader(buf);
      }
    if(numTimesteps > 1)
      {
	char buf[256];
	sprintf(buf,"RawIV format only supports 1 timestep (%d requested)",numTimesteps);
	throw InvalidRawIVHeader(buf);
      }
    if(voxelTypes.size() > 1)
      throw InvalidRawIVHeader("RawIV format only supports 1 variable and 1 timestep. (too many voxel types specified)");
    if(min_time != max_time)
      throw InvalidRawIVHeader("RawIV format does not support multiple timesteps. (min time and max time must be the same)");

    memset(&rawivHeader,0,sizeof(RawIVHeader));

    rawivHeader.min[0] = boundingBox.minx;
    rawivHeader.min[1] = boundingBox.miny;
    rawivHeader.min[2] = boundingBox.minz;
    rawivHeader.max[0] = boundingBox.maxx;
    rawivHeader.max[1] = boundingBox.maxy;
    rawivHeader.max[2] = boundingBox.maxz;
    rawivHeader.numVerts = dimension[0] * dimension[1] * dimension[2];
    rawivHeader.numCells = (dimension[0]-1) * (dimension[1]-1) * (dimension[2]-1);
    rawivHeader.dim[0] = dimension[0];
    rawivHeader.dim[1] = dimension[1];
    rawivHeader.dim[2] = dimension[2];
    rawivHeader.span[0] = (boundingBox.maxx-boundingBox.minx)/(dimension[0]-1);
    rawivHeader.span[1] = (boundingBox.maxy-boundingBox.miny)/(dimension[1]-1);
    rawivHeader.span[2] = (boundingBox.maxz-boundingBox.minz)/(dimension[2]-1);
    
    /* write extended rawiv header volumes in which min and max voxel values are encoded in the origin values */
    /* use a pointer so the compiler doesn't try to convert the float value to int */
    unsigned int *tmp_int = (unsigned int*)(&(rawivHeader.origin[0]));
    *tmp_int = 0xBAADBEEF;
    rawivHeader.origin[1] = 0.0;
    rawivHeader.origin[2] = 0.0;

    if(!big_endian())
      {
	for(i=0; i<3; i++) SWAP_32(&(rawivHeader.min[i]));
	for(i=0; i<3; i++) SWAP_32(&(rawivHeader.max[i]));
	SWAP_32(&(rawivHeader.numVerts));
	SWAP_32(&(rawivHeader.numCells));
	for(i=0; i<3; i++) SWAP_32(&(rawivHeader.dim[i]));
	for(i=0; i<3; i++) SWAP_32(&(rawivHeader.origin[i]));
	for(i=0; i<3; i++) SWAP_32(&(rawivHeader.span[i]));
      }

    if((output = fopen(filename.c_str(),"wb")) == NULL)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error opening file '" + filename + "': " + buf;
	throw WriteError(errStr);
      }

    if(fwrite(&rawivHeader,sizeof(rawivHeader),1,output) != 1)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error writing header to file '" + filename + "': " + buf;
	fclose(output);
	throw WriteError(errStr);
      }

    unsigned char * scanline = (unsigned char *)calloc(dimension[0]*VoxelTypeSizes[voxelTypes[0]],1);
    if(scanline == NULL)
      {
	fclose(output);
	throw MemoryAllocationError("Unable to allocate memory for write buffer");
      }

    // write a scanline at a time
    for(k=0; k<dimension[2]; k++)
      for(j=0; j<dimension[1]; j++)
	{
	  if(fwrite(scanline,VoxelTypeSizes[voxelTypes[0]],dimension[0],output) != dimension[0])
	    {
	      geterrstr(errno,buf,256);
	      std::string errStr = "Error writing volume data to file '" + filename + "': " + buf;
	      free(scanline);
	      fclose(output);
	      throw WriteError(errStr);
	    }
	}

    free(scanline);
    fclose(output);
  }

  void writeRawIV(const Volume& wvol, 
		  const std::string& filename,
		  unsigned int var, unsigned int time,
		  uint64 off_x, uint64 off_y, uint64 off_z)
  {
    VolumeFileInfo volinfo;
    char buf[256];
    RawIVHeader rawivHeader;
    bool creatingNewFile = false;
     
    FILE *output;
    size_t i,j,k;

    uint64 outvol_xdim, outvol_ydim, outvol_zdim;

    memset(buf,0,256);
     
    if(var > 0)
      throw IndexOutOfBounds("Variable index out of bounds.");
    if(time > 0)
      throw IndexOutOfBounds("Timestep index out of bounds.");

    Volume vol(wvol);

    //check if the file exists and we can write the specified subvolume to it
    try
      {
	volinfo.read(filename);
	//if(!(Dimension(off_x+vol.XDim(),off_y+vol.YDim(),off_z+vol.ZDim()) <= volinfo.dimension()))
	if(off_x+vol.XDim() > volinfo.dimension()[0] &&
	   off_y+vol.YDim() > volinfo.dimension()[1] &&
	   off_z+vol.ZDim() > volinfo.dimension()[2])
	  {
	    std::string errStr = "File '" + filename + "' exists but is too small to write volume at specified offset";
	    throw IndexOutOfBounds(errStr);
	  }
	vol.voxelType(volinfo.voxelType()); //change the volume's voxel type to match that of the file
      }
    catch(ReadError &e)
      {
	//create a blank file since file doesn't exist (or there was an error reading the existing file)
	BoundingBox box(vol.boundingBox());
	box.minx -= off_x * vol.XSpan();
	box.miny -= off_y * vol.YSpan();
	box.minz -= off_z * vol.ZSpan();
	Dimension dim(vol.dimension());
	dim[0] += off_x;
	dim[1] += off_y;
	dim[2] += off_z;
	createRawIV(filename,box,dim,std::vector<VoxelType>(1,vol.voxelType()),1,1,0.0,0.0);
	volinfo.read(filename);

	if(var >= volinfo.numVariables())
	  {
	    std::string errStr = "Variable index exceeds number of variables in file '" + filename + "'";
	    throw IndexOutOfBounds(errStr);
	  }
	if(time >= volinfo.numTimesteps())
	  {
	    std::string errStr = "Timestep index exceeds number of timesteps in file '" + filename + "'";
	    throw IndexOutOfBounds(errStr);
	  }

	creatingNewFile = true;
      }

    if((output = fopen(filename.c_str(),"r+b")) == NULL)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error opening file '" + filename + "': " + buf;
	throw WriteError(errStr);
      }
	
    //read the header to set the min/max values for rawiv ext.
    if(fread(&rawivHeader, sizeof(rawivHeader), 1, output) != 1)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error reading file '" + filename + "': " + buf;
        fclose(output);
        throw ReadError(errStr);
      }

    if(!big_endian())
    {
      for(i=0; i<3; i++) SWAP_32(&(rawivHeader.origin[i]));
      for(i=0; i<3; i++) SWAP_32(&(rawivHeader.dim[i]));
    }
    
    /* write extended rawiv header volumes in which min and max voxel values are encoded in the origin values */
    /* use a pointer so the compiler doesn't try to convert the float value to int */
    unsigned int *tmp_int = (unsigned int*)(&(rawivHeader.origin[0]));
    //if(*tmp_int != 0xBAADBEEF) //if the header doesn't have a min/max extension, write one
      {
	*tmp_int = 0xBAADBEEF;
	if(creatingNewFile) //the blank file is zero filled...
	  {
	    rawivHeader.origin[1] = MIN(0.0,vol.min());
	    rawivHeader.origin[2] = MAX(0.0,vol.max());
	  }
	else
	  {
	    rawivHeader.origin[1] = MIN(volinfo.min(),vol.min());
	    rawivHeader.origin[2] = MAX(volinfo.max(),vol.max());
	  }
      }

    outvol_xdim = rawivHeader.dim[0];
    outvol_ydim = rawivHeader.dim[1];
    outvol_zdim = rawivHeader.dim[2];

    if(!big_endian())
    {
      for(i=0; i<3; i++) SWAP_32(&(rawivHeader.origin[i]));
      for(i=0; i<3; i++) SWAP_32(&(rawivHeader.dim[i]));
    }
    
    if(FSEEK(output,0,SEEK_SET) == -1)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error seeking in file '" + filename + "': " + buf;
	fclose(output);
	throw ReadError(errStr);
      }

    if(fwrite(&rawivHeader,sizeof(rawivHeader),1,output) != 1)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error writing header to file '" + filename + "': " + buf;
	fclose(output);
	throw WriteError(errStr);
      }
    
    unsigned char * scanline = (unsigned char *)malloc(vol.XDim()*vol.voxelSize());
    if(scanline == NULL)
      {
	fclose(output);
	throw MemoryAllocationError("Unable to allocate memory for write buffer");
      }
	
    /*
      write the volume data
    */
    off_t file_offx, file_offy, file_offz;
    for(k=off_z; k<=(off_z+vol.ZDim()-1); k++)
      {
	file_offz = 68+k*outvol_xdim*outvol_ydim*vol.voxelSize();
	for(j=off_y; j<=(off_y+vol.YDim()-1); j++)
	  {
	    file_offy = j*outvol_xdim*vol.voxelSize();
	    file_offx = off_x*vol.voxelSize();
	    
	    //seek and write a scanline at a time
	    if(FSEEK(output,file_offx+file_offy+file_offz,SEEK_SET) == -1)
	      {
		geterrstr(errno,buf,256);
		std::string errStr = "Error seeking in file '" + filename + "': " + buf;
		fclose(output);
		throw ReadError(errStr);
	      }

	    memcpy(scanline,*vol+
		   ((k-off_z)*vol.XDim()*vol.YDim()*vol.voxelSize())+
		   ((j-off_y)*vol.XDim()*vol.voxelSize()),
		   vol.XDim()*vol.voxelSize());

	    /* swap the volume data if on little endian machine */
	    if(!big_endian())
	      {
		size_t len = vol.XDim();
		switch(vol.voxelType())
		  {
		  case UShort: for(i=0;i<len;i++) SWAP_16(scanline+i*vol.voxelSize()); break;
		  case Float:  for(i=0;i<len;i++) SWAP_32(scanline+i*vol.voxelSize()); break;
		  case Double: for(i=0;i<len;i++) SWAP_64(scanline+i*vol.voxelSize()); break;
		  default: break; /* no swapping needed for unsigned char data, and unsigned int is not defined for rawiv */
		  }
	      }

	    if(fwrite(scanline,vol.voxelSize(),vol.XDim(),output) != vol.XDim())
	      {
		geterrstr(errno,buf,256);
		std::string errStr = "Error writing volume data to file '" + filename + "': " + buf;
		free(scanline);
		fclose(output);
		throw WriteError(errStr);
	      }
	  }
      }

    free(scanline);

    //fix a bug with min/max values in header... if we use this function to replace a whole volume's set of voxels,
    //the min/max calculation above can be wrong!
    {
      Volume minmax_vol;
      readVolumeFile(minmax_vol,volinfo.filename());
      minmax_vol.unsetMinMax();
      float realminmax[2];
      realminmax[0] = minmax_vol.min();
      realminmax[1] = minmax_vol.max();
     
      ptrdiff_t offset = (unsigned char*)(&(rawivHeader.origin[1])) - (unsigned char *)(&(rawivHeader.min[0]));

      if(FSEEK(output,offset,SEEK_SET) == -1)
	{
	  geterrstr(errno,buf,256);
	  std::string errStr = "Error seeking in file '" + filename + "': " + buf;
	  fclose(output);
	  throw ReadError(errStr);
	}

      if(!big_endian())
	{
	  SWAP_32(&realminmax[0]);
	  SWAP_32(&realminmax[1]);
	}
      
      if(fwrite(realminmax,sizeof(float),2,output) != 2)
	{
	  geterrstr(errno,buf,256);
	  std::string errStr = "Error writing min/max data to file '" + filename + "': " + buf;
	  fclose(output);
	  throw WriteError(errStr);
	}
    }

    fclose(output);
  }
};
