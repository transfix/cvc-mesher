/* $Id: MRC.cpp,v 1.5 2008/08/15 21:53:04 transfix Exp $ */

#include <boost/format.hpp>
#include <boost/scoped_array.hpp>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#ifdef __SOLARIS__
#include <ieeefp.h>
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "VolMagick.h"
#include "endians.h"

///\struct MrcHeader
///\brief The header of an MRC file.
typedef struct {
	
//! # of columns ( fastest changing in the map    
	int    nx;
//! # of rows                                     
	int    ny;
//! # of sections (slowest changing in the map    
	int    nz;
	
//! data type
//! 0 = image data in bytes
//! 1 = image data in short integer
//! 2 = image data in floats
//! 3 = complex data in complex short integers
//! 4 = complex data in complex reals          
	int    mode;
	
//! number of first column in map (default = 0)   
	int    nxstart;
//! number of first row in map (default = 0)      
	int    nystart;
//! number of first ssection in map (default = 0) 
	int    nzstart;
	
//! number of intervals along X                   
	int    mx;
//! number of intervals along Y                   
	int    my;
//! number of intervals along Z                   
	int    mz;
	
//! cell dimensions in X (angstrom)               
	float  xlength;
//! cell dimensions in Y (angstrom)               
	float  ylength;
//! cell dimensions in Z (angstrom)               
	float  zlength;
	
//! cell angles between Y and Z                   
	float  alpha;
//! cell angles between X and Z                   
	float  beta;
//! cell angles between X and Y                   
	float  gamma;
	
//! number of axis corresponding to columns (X)   
	int    mapc;
//! number of axis corresponding to rows (Y)      
	int    mapr;
//! number of axis corresponding to sections (Z)  
	int    maps;
	
//! minimum density value                         
	float  amin;
//! maximum density value                         
	float  amax;
//! mean density value                            
	float  amean;
	
//! space group number (0 for images)             
	int    ispg;
//! # of bytes for symmetry operators             
	int    nsymbt;
	
//! user defined storage space                    
	int    extra[25];
	
//! X phase origin                                
	float  xorigin;
//! Y phase origin                                
	float  yorigin;
//! Z phase origin
	float  zorigin;

//! character string 'MAP '
	char   map[4];

//! machine stamp
	int    machst;

//! rms deviation of map from mean density
	float  rms;
	
//! # of labels being used in the MRC header      
	int    nlabl;
	
//! actual text labels                            
	char   label[10][80];
	
} MrcHeader;

typedef struct
{
  float aTilt;
  float bTilt;
  float xStage;
  float yStage;
  float zStage;
  float xShift;
  float yShift;
  float defocus;
  float expTime;
  float meanInt;
  float tiltAxis;
  float pixelSize;
  float imageMag;
  char filler[76];
} ExtendedMrcHeader;

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

// XXX: This is UGLY. Windows does not have this function in its math library.
#if defined(_MSC_VER)
static inline int finite(float) { return 0; }
#endif

namespace VolMagick
{
  static inline bool checkHeader(MrcHeader& header, size_t fileSize)
  {
    size_t sizes[] = { 1, 2, 4 };

    //check for details we dont support
    if(header.mode<0 || header.mode>2) return false;

    //check the fileSize
    if(header.nsymbt == 0)
      {
	if(sizes[header.mode]*header.nx*header.ny*header.nz + sizeof(MrcHeader) != fileSize)
	  return false;
      }
    else
      {
	if(sizes[header.mode]*header.nx*header.ny*header.nz + sizeof(MrcHeader) + sizeof(ExtendedMrcHeader) != fileSize)
	  return false;
      }

    return true;
  }

  static inline void swapHeader(MrcHeader& header)
  {
    SWAP_32(&(header.nx));
    SWAP_32(&(header.ny));
    SWAP_32(&(header.nz));
    SWAP_32(&(header.mode));
    SWAP_32(&(header.nxstart));
    SWAP_32(&(header.nystart));
    SWAP_32(&(header.nzstart));

    SWAP_32(&(header.mx));
    SWAP_32(&(header.my));
    SWAP_32(&(header.mz));
	
    SWAP_32(&(header.xlength));
    SWAP_32(&(header.ylength));
    SWAP_32(&(header.zlength));
	
    SWAP_32(&(header.alpha));
    SWAP_32(&(header.beta));
    SWAP_32(&(header.gamma));
	
    SWAP_32(&(header.mapc));
    SWAP_32(&(header.mapr));
    SWAP_32(&(header.maps));
	
    SWAP_32(&(header.amin));
    SWAP_32(&(header.amax));
    SWAP_32(&(header.amean));
	
    SWAP_32(&(header.ispg));
    SWAP_32(&(header.nsymbt));
	
    for(unsigned int i=0; i<25; i++) SWAP_32(&(header.extra[i]));
	
    SWAP_32(&(header.xorigin));
    SWAP_32(&(header.yorigin));
    SWAP_32(&(header.zorigin));

    SWAP_32(&(header.rms));
	
    SWAP_32(&(header.nlabl));
  }

  void VolumeFileInfo::readMRC(const std::string& file)
  {
    char buf[256];
    FILE *input;
    size_t i;
    VoxelType mrcTypes[] = { UChar, UShort, Float };
    bool mustSwap = false;

    MrcHeader header;
    ExtendedMrcHeader extheader;

    memset(buf,0,256);

    if((input = fopen(file.c_str(),"rb")) == NULL)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error opening file '" + file + "': " + buf;
	throw ReadError(errStr);
      }
    
    if(fread(&header, sizeof(MrcHeader), 1, input) != 1)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error reading file '" + file + "': " + buf;
        fclose(input);
        throw ReadError(errStr);
      }
    
    struct stat s;
    if(stat(file.c_str(), &s)==-1)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error stat-ing file '" + file + "': " + buf;
	fclose(input);
	throw ReadError(errStr);
      }

    // try to figure out the endianness of the file
    if(checkHeader(header, s.st_size))
      mustSwap = false;
    else 
      {
	// swap and try again
	swapHeader(header);
	if(checkHeader(header, s.st_size))
	  mustSwap = true;
	else 
	  {
	    // we dont support wierd or exotic endianness
	    fclose(input);
	    throw InvalidMRCFile("Cannot determine endianness");
	  }
      }

    if(!(header.map[0]=='M' && header.map[1]=='A' && header.map[2]=='P'))
      {
	//interpret old header style
	_dimension = Dimension(header.nx,header.ny,header.nz);
	_numTimesteps = 1;
	_numVariables = 1;
	_names.clear();
	_names.push_back("No Name");
	_voxelTypes.clear();
	if(header.mode > 2)
	  throw InvalidMRCFile("Unsupported datatype");
	_voxelTypes.push_back(mrcTypes[header.mode]);

	// we need to double check the meaning of xlength
	// (plus some extra paranoia)
	if (header.xlength<=0.0 || header.ylength<=0.0 || header.zlength<=0.0
	    || !finite(header.xlength) || !finite(header.ylength) || !finite(header.zlength))
	  _boundingBox = BoundingBox(0.0,0.0,0.0,
				     header.nx,header.ny,header.nz);
	else
	  _boundingBox = BoundingBox(0.0,0.0,0.0,
				     header.xlength,header.ylength,header.zlength);
      }
    else
      {
	double tmpmin[3],tmpmax[3];

	// new MRC file format
	_dimension = Dimension(header.nx,header.ny,header.nz);
	_numTimesteps = 1;
	_numVariables = 1;
	_names.clear();
	_names.push_back("No Name");
	_voxelTypes.clear();
	_voxelTypes.push_back(mrcTypes[header.mode]);

	// make sure we aren't using garbage values
	if (!finite(header.xorigin) || !finite(header.yorigin) || !finite(header.zorigin))
	  {
	    tmpmin[0] = 0.0;
	    tmpmin[1] = 0.0;
	    tmpmin[2] = 0.0;
	  }
	else 
	  {
	    tmpmin[0] = header.xorigin;
	    tmpmin[1] = header.yorigin;
	    tmpmin[2] = header.zorigin;
	  }

	// we need to double check the meaning of xlength
	// (plus some extra paranoia)
	if (header.xlength<=0.0 || header.ylength<=0.0 || header.zlength<=0.0
	    || !finite(header.xlength) || !finite(header.ylength) || !finite(header.zlength))
	  {
	    // hmm, this is wierd
	    tmpmax[0] = tmpmin[0] + header.nx;
	    tmpmax[1] = tmpmin[1] + header.ny;
	    tmpmax[2] = tmpmin[2] + header.nz;
	  }
	else 
	  {
	    tmpmax[0] = tmpmin[0] + header.xlength;
	    tmpmax[1] = tmpmin[1] + header.ylength;
	    tmpmax[2] = tmpmin[2] + header.zlength;
	  }

	_boundingBox = BoundingBox(tmpmin[0],tmpmin[1],tmpmin[2],
				   tmpmax[0],tmpmax[1],tmpmax[2]);
      }

    _tmin = _tmax = 0.0;

    /* new volume, so min/max is now unset */
    _minIsSet.clear();
    _minIsSet.resize(_numVariables); for(i=0; i<_minIsSet.size(); i++) _minIsSet[i].resize(_numTimesteps);
    _min.clear();
    _min.resize(_numVariables); for(i=0; i<_min.size(); i++) _min[i].resize(_numTimesteps);
    _maxIsSet.clear();
    _maxIsSet.resize(_numVariables); for(i=0; i<_maxIsSet.size(); i++) _maxIsSet[i].resize(_numTimesteps);
    _max.clear();
    _max.resize(_numVariables); for(i=0; i<_max.size(); i++) _max[i].resize(_numTimesteps);

    /* the min and max values are in the header */
    _min[0][0] = header.amin;
    _max[0][0] = header.amax;
    _minIsSet[0][0] = true;
    _maxIsSet[0][0] = true;

    fclose(input);
  }

  void readMRC(Volume& vol,
	       const std::string& filename, 
	       unsigned int var, unsigned int time,
	       uint64 off_x, uint64 off_y, uint64 off_z,
	       const Dimension& subvoldim)
  {
    char buf[256];

    FILE *input;
    size_t i,j,k;
    bool mustSwap = false;
    MrcHeader header;
    ExtendedMrcHeader extheader;

    memset(buf,0,256);

    if(var > 0)
      throw IndexOutOfBounds("Variable index out of bounds.");
    if(time > 0)
      throw IndexOutOfBounds("Timestep index out of bounds.");
    if(subvoldim.isNull())
      throw IndexOutOfBounds("Specified subvolume dimension is null.");

    VolumeFileInfo vfi(filename);

    if((off_x + subvoldim[0] - 1 >= vfi.XDim()) ||
       (off_y + subvoldim[1] - 1 >= vfi.YDim()) ||
       (off_z + subvoldim[2] - 1 >= vfi.ZDim()))
      {
	throw IndexOutOfBounds("Subvolume specified is outside volume dimensions");
      }

    /** errors checked, now we can start modifying the volume object */
    vol.voxelType(vfi.voxelType());
    vol.dimension(subvoldim);
    vol.boundingBox(BoundingBox(vfi.XMin()+off_x*vfi.XSpan(),
				vfi.YMin()+off_y*vfi.YSpan(),
				vfi.ZMin()+off_z*vfi.ZSpan(),
				vfi.XMin()+(off_x+subvoldim[0]-1)*vfi.XSpan(),
				vfi.YMin()+(off_y+subvoldim[1]-1)*vfi.YSpan(),
				vfi.ZMin()+(off_z+subvoldim[2]-1)*vfi.ZSpan()));
    vol.min(vfi.min());
    vol.max(vfi.max());

    /*
      read the volume data
    */
    if((input = fopen(filename.c_str(),"rb")) == NULL)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error opening file '" + filename + "': " + buf;
	throw ReadError(errStr);
      }

#if 0
    /* determine if we must swap values or not */
    struct stat s;
    if(stat(filename.c_str(), &s)==-1)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error stat-ing file '" + filename + "': " + buf;
	fclose(input);
	throw ReadError(errStr);
      }

    if(fread(&header, sizeof(MrcHeader), 1, input) != 1)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error reading file '" + filename + "': " + buf;
        fclose(input);
        throw ReadError(errStr);
      }
    if(checkHeader(header, s.st_size))
      mustSwap = false;
    else 
      {
	// swap and try again
	swapHeader(header);
	if(checkHeader(header, s.st_size))
	  mustSwap = true;
	else 
	  {
	    // we dont support wierd or exotic endianness
	    fclose(input);
	    throw InvalidMRCFile("Cannot determine endianness");
	  }
      }

    if(header.nsymbt) //if there is extended header info
      {
	if(fread(&extheader, sizeof(ExtendedMrcHeader), 1, input) != 1)
	  {
	    geterrstr(errno,buf,256);
	    std::string errStr = "Error reading file '" + filename + "': " + buf;
	    fclose(input);
	    throw ReadError(errStr);
	  }

	//no need to swap the extended header because we're ignoring it...
      }
#endif

    off_t file_offx, file_offy, file_offz;
    for(k=off_z; k<=(off_z+subvoldim[2]-1); k++)
      {
	file_offz = 1024+k*vol.XDim()*vol.YDim()*vol.voxelSize();
	for(j=off_y; j<=(off_y+subvoldim[1]-1); j++)
	  {
	    file_offy = j*vol.XDim()*vol.voxelSize();
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

    if(mustSwap)
      {
	size_t len = vol.XDim()*vol.YDim()*vol.ZDim();
	switch(vol.voxelType())
	  {
	  case UShort: for(i=0;i<len;i++) SWAP_16(*vol+i*vol.voxelSize()); break;
	  case Float:  for(i=0;i<len;i++) SWAP_32(*vol+i*vol.voxelSize()); break;
	    //	  case Double: for(i=0;i<len;i++) SWAP_64(*vol+i*vol.voxelSize()); break;
	  default: break; /* no swapping needed for unsigned char data, and unsigned int is not defined for mrc */
	  }
      }
    
    //convert signed values to unsigned since volmagick doesnt support signed
    switch(vol.voxelType())
      {
      case UChar:
	{
	  vol.min(((vol.min() - SCHAR_MIN)/(SCHAR_MAX - SCHAR_MIN))*UCHAR_MAX);
	  vol.max(((vol.max() - SCHAR_MIN)/(SCHAR_MAX - SCHAR_MIN))*UCHAR_MAX);
	  size_t len = vol.XDim()*vol.YDim()*vol.ZDim();
	  for(i=0; i<len; i++)
	    {
	      char c = *((char*)(*vol+i*vol.voxelSize()));
	      *((unsigned char*)(*vol+i*vol.voxelSize())) = ((float(c) - SCHAR_MIN)/(SCHAR_MAX - SCHAR_MIN))*UCHAR_MAX;
	    }
	}
	break;
      case UShort:
	{
	  vol.min(((vol.min() - SHRT_MIN)/(SHRT_MAX - SHRT_MIN))*USHRT_MAX);
	  vol.max(((vol.max() - SHRT_MIN)/(SHRT_MAX - SHRT_MIN))*USHRT_MAX);	  
	  size_t len = vol.XDim()*vol.YDim()*vol.ZDim();
	  for(i=0; i<len; i++)
	    {
	      short c = *((short*)(*vol+i*vol.voxelSize()));
	      *((unsigned short*)(*vol+i*vol.voxelSize())) = ((float(c) - SHRT_MIN)/(SHRT_MAX - SHRT_MIN))*USHRT_MAX;
	    }
	}
	break;
      default: break;
      }

    fclose(input);
  }

  void createMRC(const std::string& filename,
		 const BoundingBox& boundingBox,
		 const Dimension& dimension,
		 const std::vector<VoxelType>& voxelTypes,
		 unsigned int numVariables, unsigned int numTimesteps,
		 double min_time, double max_time)
  {
    using namespace boost;
    MrcHeader mrcHeader;
     
    FILE *output;
    size_t i,j,k;

    if(boundingBox.isNull())
      throw InvalidBoundingBox("Bounding box must not be null");
    if(dimension.isNull())
      throw InvalidBoundingBox("Dimension must not be null");
    if(numVariables > 1)
      throw InvalidMRCHeader(str(format("MRC format only supports 1 variable (%1% requested)") % numVariables));
    if(numTimesteps > 1)
      throw InvalidMRCHeader(str(format("MRC format only supports 1 timestep (%1% requested)") % numTimesteps));
    if(voxelTypes.size() > 1)
      throw InvalidMRCHeader("MRC format only supports 1 variable and 1 timestep. (too many voxel types specified)");
    if(min_time != max_time)
      throw InvalidMRCHeader("MRC format does not support multiple timesteps. (min time and max time must be the same)");

    //check for unsupported type
    if(voxelTypes[0] == UInt ||
       voxelTypes[0] == Double ||
       voxelTypes[0] == UInt64)
      throw InvalidMRCHeader(str(format("Unsupported type: %1%") % VoxelTypeStrings[voxelTypes[0]]));

    memset(&mrcHeader,0,sizeof(MrcHeader));

    int type_conv[] = { 0, 1, 2, 2, 2, 2 };
    int type_sizes[] = { 1, 2, 4, 4, 4, 4 };
    mrcHeader.nx = dimension[0];
    mrcHeader.ny = dimension[1];
    mrcHeader.nz = dimension[2];
    mrcHeader.mode = type_conv[int(voxelTypes[0])];
    strcpy(mrcHeader.map,"MAP");
    mrcHeader.xorigin = boundingBox.minx;
    mrcHeader.yorigin = boundingBox.miny;
    mrcHeader.zorigin = boundingBox.minz;
    mrcHeader.xlength = boundingBox.maxx;
    mrcHeader.ylength = boundingBox.maxy;
    mrcHeader.zlength = boundingBox.maxz;
    if(!big_endian())
      swapHeader(mrcHeader);
    
    if((output = fopen(filename.c_str(),"wb")) == NULL)
      {
	char buf[256] = { 0 };
	geterrstr(errno,buf,256);
	std::string errStr = "Error opening file '" + filename + "': " + buf;
	throw WriteError(errStr);
      }

    if(fwrite(&mrcHeader,sizeof(mrcHeader),1,output) != 1)
      {
	char buf[256] = { 0 };
	geterrstr(errno,buf,256);
	std::string errStr = "Error writing header to file '" + filename + "': " + buf;
	fclose(output);
	throw WriteError(errStr);
      }

    scoped_array<unsigned char> scanline;
    try
      {
	scanline.reset(new unsigned char[dimension[0]*VoxelTypeSizes[voxelTypes[0]]]);
      }
    catch(std::bad_alloc& e)
      {
	fclose(output);
	throw MemoryAllocationError("Unable to allocate memory for write buffer");
      }
    memset(scanline.get(),0,dimension[0]*VoxelTypeSizes[voxelTypes[0]]*sizeof(unsigned char));
    // write a scanline at a time
    for(k=0; k<dimension[2]; k++)
      for(j=0; j<dimension[1]; j++)
	{
	  if(fwrite(scanline.get(),VoxelTypeSizes[voxelTypes[0]],dimension[0],output) != dimension[0])
	    {
	      char buf[256] = { 0 };
	      geterrstr(errno,buf,256);
	      std::string errStr = "Error writing volume data to file '" + filename + "': " + buf;
	      fclose(output);
	      throw WriteError(errStr);
	    }
	}

    fclose(output);
  }

  void writeMRC(const Volume& wvol, 
		const std::string& filename,
		unsigned int var, unsigned int time,
		uint64 off_x, uint64 off_y, uint64 off_z)
  {
    using namespace boost;
    VolumeFileInfo volinfo;
    char buf[256];
    MrcHeader header;
    ExtendedMrcHeader extheader;
    bool creatingNewFile = false;
    bool mustSwap = false;
     
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
    catch(ReadError e)
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
	createMRC(filename,box,dim,std::vector<VoxelType>(1,vol.voxelType()),1,1,0.0,0.0);
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
	
    /* determine if we must swap values or not */
    struct stat s;
    if(stat(filename.c_str(), &s)==-1)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error stat-ing file '" + filename + "': " + buf;
	fclose(output);
	throw ReadError(errStr);
      }

    if(fread(&header, sizeof(MrcHeader), 1, output) != 1)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error reading file '" + filename + "': " + buf;
        fclose(output);
        throw ReadError(errStr);
      }
    if(checkHeader(header, s.st_size))
      mustSwap = false;
    else 
      {
	// swap and try again
	swapHeader(header);
	if(checkHeader(header, s.st_size))
	  mustSwap = true;
	else 
	  {
	    // we dont support wierd or exotic endianness
	    fclose(output);
	    throw InvalidMRCFile("Cannot determine endianness");
	  }
      }

    if(header.nsymbt) //if there is extended header info
      {
	if(fread(&extheader, sizeof(ExtendedMrcHeader), 1, output) != 1)
	  {
	    geterrstr(errno,buf,256);
	    std::string errStr = "Error reading file '" + filename + "': " + buf;
	    fclose(output);
	    throw ReadError(errStr);
	  }

	//no need to swap the extended header because we're ignoring it...
      }

#warning TODO: correctly calculate the new min/max values if type is UChar or UShort
    //set the header's min/max values
    if(creatingNewFile)
      {
	header.amin = MIN(0.0,vol.min());
	header.amax = MAX(0.0,vol.max());
      }
    else
      {
	header.amin = MIN(volinfo.min(),vol.min());
	header.amax = MAX(volinfo.max(),vol.max());
      }

    outvol_xdim = header.nx;
    outvol_ydim = header.ny;
    outvol_zdim = header.nz;

    if(!big_endian())
      swapHeader(header); //always write big endian data so it's similar to rawiv
    
    if(FSEEK(output,0,SEEK_SET) == -1)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error seeking in file '" + filename + "': " + buf;
	fclose(output);
	throw ReadError(errStr);
      }

    if(fwrite(&header,sizeof(header),1,output) != 1)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error writing header to file '" + filename + "': " + buf;
	fclose(output);
	throw WriteError(errStr);
      }

    scoped_array<unsigned char> scanline;
    try
      {
	scanline.reset(new unsigned char[vol.XDim()*vol.voxelSize()]);
      }
    catch(std::bad_alloc& e)
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

	    memcpy(scanline.get(),*vol+
		   ((k-off_z)*vol.XDim()*vol.YDim()*vol.voxelSize())+
		   ((j-off_y)*vol.XDim()*vol.voxelSize()),
		   vol.XDim()*vol.voxelSize());

	    switch(vol.voxelType())
	      {
	      case UChar: break;
	      case UShort: break;
	      default: break;
	      }

	    /* swap the volume data if on little endian machine */
	    if(!big_endian())
	      {
		size_t len = vol.XDim();
		switch(vol.voxelType())
		  {
		  case UShort: for(i=0;i<len;i++) SWAP_16(scanline.get()+i*vol.voxelSize()); break;
		  case Float:  for(i=0;i<len;i++) SWAP_32(scanline.get()+i*vol.voxelSize()); break;
		  case Double: for(i=0;i<len;i++) SWAP_64(scanline.get()+i*vol.voxelSize()); break;
		  default: break; /* no swapping needed for unsigned char data, and unsigned int is not defined for rawiv */
		  }
	      }

	    if(fwrite(scanline.get(),vol.voxelSize(),vol.XDim(),output) != vol.XDim())
	      {
		geterrstr(errno,buf,256);
		std::string errStr = "Error writing volume data to file '" + filename + "': " + buf;
		fclose(output);
		throw WriteError(errStr);
	      }
	  }
      }

    fclose(output);
  }
};
