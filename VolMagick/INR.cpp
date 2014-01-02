/* $Id: INR.cpp,v 1.3 2008/08/06 16:49:32 transfix Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include <string>

#define PHOENIX_LIMIT 10
#define BOOST_SPIRIT_SELECT_LIMIT 10

#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/spirit.hpp>
#include <boost/spirit/dynamic/select.hpp>

#if 0
#include <boost/spirit/core.hpp>
#include <boost/spirit/utility.hpp>
#include <boost/spirit/symbols.hpp>
#include <boost/spirit/error_handling/exceptions.hpp>
#endif

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

using namespace boost::spirit;

/*
Inrimage format

This is the format of the FOVEA 3D images that was developed at INRIA. It is composed of a 256 characters header followed by the image raw data.

The header contains the following information:
¥ the image dimensions
¥ number of values per pixel (or voxel)
¥ the type of coding (integer, float)
¥ the size of the coding in bits
¥ the type of machine that coded the information: Sun, Dec, Intel.

An example of an image file:

      #INRIMAGE-4#{

      XDIM=128 // x dimension

      YDIM=128 // y dimension

      ZDIM=128 // z dimension

      VDIM=1 // number of scalar per voxel (1 = scalar image, 3 = 3D image of vectors)

      VX=0.66 // voxel size in x

      VY=0.66 // voxel size in y

      VZ=1 // voxel size in z

      TYPE=unsigned fixed // float, signed fixed, or unsigned fixed

      PIXSIZE=16 bits // 8, 16, 32, or 64

      SCALE=2**0 // not used

      CPU=decm // decm, alpha, pc, sun, sgi

      **}

      // little endianness : decm, alpha, pc; big endianness :sun, sgi

      // fill with carriage return or with any other information

      // until the total size of the header is 256 characters (including final newline character)

      // raw data, size=XDIM*YDIM*ZDIM*VDIM*PIXSIZE/8 
*/

namespace VolMagick
{
  enum Endianness_t { LittleEndian, BigEndian };

  VOLMAGICK_DEF_EXCEPTION(InvalidINRHeader);

  namespace INR
  {
    //semantic actions
    struct set_xdim_t
    {
      set_xdim_t(VolumeFileInfo& _vol) : vol(_vol) {}
      void operator()(unsigned int val) const { vol.dimension(Dimension(val,vol.YDim(),vol.ZDim())); }
      VolumeFileInfo& vol;
    };

    struct set_ydim_t
    {
      set_ydim_t(VolumeFileInfo& _vol) : vol(_vol) {}
      void operator()(unsigned int val) const { vol.dimension(Dimension(vol.XDim(),val,vol.ZDim())); }
      VolumeFileInfo& vol;
    };

    struct set_zdim_t
    {
      set_zdim_t(VolumeFileInfo& _vol) : vol(_vol) {}
      void operator()(unsigned int val) const { vol.dimension(Dimension(vol.XDim(),vol.YDim(),val)); }
      VolumeFileInfo& vol;
    };

    struct set_vdim_t
    {
      set_vdim_t(VolumeFileInfo& _vol) : vol(_vol) {}
      void operator()(unsigned int val) const { vol.numVariables(val); }
      VolumeFileInfo& vol;
    };

    struct set_vx_t
    {
      set_vx_t(VolumeFileInfo& _vol) : vol(_vol) {}
      void operator()(double val) const { vol.boundingBox(BoundingBox(0.0,0.0,0.0,
								      val*(vol.XDim()-1),vol.YMax(),vol.ZMax())); }
      VolumeFileInfo& vol;
    };

    struct set_vy_t
    {
      set_vy_t(VolumeFileInfo& _vol) : vol(_vol) {}
      void operator()(double val) const { vol.boundingBox(BoundingBox(0.0,0.0,0.0,
								      vol.XMax(),val*(vol.YDim()-1),vol.ZMax())); }
      VolumeFileInfo& vol;
    };

    struct set_vz_t
    {
      set_vz_t(VolumeFileInfo& _vol) : vol(_vol) {}
      void operator()(double val) const { vol.boundingBox(BoundingBox(0.0,0.0,0.0,
								      vol.XMax(),vol.YMax(),val*(vol.ZDim()-1))); }
      VolumeFileInfo& vol;
    };

    struct set_type_t
    {
      set_type_t(VolumeFileInfo& _vol, std::string& _dt, std::string& _ps) : 
	vol(_vol), datatype_str(_dt), pixsize_str(_ps) {}
      void operator()(const char *beg, const char *end) const
      {
	datatype_str = std::string(beg, end);

	if(datatype_str == "unsigned fixed" && pixsize_str == "8 bits") // UChar
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = UChar;
	else if(datatype_str == "unsigned fixed" && pixsize_str == "16 bits") // UShort
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = UShort;
	else if(datatype_str == "unsigned fixed" && pixsize_str == "32 bits") // UInt
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = UInt;
	else if(datatype_str == "float" && pixsize_str == "32 bits") // Float
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = Float;
	else if(datatype_str == "float" && pixsize_str == "64 bits") // Double
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = Double;
	else if(datatype_str == "unsigned fixed" && pixsize_str == "64 bits") // UInt64
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = UInt64;
	else throw_(beg,InvalidINRHeader("Unsupported type/pixsize combination: "
					 "No 'signed fixed' types, or pixsizes < 32 bits for 'float' types"));
      }

      VolumeFileInfo& vol;
      std::string& datatype_str;
      std::string& pixsize_str;
    };

    struct set_pixsize_t
    {
      set_pixsize_t(VolumeFileInfo& _vol, std::string& _dt, std::string& _ps) : 
	vol(_vol), datatype_str(_dt), pixsize_str(_ps) {}
      void operator()(const char *beg, const char *end) const
      {
	pixsize_str = std::string(beg, end);

	if(datatype_str == "unsigned fixed" && pixsize_str == "8 bits") // UChar
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = UChar;
	else if(datatype_str == "unsigned fixed" && pixsize_str == "16 bits") // UShort
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = UShort;
	else if(datatype_str == "unsigned fixed" && pixsize_str == "32 bits") // UInt
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = UInt;
	else if(datatype_str == "float" && pixsize_str == "32 bits") // Float
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = Float;
	else if(datatype_str == "float" && pixsize_str == "64 bits") // Double
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = Double;
	else if(datatype_str == "unsigned fixed" && pixsize_str == "64 bits") // UInt64
	  for(unsigned int i=0; i<vol.numVariables(); i++)
	    vol.voxelTypes(i) = UInt64;
	else throw_(beg,InvalidINRHeader("Unsupported type/pixsize combination: "
					 "No 'signed fixed' types, or pixsizes < 32 bits for 'float' types"));
      }

      VolumeFileInfo& vol;
      std::string& datatype_str;
      std::string& pixsize_str;
    };

    struct set_cpu_t
    {
      set_cpu_t(Endianness_t& _e) : endianness(_e) {}
      void operator()(const char *beg, const char *end) const
      {
	std::string cpu_str(beg, end);
	if(cpu_str == "decm" || cpu_str == "alpha" || cpu_str == "pc")
	  endianness = LittleEndian;
	else if(cpu_str == "sun" || cpu_str == "sgi")
	  endianness = BigEndian;
	else throw_(beg,InvalidINRHeader("Invalid CPU type.  Choose one of decm, alpha, pc, sun, or sgi"));
      }
      Endianness_t& endianness;
    };

    //comment rule
    struct skip_grammar : public grammar<skip_grammar>
    {
      template <typename ScannerT>
      struct definition
      {
	definition(skip_grammar const& /*self*/)
	{
	  skip 
	    = space_p
	    | comment_p("//")
	    | comment_p("#")
	    ;
	}

	rule<ScannerT> skip;
	rule<ScannerT> const& start() const { return skip; }
      };
    };

    struct handler
    {
      template <typename ScannerT, typename ErrorT>
      error_status<>
      operator()(ScannerT const& /*scan*/, ErrorT const& /*error*/) const
      {
	return error_status<>(error_status<>::fail/*,-1,error*/);
      }
    };

    //main grammar
    struct inr_grammar : public grammar<inr_grammar>
    {
      inr_grammar(VolumeFileInfo& _vol) : 
	vol(_vol), set_xdim(_vol), set_ydim(_vol), set_zdim(_vol),
	set_vdim(_vol), set_vx(_vol), set_vy(_vol), set_vz(_vol),
	set_type(_vol,datatype_str,pixsize_str),
	set_pixsize(_vol,datatype_str,pixsize_str), set_cpu(endianness),
        datatype_str("unsigned fixed"), pixsize_str("8 bits") {}

      VolumeFileInfo& vol;
      set_xdim_t set_xdim;
      set_ydim_t set_ydim;
      set_zdim_t set_zdim;
      set_vdim_t set_vdim;
      set_vx_t set_vx;
      set_vy_t set_vy;
      set_vz_t set_vz;
      set_type_t set_type;
      set_pixsize_t set_pixsize;
      set_cpu_t set_cpu;
      
      std::string datatype_str;
      std::string pixsize_str;
      Endianness_t endianness; //endianness of the data, not the runtime machine

      template <typename ScannerT>
      struct definition
      {
	definition(inr_grammar const& self)
	{
	  XDIM = strlit<>("XDIM");
	  YDIM = strlit<>("YDIM");
	  ZDIM = strlit<>("ZDIM");
	  VDIM = strlit<>("VDIM");
	  VX = strlit<>("VX");
	  VY = strlit<>("VY");
	  VZ = strlit<>("VZ");
	  DATATYPE = strlit<>("TYPE");
	  PIXSIZE = strlit<>("PIXSIZE");
	  SCALE = strlit<>("SCALE");
	  CPU = strlit<>("CPU");

#if 0
	  xdim = as_lower_d["xdim"] >> ch_p('=') >> uint_p[self.set_xdim];
	  ydim = as_lower_d["ydim"] >> ch_p('=') >> uint_p[self.set_ydim];
	  zdim = as_lower_d["zdim"] >> ch_p('=') >> uint_p[self.set_zdim];
	  vdim = as_lower_d["vdim"] >> ch_p('=') >> uint_p[self.set_vdim];
	  vx = as_lower_d["vx"] >> ch_p('=') >> ureal_p[self.set_vx];
	  vy = as_lower_d["vy"] >> ch_p('=') >> ureal_p[self.set_vy];
	  vz = as_lower_d["vz"] >> ch_p('=') >> ureal_p[self.set_vz];
	  datatype = as_lower_d["type"] >> ch_p('=') >> lexeme_d[ as_lower_d["float"] 
								  | as_lower_d["signed fixed"] 
								  | as_lower_d["unsigned fixed"] ][self.set_type];
	  pixsize = as_lower_d["pixsize"] >> ch_p('=') >> lexeme_d[ as_lower_d["8 bits"]
								    | as_lower_d["16 bits"]
								    | as_lower_d["32 bits"]
								    | as_lower_d["64 bits"] ][self.set_pixsize];
	  scale = as_lower_d["scale"] >> ch_p('=') >> *anychar_p; //ignore this
	  cpu = as_lower_d["cpu"] >> ch_p('=') >> (*anychar_p)[self.set_cpu];
#endif 

	  var = XDIM
	    | YDIM
	    | ZDIM
	    | VDIM
	    | VX
	    | VY
	    | VZ
	    | DATATYPE
	    | PIXSIZE
	    | SCALE
	    | CPU
	    ;

	  /*value = uint_p 
	    | ureal_p
	    | (as_lower_d["float"]
	       | ((as_lower_d["signed"] | as_lower_d["unsigned"]) >> as_lower_d["fixed"]))
	    | (uint_p >> as_lower_d["bits"])
	    | *anychar_p //catch all for SCALE
	    ;*/

	  value = select_p(uint_p,
			   ureal_p,
			   (as_lower_d["float"]
			    | ((as_lower_d["signed"] | as_lower_d["unsigned"]) >> as_lower_d["fixed"])),
			   (uint_p >> as_lower_d["bits"]),
			   *anychar_p);
	    
	  statement = var >> '=' >> value;

	  translation_unit = my_guard( *statement )[handler()];
	  //parser start symbol
	  //translation_unit = my_guard( xdim | ydim /* *(xdim | ydim | zdim | vdim | datatype |
	  //				  pixsize | cpu | vx | vy | vz | scale)*/ )[handler()];
	}

	//defined variables
	rule<ScannerT> XDIM, YDIM, ZDIM, VDIM, VX, VY, VZ, DATATYPE, PIXSIZE, SCALE, CPU;

	rule<ScannerT>
          statement, var, value;

#if 0
	rule<ScannerT> 
	  xdim, ydim, zdim, vdim, datatype, 
	  pixsize, cpu, vx, vy, vz, scale;
#endif

	guard<InvalidINRHeader> my_guard;

	rule<ScannerT> translation_unit;
	rule<ScannerT> const& start() const 
	{ 
	  //return my_guard(translation_unit);
	  //return my_guard(translation_unit)[handler()];
	  return translation_unit;
	}
      };
    };
  };

  void VolumeFileInfo::readINR(const std::string& file)
  {
    char header_c[256], buf[256];
    FILE *input;
    INR::inr_grammar g(*this);
    INR::skip_grammar skip;

    if((input = fopen(file.c_str(),"rb")) == NULL)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error opening file '" + file + "': " + buf;
	throw ReadError(errStr);
      }
    
    if(fread(header_c, 256, 1, input) != 1)
      {
	geterrstr(errno,buf,256);
        std::string errStr = "Error reading file '" + file + "': " + buf;
        fclose(input);
        throw ReadError(errStr);
      }

    std::string header(header_c,header_c+256);
    parse_info<> info = parse(header.c_str(), g, skip);

    if(!info.full)
      {
	std::cout << "Failed at this point: " << info.stop << std::endl;
	throw InvalidINRHeader("Failed to parse header! Check that it is formatted correctly.");
      }

    /*
    std::string::const_iterator start = header.begin();
    std::string::const_iterator end = header.end();
    
    parse_info<std::string::const_iterator> result = 
      parse(start, end, g, skip);
    
    if(!result.full) //throw if parse fails
      {
	std::string stopstring(result.stop,end);
	stopstring.resize(20);
	throw InvalidINRHeader("Failed to parse header starting here: " + stopstring);
      }
    else
      {
	
      }*/
  }

  void readINR(Volume& vol,
	       const std::string& filename, 
	       unsigned int var, unsigned int time,
	       uint64 off_x, uint64 off_y, uint64 off_z,
	       const Dimension& subvoldim)
  {
    VolumeFileInfo volinfo;
    INR::inr_grammar g(volinfo);
  }

  void createINR(const std::string& filename,
		 const BoundingBox& boundingBox,
		 const Dimension& dimension,
		 const std::vector<VoxelType>& voxelTypes,
		 unsigned int numVariables, unsigned int numTimesteps,
		 double min_time, double max_time)
  {
    char buf[256];
    char header[256];

    FILE *output;
    size_t i,j,k,v;

    memset(buf,0,256);
    memset(header,'\n',256);

    if(boundingBox.isNull())
      throw InvalidBoundingBox("Bounding box must not be null");
    if(dimension.isNull())
      throw InvalidBoundingBox("Dimension must not be null");
    if(numTimesteps > 1)
      throw InvalidINRHeader("INR format only supports 1 timestep.");
    if(min_time != max_time)
      throw InvalidINRHeader("INR format does not support multiple timesteps.");
    if(numVariables == 0)
      throw InvalidINRHeader("Must have at least 1 variable");
    if(numVariables != voxelTypes.size())
      throw InvalidINRHeader("voxelTypes vector must specify a type for every variable.");

    for(std::vector<VoxelType>::const_iterator i = voxelTypes.begin();
	i != voxelTypes.end();
	i++)
      if(*i != *(voxelTypes.begin()))
	throw InvalidINRHeader("INR format requires all variables to have the same datatype.");

    std::string datatype, pixsize;
    switch(voxelTypes[0])
      {
      default:
      case UChar: datatype = "unsigned fixed"; pixsize = "8 bits"; break;
      case UShort: datatype = "unsigned fixed"; pixsize = "16 bits"; break;
      case UInt: datatype = "unsigned fixed"; pixsize = "32 bits"; break;
      case Float: datatype = "float"; pixsize = "32 bits"; break;
      case Double: datatype = "float"; pixsize = "64 bits"; break;
      case UInt64: datatype = "unsigned fixed"; pixsize = "64 bits"; break;
      }

    strncpy(header,boost::str(boost::format("#INRIMAGE-4#{\n"
					    "XDIM=%1%\n"
					    "YDIM=%2%\n"
					    "ZDIM=%3%\n"
					    "VDIM=%4%\n"
					    "VX=%5%\n"
					    "VY=%6%\n"
					    "VZ=%7%\n"
					    "TYPE=%8%\n"
					    "PIXSIZE=%9%\n"
					    "CPU=%10%\n")
			      % dimension[0]
			      % dimension[1]
			      % dimension[2]
			      % numVariables
			      % ((boundingBox.maxx - boundingBox.minx)/(dimension[0]-1))
			      % ((boundingBox.maxy - boundingBox.miny)/(dimension[1]-1))
			      % ((boundingBox.maxz - boundingBox.minz)/(dimension[2]-1))
			      % datatype
			      % pixsize
			      % (big_endian()?"sun":"pc")).c_str(),
	    256);
    memcpy(header+256-4,"##}\n",4);
    
    if((output = fopen(filename.c_str(),"wb")) == NULL)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error opening file '" + filename + "': " + buf;
	throw WriteError(errStr);
      }

    if(fwrite(header,256,1,output) != 1)
      {
	geterrstr(errno,buf,256);
	std::string errStr = "Error writing header to file '" + filename + "': " + buf;
	fclose(output);
	throw WriteError(errStr);
      }

     //write a scanline at a time
    for(v=0; v<numVariables; v++)
      {
	// each variable may have its own type, so scanline size may differ
	unsigned char * scanline = (unsigned char *)calloc(dimension[0]*VoxelTypeSizes[voxelTypes[v]],1);
	if(scanline == NULL)
	  {
	    fclose(output);
	    throw MemoryAllocationError("Unable to allocate memory for write buffer");
	  }

	for(k=0; k<dimension[2]; k++)
	  for(j=0; j<dimension[1]; j++)
	    {
	      if(fwrite(scanline,VoxelTypeSizes[voxelTypes[v]],dimension[0],output) != dimension[0])
		{
		  geterrstr(errno,buf,256);
		  std::string errStr = "Error writing volume data to file '" + filename + "': " + buf;
		  free(scanline);
		  fclose(output);
		  throw WriteError(errStr);
		}
	    }

	free(scanline);
      }

    fclose(output);
  }

  void writeINR(const Volume& wvol, 
		const std::string& filename,
		unsigned int var, unsigned int time,
		uint64 off_x, uint64 off_y, uint64 off_z)
  {
    
  }
};
