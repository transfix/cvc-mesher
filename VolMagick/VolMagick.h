/* $Id: VolMagick.h,v 1.4 2008/08/15 21:53:04 transfix Exp $ */

#ifndef __VOLMAGICK_H__
#define __VOLMAGICK_H__

#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <cmath>
//#include <string.h>
//#include <math.h>

#include <boost/shared_array.hpp>

#ifdef min
#ifndef _MSC_VER
#warning The macro 'min' was defined.  The volmagick library uses 'min' in the VolMagick::Voxels class and must undefine it for the definition of the class to compile correctly! Sorry!
#endif
#undef min
#endif

#ifdef max
#ifndef _MSC_VER
#warning The macro 'max' was defined.  The volmagick library uses 'max' in the VolMagick::Voxels class and must undefine it for the definition of the class to compile correctly! Sorry!
#endif
#undef max
#endif

#include "Exceptions.h"
#include "Dimension.h"
#include "BoundingBox.h"

//use a safer min/max
#if 0
#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif
#endif

namespace VolMagick
{
  template <class T> const T& MIN(const T& a, const T& b) { return std::min(a,b); }
  template <class T> const T& MAX(const T& a, const T& b) { return std::max(a,b); }

#ifndef _MSC_VER
  typedef unsigned long long uint64;
  typedef long long int64;
#else
  typedef unsigned __int64 uint64;
  typedef __int64 int64;
#endif

  static inline uint64 upToPowerOfTwo(uint64 value)
    {
      uint64 c = 0;
      uint64 v = value;
      
      /* round down to nearest power of two */
      while (v>1)
	{
	  v = v>>1;
	  c++;
	}
      
      /* if that isn't exactly the original value */
      if ((v<<c)!=value)
	{
	  /* return the next power of two */
	  return (v<<(c+1));
	}
      else
	{
	  /* return this power of two */
	  return (v<<c);
	}
    }

  enum VoxelType { UChar, UShort, UInt, Float, Double, UInt64 };
  static const unsigned int VoxelTypeSizes[] = { 1, 2, 4, 4, 8, 8};
  static const char * VoxelTypeStrings[] = {"unsigned char","unsigned short","unsigned int","float","double","unsigned int64"};

  class Voxels;

  //this class is used to provide calling code a means of getting periodic control during long operations
  class VoxelOperationStatusMessenger
  {
  public:
    enum Operation { CalculatingMinMax, 
		     CalculatingMin, CalculatingMax, SubvolumeExtraction,
		     Fill, Map, Resize, Composite, BilateralFilter, 
		     ContrastEnhancement, AnisotropicDiffusion, CombineWith,
                     ReadVolumeFile, WriteVolumeFile, CreateVolumeFile, CalcGradient };
    
    VoxelOperationStatusMessenger() {}
    virtual ~VoxelOperationStatusMessenger() {}

    virtual void start(const Voxels *vox, Operation op, uint64 numSteps) const = 0;
    virtual void step(const Voxels *vox, Operation op, uint64 curStep) const = 0;
    virtual void end(const Voxels *vox, Operation op) const = 0;
  };

  //sets the default messenger for all voxel objects to use
  void setDefaultMessenger(const VoxelOperationStatusMessenger* vosm);

  class CompositeFunction;

  class Voxels
  {
  public:
    Voxels(const Dimension& d = Dimension(4,4,4), VoxelType vt = UChar);
    Voxels(const void *v, const Dimension& d, VoxelType vt);
    Voxels(const Voxels& v);
    virtual ~Voxels();

    /*
      Voxels Dimensions
    */
    Dimension& dimension() { return _dimension; }
    const Dimension& dimension() const { return _dimension; }
    virtual void dimension(const Dimension& d);
    uint64 XDim() const { return dimension().xdim; }
    uint64 YDim() const { return dimension().ydim; }
    uint64 ZDim() const { return dimension().zdim; }

    /*
      Voxel I/O
    */
    double operator()(uint64 i) const /* reading a voxel value */
    {
      if(i >= XDim()*YDim()*ZDim()) 
	throw IndexOutOfBounds("");
      
      switch(voxelType())
	{
	case UChar:
	  return double(*((unsigned char *)(_voxels.get()+i*voxelSize())));
	case UShort:
	  return double(*((unsigned short *)(_voxels.get()+i*voxelSize())));
	case UInt:
	  return double(*((unsigned int *)(_voxels.get()+i*voxelSize())));
	case Float:
	  return double(*((float *)(_voxels.get()+i*voxelSize())));
	case Double:
	  return double(*((double *)(_voxels.get()+i*voxelSize())));
	case UInt64:
	  return double(*((uint64 *)(_voxels.get()+i*voxelSize())));
	}
      return 0;
    }
    double operator()(uint64 i, uint64 j, uint64 k) const /* reading a voxel value */
    {
      return (*this)(i+j*XDim()+k*XDim()*YDim());
    }
    
    void operator()(uint64 i, double val) /* writing a voxel value */
    {
      if(i >= XDim()*YDim()*ZDim()) 
	throw IndexOutOfBounds("");

      preWriteCopy();

      switch(voxelType())
	{
	case UChar:
	  *((unsigned char *)(_voxels.get()+i*voxelSize())) = (unsigned char)(val);
	  break;
	case UShort:
	  *((unsigned short *)(_voxels.get()+i*voxelSize())) = (unsigned short)(val);
	  break;
	case UInt:
	  *((unsigned int *)(_voxels.get()+i*voxelSize())) = (unsigned int)(val);
	  break;
	case Float:
	  *((float *)(_voxels.get()+i*voxelSize())) = float(val);
	  break;
	case Double:
	  *((double *)(_voxels.get()+i*voxelSize())) = double(val);
	  break;
	case UInt64:
	  *((uint64 *)(_voxels.get()+i*voxelSize())) = uint64(val);
	}

      //NOTE: we cant modify min/max here because it would mess up a map() operation, and perhaps other things
      //if(_minIsSet && val < min()) min(val);
      //if(_maxIsSet && val > max()) max(val);
    }
    void operator()(uint64 i, uint64 j, uint64 k, double val) /* writing a voxel value */
    {
      (*this)(i+j*XDim()+k*XDim()*YDim(),val);
    }
    
    unsigned char * operator*() { return _voxels.get(); }
    const unsigned char * operator*() const { return _voxels.get(); }

    VoxelType voxelType() const { return _voxelType; }
    void voxelType(VoxelType);
    uint64 voxelSize() const { return VoxelTypeSizes[voxelType()]; }
    const char * voxelTypeStr() const { return VoxelTypeStrings[voxelType()]; }

     /* min and max values */
    double min() const { if(!_minIsSet) calcMinMax(); return _min; }
    void min(double m) { _min = m; _minIsSet = true; }
    double max() const { if(!_maxIsSet) calcMinMax(); return _max; }
    void max(double m) { _max = m; _maxIsSet = true; }
    void unsetMinMax() { _minIsSet = _maxIsSet = false; }
    bool minIsSet() const { return _minIsSet; }
    bool maxIsSet() const { return _maxIsSet; }

    /* calculate min and max values for selected subvolumes */
    double min(uint64 off_x, uint64 off_y, uint64 off_z,
	       const Dimension& dim) const;
    double max(uint64 off_x, uint64 off_y, uint64 off_z,
	       const Dimension& dim) const;

    Voxels& operator=(const Voxels& vox) { copy(vox); return *this; }

    void messenger(const VoxelOperationStatusMessenger* vosm) { _vosm = vosm; }
    const VoxelOperationStatusMessenger* messenger() const { return _vosm; }

    /*
      operations!
    */
    virtual Voxels& copy(const Voxels& vox); //turns this object into a copy of vox
    //subvolume extraction: removes voxels outside of the subvolume specified
    virtual Voxels& sub(uint64 off_x, uint64 off_y, uint64 off_z,
			const Dimension& subvoldim);
    Voxels& fill(double val); //set all voxels to the specified value
    Voxels& fillsub(uint64 off_x, uint64 off_y, uint64 off_z,
		    const Dimension& subvoldim, double val); //set all voxels in specified subvolume to val
    Voxels& map(double min_, double max_); //maps voxels from min to max
    Voxels& resize(const Dimension& newdim); //resizes this object to the specified dimension using trilinear interpolation
    Voxels& bilateralFilter(double radiometricSigma = 200.0, double spatialSigma = 1.5, unsigned int filterRadius = 2);
    //Voxels& rotate(double deg_x, double deg_y, double deg_z); //rotates the object about the x,y,z axis
    /*
      compose vox into this object using the specified composite function.  Yes, the offset may be negative.
      Only the voxels that overlap will be subject to the composition function.
    */
    virtual Voxels& composite(const Voxels& compVox, int64 off_x, int64 off_y, int64 off_z, const CompositeFunction& func);
    
    /*
     * Contrast enhancement: enhances contrast between voxel values. 'resistor' must be a value between 0.0 and 1.0
     * Requres memory to hold the original volume + 6x the original volume using float values for voxels...
     */
    virtual Voxels& contrastEnhancement(double resistor = 0.95);

    virtual Voxels& anisotropicDiffusion(unsigned int iterations = 20);

    //special access to the shared array - careful with this!
    // 01/11/2014 - Joe R. - creation
    const boost::shared_array<unsigned char>& data() const { return _voxels; }
    boost::shared_array<unsigned char>& data() { return _voxels; }

  protected:
    void calcMinMax() const;
    void preWriteCopy()
    {
      if(_voxels.unique()) return; //nothing to copy if our voxels are already unique

      try
	{
	  boost::shared_array<unsigned char> tmp(_voxels);
	  _voxels.reset(new unsigned char[XDim()*YDim()*ZDim()*voxelSize()]);
	  memcpy(_voxels.get(),tmp.get(),XDim()*YDim()*ZDim()*voxelSize());
	}
      catch(std::bad_alloc& e)
	{
	  throw MemoryAllocationError("Could not allocate memory for voxels during copy-on-write!");
	}
    }

    //unsigned char *_voxels;
    boost::shared_array<unsigned char> _voxels;

    Dimension _dimension;

    VoxelType _voxelType;

    mutable bool _minIsSet;
    mutable double _min;
    mutable bool _maxIsSet;
    mutable double _max;

    const VoxelOperationStatusMessenger* _vosm;
  };

  class CompositeFunction
  {
  public:
    CompositeFunction() {}
    virtual ~CompositeFunction() {}

    /*
      in_vox - input voxels
      in_i,j,k - input indices specifying the current input voxel being composited
      this_vox - the destination voxels object where the result of the composition will be stored
      this_i,j,k - the destination voxel indices
      returns - the result of the composition
    */
    virtual double operator()(const Voxels& in_vox, uint64 in_i, uint64 in_j, uint64 in_k,
			      const Voxels& this_vox,uint64 this_i, uint64 this_j, uint64 this_k) const = 0;
  };

  /*
    Replaces the destination voxel with the input voxel
  */
  class CopyFunc : public CompositeFunction
  {
  public:
    CopyFunc() {}
    virtual ~CopyFunc() {}

    virtual double operator()(const Voxels& in_vox, uint64 in_i, uint64 in_j, uint64 in_k,
			      const Voxels& this_vox,uint64 this_i, uint64 this_j, uint64 this_k) const
    {
      return in_vox(in_i,in_j,in_k); /* make compiler happy.... */ this_vox(0); this_i=0; this_j=0; this_k=0;
    }
  };

  /*
    Adds the input voxel to the destination voxel;
  */
  class AddFunc : public CompositeFunction
  {
  public:
    AddFunc() {}
    virtual ~AddFunc() {}

    virtual double operator()(const Voxels& in_vox, uint64 in_i, uint64 in_j, uint64 in_k,
			      const Voxels& this_vox,uint64 this_i, uint64 this_j, uint64 this_k) const
    {
      return this_vox(this_i,this_j,this_k) + in_vox(in_i,in_j,in_k);
    }
  };

  /*
    Subtracts the destination voxel with the input voxel
  */
  class SubtractFunc : public CompositeFunction
  {
  public:
    SubtractFunc() {}
    virtual ~SubtractFunc() {}

    virtual double operator()(const Voxels& in_vox, uint64 in_i, uint64 in_j, uint64 in_k,
			      const Voxels& this_vox,uint64 this_i, uint64 this_j, uint64 this_k) const
    {
      return this_vox(this_i,this_j,this_k) - in_vox(in_i,in_j,in_k);
    }
  };

  class Volume : public Voxels
  {
  public:
    Volume(const Dimension& d = Dimension(4,4,4), 
	   VoxelType vt = UChar, 
	   const BoundingBox& box = BoundingBox(-0.5,-0.5,-0.5,0.5,0.5,0.5)) 
      : Voxels(d,vt), _boundingBox(box) {}
    Volume(const unsigned char *v, 
	   const Dimension& d, 
	   VoxelType vt, 
	   const BoundingBox& box = BoundingBox(-0.5,-0.5,-0.5,0.5,0.5,0.5))
      : Voxels(v,d,vt), _boundingBox(box) {}
    Volume(const Voxels& vox,
	   const BoundingBox& box = BoundingBox(-0.5,-0.5,-0.5,0.5,0.5,0.5))
      : Voxels(vox), _boundingBox(box) {}
    Volume(const Volume& vol)
      : Voxels(vol), _boundingBox(vol.boundingBox()), _desc(vol.desc()) {}
    ~Volume() {}

    /*
      Bounding box in object space
     */
    BoundingBox& boundingBox() { return _boundingBox; }
    const BoundingBox& boundingBox() const { return _boundingBox; }
    void boundingBox(const BoundingBox& box) { _boundingBox = box; }

    double XMin() const { return boundingBox().minx; }
    double XMax() const { return boundingBox().maxx; }
    double YMin() const { return boundingBox().miny; }
    double YMax() const { return boundingBox().maxy; }
    double ZMin() const { return boundingBox().minz; }
    double ZMax() const { return boundingBox().maxz; }

    double XSpan() const { return XDim()-1 == 0 ? 0.0 : (boundingBox().maxx-boundingBox().minx)/(XDim()-1); }
    double YSpan() const { return YDim()-1 == 0 ? 0.0 : (boundingBox().maxy-boundingBox().miny)/(YDim()-1); }
    double ZSpan() const { return ZDim()-1 == 0 ? 0.0 : (boundingBox().maxz-boundingBox().minz)/(ZDim()-1); }

    /*
      Volume description (used when this object is being saved and the volume
      format supports volume descriptions)
    */
    const std::string& desc() const { return _desc; }
    void desc(const std::string& d) { _desc = d; }

    Volume& operator=(const Volume& vol) { copy(vol); return *this; }

    /*
      Operations!
    */
    virtual Volume& copy(const Volume& vol); // makes this a copy of vol
    virtual Volume& sub(uint64 off_x, uint64 off_y, uint64 off_z,
			const Dimension& subvoldim
#ifdef _MSC_VER
			, int brain_damage = 1 //avoiding VC++ error C2555
#endif
			);
    /*
      compose volumes using object space coordinates.  Makes a duplicate of compVol and resizes it to match
      the grid resolution of this volume, then does normal voxel composition.
    */
    //virtual Volume& compositeObj(const Volume& compVol, double off_x, double off_y, double off_z, const CompositeFunction& func);

    virtual Volume& sub(const BoundingBox& subvolbox); //Gets a subvolume from a bounding box.
                                                       //Aims to keep the span of the subvolume
                                                       //as close as possible to the original.

    //Creates a subvolume with a bounding box == subvolbox, and a dimension == subvoldim
    virtual Volume& sub(const BoundingBox& subvolbox, const Dimension& subvoldim);

    //returns a linearly interpolated voxel value for the object coordinates supplied.  The coordinates must
    //be inside the bounding box, or an exception is thrown.
    double interpolate(double obj_x, double obj_y, double obj_z) const;

    //makes this volume into a new volume that contains both this volume and the volume specified, bounding box and all
    //If dimension is specified, this volume will be resized to that dimension
    Volume& combineWith(const Volume& vol, const Dimension& dim);
    Volume& combineWith(const Volume& vol);

  protected:
    BoundingBox _boundingBox;
    std::string _desc;
  };

  class VolumeFileInfo
  {
  public:
    VolumeFileInfo() : 
      _numVariables(0), _numTimesteps(0) {}
    VolumeFileInfo(const std::string& file) { VolumeFileInfo(); read(file); }
    VolumeFileInfo(const VolumeFileInfo& vfi) :
      _dimension(vfi.dimension()), _boundingBox(vfi.boundingBox()),
      _minIsSet(vfi._minIsSet), _min(vfi._min), 
      _maxIsSet(vfi._maxIsSet), _max(vfi._max),
      _numVariables(vfi.numVariables()), _numTimesteps(vfi.numTimesteps()),
      _voxelTypes(vfi._voxelTypes), _filename(vfi.filename()), _names(vfi._names),
      _tmin(0.0), _tmax(0.0)
      {}
    ~VolumeFileInfo() {}

    VolumeFileInfo& operator=(const VolumeFileInfo& vfi)
      {
	_dimension = vfi.dimension();
	_boundingBox = vfi.boundingBox();
	_minIsSet = vfi._minIsSet;
	_min = vfi._min;
	_maxIsSet = vfi._maxIsSet;
	_max = vfi._max;
	_numVariables = vfi.numVariables();
	_numTimesteps = vfi.numTimesteps();
	_voxelTypes = vfi._voxelTypes;
	_filename = vfi.filename();
	_names = vfi._names;
	_tmin = vfi.TMin();
	_tmax = vfi.TMax();
	return *this;
      }

    /*
      call VolumeFileInfo::read() to fill out this object from
      the info in the supplied file header.
    */
    void read(const std::string& file);

    /***** Volume info accessors *****/
    /*
      Volume Dimensions
    */
    Dimension& dimension() { return _dimension; }
    const Dimension& dimension() const { return _dimension; }
    void dimension(const Dimension& d) { _dimension = d; }
    uint64 XDim() const { return dimension().xdim; }
    uint64 YDim() const { return dimension().ydim; }
    uint64 ZDim() const { return dimension().zdim; }

    /*
      Bounding box in object space
     */
    BoundingBox& boundingBox() { return _boundingBox; }
    const BoundingBox& boundingBox() const { return _boundingBox; }
    void boundingBox(const BoundingBox& box) { _boundingBox = box; }
    
    double XMin() const { return boundingBox().minx; }
    double XMax() const { return boundingBox().maxx; }
    double YMin() const { return boundingBox().miny; }
    double YMax() const { return boundingBox().maxy; }
    double ZMin() const { return boundingBox().minz; }
    double ZMax() const { return boundingBox().maxz; }
    void TMin(double t) { _tmin = t; }
    double TMin() const { return _tmin; }
    void TMax(double t) { _tmax = t; }
    double TMax() const { return _tmax; }

    double XSpan() const { return XDim()-1 == 0 ? 1.0 : (boundingBox().maxx-boundingBox().minx)/(XDim()-1); }
    double YSpan() const { return YDim()-1 == 0 ? 1.0 : (boundingBox().maxy-boundingBox().miny)/(YDim()-1); }
    double ZSpan() const { return ZDim()-1 == 0 ? 1.0 : (boundingBox().maxz-boundingBox().minz)/(ZDim()-1); }
    double TSpan() const { return (TMax()-TMin())/numTimesteps(); }

    /* min and max voxel values */
    double min(unsigned int var = 0, unsigned int time = 0) const 
    { if(!_minIsSet[var][time]) calcMinMax(var,time); return _min[var][time]; }
    void min(double val, unsigned int var, unsigned int time)
    { _min[var][time] = val; _minIsSet[var][time] = true; }
    double max(unsigned int var = 0, unsigned int time = 0) const 
    { if(!_maxIsSet[var][time]) calcMinMax(var,time); return _max[var][time]; }
    void max(double val, unsigned int var, unsigned int time)
    { _max[var][time] = val; _maxIsSet[var][time] = true; }

    void numVariables(unsigned int vars) { _numVariables = vars; _voxelTypes.resize(vars); }
    unsigned int numVariables() const { return _numVariables; }
    void numTimesteps(unsigned int times) { _numTimesteps = times; }
    unsigned int numTimesteps() const { return _numTimesteps; }
    std::vector<VoxelType> voxelTypes() const { return _voxelTypes; }
    std::vector<VoxelType>& voxelTypes() { return _voxelTypes; }
    VoxelType voxelTypes(unsigned int vt) const { return _voxelTypes[vt]; }
    VoxelType& voxelTypes(unsigned int vt) { return _voxelTypes[vt]; }
    VoxelType voxelType() const { return voxelTypes(0); }
    std::string voxelTypeStr(unsigned vt = 0) const { return std::string(VoxelTypeStrings[voxelTypes(vt)]); }
    uint64 voxelSizes(unsigned int vt = 0) const { return VoxelTypeSizes[voxelTypes(vt)]; }

    std::string filename() const { return _filename; }

    std::string name(unsigned int var = 0) const { return _names[var]; }

    bool isSet() const { return !filename().empty(); }

  private:
    void calcMinMax(unsigned int var = 0, unsigned int time = 0) const;

    /*
      header i/o functions called by VolumeFileInfo::read()
    */
    void readRawIV(const std::string& file);
    void readRawV(const std::string& file);
    void readMRC(const std::string& file);
    void readINR(const std::string& file);
    void readSpider(const std::string& file);

    Dimension _dimension;
    BoundingBox _boundingBox;
    mutable std::vector<std::vector<bool> > _minIsSet;
    mutable std::vector<std::vector<double> > _min;
    mutable std::vector<std::vector<bool> > _maxIsSet;
    mutable std::vector<std::vector<double> > _max;
    unsigned int _numVariables;
    unsigned int _numTimesteps;
    std::vector<VoxelType> _voxelTypes;
    std::string _filename;
    std::vector<std::string> _names;

    double _tmin, _tmax;
  };

  /*
    read the specified subvolume from the specified file and copy it to the object
    vol.
  */
  void readVolumeFile(Volume& vol,
		      const std::string& filename, 
		      unsigned int var, unsigned int time,
		      uint64 off_x, uint64 off_y, uint64 off_z,
		      const Dimension& subvoldim);
  /*
    read the entire volume from the specified file and copy it to the object vol.
  */
  void readVolumeFile(Volume& vol, 
		      const std::string& filename,
		      unsigned int var = 0, unsigned int time = 0);

  /*
    Read multi-volume file and add each volume to the vector vols.
  */
  void readVolumeFile(std::vector<Volume>& vols,
		      const std::string& filename);

  /*
    write the specified volume to the specified offset in file 'filename'
  */
  void writeVolumeFile(const Volume& vol, 
		       const std::string& filename,
		       unsigned int var = 0, unsigned int time = 0,
		       uint64 off_x = 0, uint64 off_y = 0, uint64 off_z = 0);

  /*
    Writes the vector 'vols' to the specified file.  Make sure that the file extension
    specified is for a volume file type that supports multi-volumes. Assumes 1 timestep.
  */
  void writeVolumeFile(const std::vector<Volume>& vols,
		       const std::string& filename);

  /*
    Creates a volume file using the specified information.
  */
  void createVolumeFile(const std::string& filename,
			const BoundingBox& boundingBox,
			const Dimension& dimension,
			const std::vector<VoxelType>& voxelTypes = std::vector<VoxelType>(1, UChar),
			unsigned int numVariables = 1, unsigned int numTimesteps = 1,
			double min_time = 0.0, double max_time = 0.0);

  /*
    Calculates the gradient vector field of the input voxels, and returns the xyz vector values as 3 volumes in 'grad'
    'vt' is the voxel type of the gradient volumes.  If 'vt' is of integral type (UChar, UShort, UInt), the first
    half of the set of integers maps to [-1.0,0) and the last half maps to (0,1.0].
  */
  void calcGradient(std::vector<Volume>& grad, const Volume& vol, VoxelType vt = Float);

  /*
    Copies a subvolume of vol to dest.
  */
  void sub(Volume& dest, const Volume& vol, 
	   uint64 off_x, uint64 off_y, uint64 off_z,
	   const Dimension& subvoldim);
};

#endif

