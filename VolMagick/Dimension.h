/* $Id: Dimension.h,v 1.2 2008/02/01 20:12:12 transfix Exp $ */

#ifndef __VOLMAGICK_DIMENSION_H__
#define __VOLMAGICK_DIMENSION_H__

// If your compiler complains the "The class "VolMagick::Dimension" has no member "xdim"."
// Add your architecture Q_OS_XXXX flag (see qglobal.h) in this list.
//#if defined (Q_OS_IRIX) || defined (Q_OS_AIX) || defined (Q_OS_HPUX)
# define UNION_NOT_SUPPORTED
//#endif

namespace VolMagick
{
#ifndef _MSC_VER
  typedef unsigned long long uint64;
  typedef long long int64;
#else
  typedef unsigned __int64 uint64;
  typedef __int64 int64;
#endif

  class Dimension
  {
  public:
    /* The internal data representation is public. */
#if defined (DOXYGEN) || defined (UNION_NOT_SUPPORTED)
    uint64 xdim, ydim, zdim;
#else
    union
    {
      struct { uint64 xdim, ydim, zdim; };
      uint64 dim_[3];
    };
#endif
    
    /* Default constructor */
    Dimension() : xdim(0), ydim(0), zdim(0) {}
      
    /* Standard constructor */
    Dimension(uint64 x, uint64 y, uint64 z) :
      xdim(x), ydim(y), zdim(z) {}

    /*
      Universal explicit converter from any class to Dimension (as long as that class implements
      operator[]).
    */
    template <class C> explicit Dimension(const C& m) : 
      xdim(m[0]), ydim(m[1]), zdim(m[2]) {}

    Dimension& operator=(const Dimension& d)
    {
      xdim = d.xdim; ydim = d.ydim; zdim = d.zdim;
      return *this;
    }

    bool operator==(const Dimension& d) const
    {
      return (xdim == d.xdim) && (ydim == d.ydim) && (zdim == d.zdim);
    }

    bool operator!=(const Dimension& d) const
    {
      return !((*this)==d);
    }

    bool operator<(const Dimension& d) const
    {
      return (*this <= d) && (*this != d);
    }

    bool operator>(const Dimension& d) const
    {
      return (*this >= d) && (*this != d);
    }

    bool operator<=(const Dimension& d) const
    {
      return (xdim <= d.xdim) && (ydim <= d.ydim) && (zdim <= d.zdim);
    }

    bool operator>=(const Dimension& d) const
    {
      return (xdim >= d.xdim) && (ydim >= d.ydim) && (zdim >= d.zdim);
    }

    void setDim(uint64 x, uint64 y, uint64 z) { xdim = x; ydim = y; zdim = z; }

    /* Bracket operator with a constant return value. */
    uint64 operator[](int i) const 
      {
#ifdef UNION_NOT_SUPPORTED
	return (&xdim)[i];
#else
	return dim_[i]; 
#endif      
      }

    /* Bracket operator returning an l-value. */
    uint64& operator[](int i) 
      {
#ifdef UNION_NOT_SUPPORTED
	return (&xdim)[i];
#else 
	return dim_[i];
#endif 
      }

    bool isNull() const { return xdim == 0 && ydim == 0 && zdim == 0; }

    //returns the number of voxels for this dimension
    uint64 size() const { return xdim*ydim*zdim; }
  };
};

#endif
