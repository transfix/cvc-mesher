/* $Id: BoundingBox.h,v 1.3 2008/07/18 14:02:55 transfix Exp $ */

#ifndef __VOLMAGICK_BOUNDINGBOX_H__
#define __VOLMAGICK_BOUNDINGBOX_H__

#include <stdio.h>
#include <algorithm>

#ifdef WIN32
#define SNPRINTF _snprintf
#else
#define SNPRINTF snprintf
#endif

// If your compiler complains the "The class "VolMagick::GenericBoundingBox<T>" has no member "minx"."
// Add your architecture Q_OS_XXXX flag (see qglobal.h) in this list.
//#if defined (Q_OS_IRIX) || defined (Q_OS_AIX) || defined (Q_OS_HPUX)
# define UNION_NOT_SUPPORTED
//#endif

namespace VolMagick
{
  /*
    GenericBoundingBox min/max values are inclusive.
  */
  template <typename T>
  class GenericBoundingBox
  {
  public:
    /* The internal data representation is public. */

#if defined (DOXYGEN) || defined (UNION_NOT_SUPPORTED)
    T minx, miny, minz;
#else
    union
    {
      struct { T minx, miny, minz; };
      T min_[3];
    };
#endif
    
#if defined (DOXYGEN) || defined (UNION_NOT_SUPPORTED)
    T maxx, maxy, maxz;
#else
    union
    {
      struct { T maxx, maxy, maxz; };
      T max_[3];
    };
#endif    
    
    /* default constructor */
    GenericBoundingBox()
      {
	setMin(T(0),T(0),T(0));
	setMax(T(0),T(0),T(0));
      }
      
      /* standard constructor */
    GenericBoundingBox(T minx_, T miny_, T minz_, T maxx_, T maxy_, T maxz_)
      { 
	setMin(minx_,miny_,minz_);
	setMax(maxx_,maxy_,maxz_);
	checkBounds();
      }

    GenericBoundingBox(const Dimension& dimension)
      { 
	setMin(T(0),T(0),T(0));
	setMax(T(dimension.xdim-1),T(dimension.ydim-1),T(dimension.zdim-1));
	checkBounds();
      }
	
    /*
      Universal explicit converter from any class to GenericBoundingBox (as long as that class implements
      operator[]).
    */
    template <class C> explicit GenericBoundingBox(const C& m)
      { 
	setMin(T(m[0]),T(m[1]),T(m[2]));
	setMax(T(m[3]),T(m[4]),T(m[5]));
	checkBounds(); 
      }

    GenericBoundingBox<T>& operator=(const GenericBoundingBox<T>& b)
    {
      setMin(b.minx,b.miny,b.minz);
      setMax(b.maxx,b.maxy,b.maxz);
      return *this;
    }

    void setMin(T minx_, T miny_, T minz_)
    { 
      minx=minx_; miny=miny_; minz=minz_;
    }

    void setMax(T maxx_, T maxy_, T maxz_)
    {
      maxx=maxx_; maxy=maxy_; maxz=maxz_;
    }

    
#ifdef UNION_NOT_SUPPORTED
# define REALMAX (&maxx)
# define REALMIN (&minx)
#else
# define REALMAX max_
# define REALMIN min_
#endif

    /* Bracket operator with a constant return value. */
    T operator[](int i) const { return i>=3 ? REALMAX[i-3] : REALMIN[i]; }

    /* Bracket operator returning an l-value. */
    T& operator[](int i) { return i>=3 ? REALMAX[i-3] : REALMIN[i]; }

#undef REALMAX
#undef REALMIN

    /*
      Union operator
    */

    GenericBoundingBox<T> operator+(const GenericBoundingBox<T>& rhs) const
      {
	if(rhs.isNull() && !isNull()) return *this; /* if one of the boxes are null, 
						       the result of this operation should be the non-null box */
	if(isNull() && !rhs.isNull()) return rhs;
	if(rhs.isNull() && isNull()) return *this;

	GenericBoundingBox<T> ret;
	ret.minx = MIN(minx,rhs.minx);
	ret.miny = MIN(miny,rhs.miny);
	ret.minz = MIN(minz,rhs.minz);
	ret.maxx = MAX(maxx,rhs.maxx);
	ret.maxy = MAX(maxy,rhs.maxy);
	ret.maxz = MAX(maxz,rhs.maxz);

	return ret;
      }

    GenericBoundingBox<T>& operator+=(const GenericBoundingBox<T>& rhs)
      {
	if(rhs.isNull() && !isNull()) return *this;
	if(isNull() && !rhs.isNull())
	  {
	    (*this) = rhs;
	    return *this;
	  }
	if(rhs.isNull() && isNull()) return *this;

	minx = MIN(minx,rhs.minx);
	miny = MIN(miny,rhs.miny);
	minz = MIN(minz,rhs.minz);
	maxx = MAX(maxx,rhs.maxx);
	maxy = MAX(maxy,rhs.maxy);
	maxz = MAX(maxz,rhs.maxz);
	
	return *this;
      }

    /*
      Intersection operator
    */
    GenericBoundingBox<T> operator-(const GenericBoundingBox<T>& rhs) const
      {
	if(rhs.isNull()) return rhs; /* if one of the boxes are null, 
					the result of this operation should be null */
	if(isNull()) return *this;
	GenericBoundingBox<T> ret;
	ret.minx = MAX(minx,rhs.minx);
	ret.miny = MAX(miny,rhs.miny);
	ret.minz = MAX(minz,rhs.minz);
	ret.maxx = MIN(maxx,rhs.maxx);
	ret.maxy = MIN(maxy,rhs.maxy);
	ret.maxz = MIN(maxz,rhs.maxz);

	/*
	  check to see if there is no intersection.  If there isn't, set ret
	  to null (else it will cause an exception on future operations).
	*/
	if(ret.minx > ret.maxx ||
	   ret.miny > ret.maxy ||
	   ret.minz > ret.maxz)
	  ret.minx = ret.maxx = ret.miny = ret.maxy = ret.minz = ret.maxz = T(0);

	return ret;
      }

    GenericBoundingBox<T>& operator-=(const GenericBoundingBox<T>& rhs)
      {
	if(rhs.isNull())
	  {
	    *this = rhs;
	    return *this;
	  }
	if(isNull()) return *this;
	minx = MAX(minx,rhs.minx);
	miny = MAX(miny,rhs.miny);
	minz = MAX(minz,rhs.minz);
	maxx = MIN(maxx,rhs.maxx);
	maxy = MIN(maxy,rhs.maxy);
	maxz = MIN(maxz,rhs.maxz);

	/*
	  check to see if there is no intersection.  If there isn't, set ret
	  to null (else it will cause an exception on future operations).
	*/
	if(minx > maxx ||
	   miny > maxy ||
	   minz > maxz)
	  minx = maxx = miny = maxy = minz = maxz = T(0);

	return *this;
      }

    bool isWithin(const GenericBoundingBox<T>& b) const
    {
      if(b.isNull()) return false;
      if(minx >= b.minx && miny >= b.miny && minz >= b.minz &&
	 maxx <= b.maxx && maxy <= b.maxy && maxz <= b.maxz)
	return true;
      return false;
    }

    bool contains(double x, double y, double z) const
    {
      if(isNull()) return false;
      if(x >= minx && y >= miny && z >= minz &&
	 x <= maxx && y <= maxy && z <= maxz)
	return true;
      return false;
    }

    bool operator==(const GenericBoundingBox<T>& b) const
    {
      return isWithin(b) && b.isWithin(*this);
    }

    bool operator!=(const GenericBoundingBox<T>& b) const
    {
      return !((*this)==b);
    }

    bool isNull() const
    {
      return fabs(volume()) >= 0 && fabs(volume()) <= 0;
    }

    double volume() const
    {
      checkBounds();
      return double(maxx-minx)*double(maxy-miny)*double(maxz-minz);
    }

    T XMax() const { return maxx; }
    T XMin() const { return minx; }
    T YMax() const { return maxy; }
    T YMin() const { return miny; }
    T ZMax() const { return maxz; }
    T ZMin() const { return minz; }

    void normalize()
    {
      *this = GenericBoundingBox<T>(std::min(minx,maxx),
				    std::min(miny,maxy),
				    std::min(minz,maxz),
				    std::max(minx,maxx),
				    std::max(miny,maxy),
				    std::max(minz,maxz));
    }

  private:
    void checkBounds() const throw(InvalidBoundingBox)
    {
      char buf[256];
      if(minx > maxx) 
	{
	  SNPRINTF(buf,256,"minx: %f, maxx: %f",double(minx),double(maxx));
	  throw InvalidBoundingBox(buf);
	}
      else if(miny > maxy)
	{
	  SNPRINTF(buf,256,"miny: %f, maxy: %f",double(miny),double(maxy));
	  throw InvalidBoundingBox(buf);
	}
      else if(minz > maxz)
	{
	  SNPRINTF(buf,256,"minz: %f, maxz: %f",double(minz),double(maxz));
	  throw InvalidBoundingBox(buf);
	}
    }
  };

  typedef GenericBoundingBox<double> BoundingBox; // object space
};

#endif
