/* $Id: Exceptions.h,v 1.4 2008/08/15 21:53:04 transfix Exp $ */

#ifndef __VOLMAGICK_EXCEPTIONS_H__
#define __VOLMAGICK_EXCEPTIONS_H__

#include <boost/format.hpp>
#include <exception>
#include <string>

namespace VolMagick
{
  /***** Exceptions ****/
  class Exception : public std::exception
  {
  public:
    Exception() {}
    virtual ~Exception() throw() {}
    virtual const std::string& what_str() const throw () = 0;
    virtual const char *what () const throw()
    {
      return what_str().c_str();
    }
  };
  
#define VOLMAGICK_DEF_EXCEPTION(name) \
  class name : public VolMagick::Exception \
  { \
  public: \
    name () : _msg("VolMagick::"#name) {} \
    name (const std::string& msg) : \
      _msg(boost::str(boost::format("VolMagick::" #name " exception: %1%") % msg)) {} \
    virtual ~name() throw() {} \
    virtual const std::string& what_str() const throw() { return _msg; } \
  private: \
    std::string _msg; \
  }

  VOLMAGICK_DEF_EXCEPTION(ReadError);
  VOLMAGICK_DEF_EXCEPTION(WriteError);
  VOLMAGICK_DEF_EXCEPTION(InvalidRawIVHeader);
  VOLMAGICK_DEF_EXCEPTION(InvalidRawVHeader);
  VOLMAGICK_DEF_EXCEPTION(InvalidMRCHeader);
  VOLMAGICK_DEF_EXCEPTION(InvalidMRCFile);
  VOLMAGICK_DEF_EXCEPTION(InvalidSpiderFile);
  VOLMAGICK_DEF_EXCEPTION(UnsupportedMRCFile);
  VOLMAGICK_DEF_EXCEPTION(MemoryAllocationError);
  VOLMAGICK_DEF_EXCEPTION(SubVolumeOutOfBounds);
  VOLMAGICK_DEF_EXCEPTION(InvalidBoundingBox);
  VOLMAGICK_DEF_EXCEPTION(UnsupportedVolumeFileType);
  VOLMAGICK_DEF_EXCEPTION(IndexOutOfBounds);
  VOLMAGICK_DEF_EXCEPTION(NullDimension);
  VOLMAGICK_DEF_EXCEPTION(VolumePropertiesMismatch);
};

#endif

