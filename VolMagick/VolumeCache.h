/* $Id: VolumeCache.h,v 1.2 2008/02/01 20:12:12 transfix Exp $ */

#ifndef __VOLUMECACHE_H__
#define __VOLUMECACHE_H__

#include <set>

#include "VolMagick.h"

namespace VolMagick
{
  /*
    Compare the volume of each set of voxels
  */
  struct dimcmp
  {
    bool operator()(const VolumeFileInfo& vfi1, const VolumeFileInfo& vfi2)
    {
      return vfi1.dimension().size() < vfi2.dimension().size();
    }
  };

  class VolumeCache
  {
  public:
    VolumeCache(const Dimension& d = Dimension(128,128,128)) 
      : _maxDimension(d) {}
    
    VolumeCache(const VolumeCache& vc) 
      : _volCacheInfo(vc._volCacheInfo), _maxDimension(vc.maxDimension()) {}
    
    ~VolumeCache() {}

    VolumeCache& operator=(const VolumeCache& vc)
      {
	_volCacheInfo = vc._volCacheInfo;
	return *this;
      }

    void maxDimension(const Dimension& d) { _maxDimension = d; }
    Dimension& maxDimension() { return _maxDimension; }
    const Dimension& maxDimension() const { return _maxDimension; }

    void add(const VolumeFileInfo& vfi); //add a volume to the cache

    void clear() { _volCacheInfo.clear(); } //clear the cache

    unsigned int size() const { return _volCacheInfo.size(); } //return the number of volumes in the cache

    BoundingBox boundingBox() const
    {
      if(size()>0)
	return (*(_volCacheInfo.begin())).boundingBox();
      else
	return BoundingBox();
    }

    /*
      Return a volume for the request region.
    */
    Volume get(const BoundingBox& requestRegion, unsigned int var = 0, unsigned int time = 0) const;

  private:
    std::set<VolumeFileInfo, dimcmp> _volCacheInfo;
    Dimension _maxDimension;
  };
};

#endif
