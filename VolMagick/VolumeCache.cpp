/* $Id: VolumeCache.cpp,v 1.4 2008/07/24 20:49:01 transfix Exp $ */

#include <stdlib.h>
#include <math.h>

#include "VolMagick.h"
#include "VolumeCache.h"

namespace VolMagick
{
  void VolumeCache::add(const VolumeFileInfo& vfi)
  {
    if(!vfi.isSet())
      throw ReadError("Cannot add an uninitialized VolumeFileInfo object to the cache.");
    
    /*
      Make sure that it has the same bounding box as the rest of the cache volumes.
      Else throw an exception.
    */
    for(std::set<VolumeFileInfo, dimcmp>::iterator cur = _volCacheInfo.begin();
	cur != _volCacheInfo.end(); cur++)
      if((*cur).boundingBox() != vfi.boundingBox())
	throw VolumePropertiesMismatch("All VolumeFileInfo objects in a VolumeCache object must have the same bounding box.");
    
    _volCacheInfo.insert(vfi);
  }

  Volume VolumeCache::get(const BoundingBox& requestRegion, unsigned int var, unsigned int time) const
  {
    Dimension curDim, closestDim;
    VolumeFileInfo closestCacheFile;
    Volume retVol;
    
    if(!requestRegion.isWithin(boundingBox()))
      throw SubVolumeOutOfBounds("Request region bounding box must be within the bounding box defined by the cache.");
    
    if(size()>0)
      {
	//initialize the closest cache file info to be the first volume in the cache
	closestCacheFile = *(_volCacheInfo.begin());
	/*
	closestDim = Dimension(uint64(((requestRegion.maxx-requestRegion.minx)/closestCacheFile.XSpan())+1.0),
			       uint64(((requestRegion.maxy-requestRegion.miny)/closestCacheFile.YSpan())+1.0),
			       uint64(((requestRegion.maxz-requestRegion.minz)/closestCacheFile.ZSpan())+1.0));
	*/
	/*
	closestDim = Dimension(ceil((requestRegion.maxx-requestRegion.minx)/closestCacheFile.XSpan())+1,
			       ceil((requestRegion.maxy-requestRegion.miny)/closestCacheFile.YSpan())+1,
			       ceil((requestRegion.maxz-requestRegion.minz)/closestCacheFile.ZSpan())+1);
	*/
	closestDim = Dimension(((requestRegion.maxx-requestRegion.minx)/closestCacheFile.XSpan())+1,
			       ((requestRegion.maxy-requestRegion.miny)/closestCacheFile.YSpan())+1,
			       ((requestRegion.maxz-requestRegion.minz)/closestCacheFile.ZSpan())+1);
	
	for(std::set<VolumeFileInfo, dimcmp>::iterator cur = _volCacheInfo.begin();
	    cur != _volCacheInfo.end(); cur++)
	  {
	    /*
	    curDim = Dimension(uint64(((requestRegion.maxx-requestRegion.minx)/(*cur).XSpan())+1.0),
			       uint64(((requestRegion.maxy-requestRegion.miny)/(*cur).YSpan())+1.0),
			       uint64(((requestRegion.maxz-requestRegion.minz)/(*cur).ZSpan())+1.0));
	    */
	    /*
	    curDim = Dimension(ceil((requestRegion.maxx-requestRegion.minx)/cur->XSpan())+1,
			       ceil((requestRegion.maxy-requestRegion.miny)/cur->YSpan())+1,
			       ceil((requestRegion.maxz-requestRegion.minz)/cur->ZSpan())+1);
	    */
	    curDim = Dimension(((requestRegion.maxx-requestRegion.minx)/cur->XSpan())+1,
			       ((requestRegion.maxy-requestRegion.miny)/cur->YSpan())+1,
			       ((requestRegion.maxz-requestRegion.minz)/cur->ZSpan())+1);
	    
	    //get the closest cache file to the _maxDimension
	    if(abs(int64(curDim.size()) - int64(_maxDimension.size())) <
	       abs(int64(closestDim.size()) - int64(_maxDimension.size())))
	      {
		closestCacheFile = *cur;
		closestDim = curDim;
	      }
	  }
	
	//extract the smallest subvolume within voxel boundaries
	readVolumeFile(retVol,
		       closestCacheFile.filename(),
		       var, time,
		       /*floor*/((requestRegion.minx-boundingBox().minx)/closestCacheFile.XSpan()),
		       /*floor*/((requestRegion.miny-boundingBox().miny)/closestCacheFile.YSpan()),
		       /*floor*/((requestRegion.minz-boundingBox().minz)/closestCacheFile.ZSpan()),
		       closestDim);
	
	//extract the exact subvolume and return a volume exactly _maxDimension in size
	if(_maxDimension.size() < closestDim.size())
	  retVol.sub(requestRegion,_maxDimension);
      }
    
    return retVol;
  }
};
