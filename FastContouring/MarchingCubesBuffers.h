/*
  Copyright 2002-2003 The University of Texas at Austin
  
	Authors: Anthony Thane <thanea@ices.utexas.edu>
	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of Volume Rover.

  Volume Rover is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Volume Rover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with iotree; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

// MarchingCubesBuffers.h: interface for the MarchingCubesBuffers class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MARCHINGCUBESBUFFERS_H__2F1981CF_3DB7_402D_B078_AB961D4B6189__INCLUDED_)
#define AFX_MARCHINGCUBESBUFFERS_H__2F1981CF_3DB7_402D_B078_AB961D4B6189__INCLUDED_

#include <boost/array.hpp>
#include <boost/shared_array.hpp>

namespace FastContouring
{
///\class MarchingCubesBuffers MarchingCubesBuffers.h
///\author Anthony Thane
///\brief The MarchingCubesBuffers class is used for caching during isocontour
/// extraction.
  class MarchingCubesBuffers  
  {
  public:
    MarchingCubesBuffers();
    MarchingCubesBuffers(const MarchingCubesBuffers& copy);
    virtual ~MarchingCubesBuffers();

    MarchingCubesBuffers& operator=(const MarchingCubesBuffers& copy);

///\fn bool allocateEdgeBuffers(unsigned int width, unsigned int height)
///\brief Allocates memory for a cell layer's worth of edges.
///\param width The width of the cell layer
///\param height The height of the cell layer
///\return A bool indicating success or failure
    bool allocateEdgeBuffers(unsigned int width, unsigned int height);
///\fn void destroyEdgeBuffers()
///\brief Deallocates previously allocated edge buffers
    void destroyEdgeBuffers();
///\fn void swapEdgeBuffers()
///\brief Swaps edge buffers (There are two sets of vertices for each layer, one
/// of them gets reused on the next layer of the volume. This function does the
/// swapping for that)
    void swapEdgeBuffers();
// cached edges to avoid recomputation
    //unsigned int* m_EdgeCaches[5];
    typedef boost::array<boost::shared_array<unsigned int>, 5> EdgeCacheArray;
    EdgeCacheArray m_EdgeCaches;
    //unsigned int* m_VertClassifications[2];
    typedef boost::array<boost::shared_array<unsigned int>, 2> VertexClassificationsArray;
    VertexClassificationsArray m_VertClassifications;

  protected:
    void setDefaults();
    
    bool forceAllocateEdgeBuffers(unsigned int width, unsigned int height);
    
    unsigned int m_AmountAllocated;
  };
}

#endif // !defined(AFX_MARCHINGCUBESBUFFERS_H__2F1981CF_3DB7_402D_B078_AB961D4B6189__INCLUDED_)
