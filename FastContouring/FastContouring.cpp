#include "FastContouring.h"
#include "cubes.h"
#include "ContourGeometry.h"

namespace FastContouring
{

  void ContourExtractor::setVolume(const VolMagick::Volume& vol)
  {
    m_Data = vol;
    m_Buffers.allocateEdgeBuffers(vol.XDim(), vol.YDim());
    
    // set up the save matrix
    m_SaveMatrix.reset();
    
    // The vertices output by the extraction algorithm are in index space.  The save matrix
    // converts vertices from index space to object space... point = min + idx * ((max - min)/(dim-1))
    m_SaveMatrix.preMultiplication(Matrix::scale((vol.XMax() - vol.XMin())/double(vol.XDim()-1),
						 (vol.YMax() - vol.YMin())/double(vol.YDim()-1),
						 (vol.ZMax() - vol.ZMin())/double(vol.ZDim()-1)));
    
    m_SaveMatrix.preMultiplication(Matrix::translation(vol.XMin(),
						       vol.YMin(),
						       vol.ZMin()));
  }

  /*********** Utility functions for extractContour() **************/
  void ContourExtractor::classifyVertices(unsigned int k, unsigned int* cacheMemory, float isovalue) const
  {
    unsigned int i;
    unsigned int sourceOffset = k*m_Data.XDim()*m_Data.YDim();
    unsigned int destOffset = 0;
    
    for (i=0; i<m_Data.XDim()*m_Data.YDim(); i++) {
      cacheMemory[destOffset++] = ((float)m_Data(sourceOffset++) < isovalue);
    }
  }

  void ContourExtractor::getNormal(unsigned int i, unsigned int j, unsigned int k,
				   float& nx, float& ny, float& nz) const
  {
    double dx,dy,dz,length;
    int negZOffset, negYOffset, negXOffset;
    int posZOffset, posYOffset, posXOffset;
    double scaleX = 1.0, scaleY = 1.0, scaleZ = 1.0;
    unsigned int inputindex;
    if (k==0) { // border offset
      negZOffset = 0; posZOffset =  m_Data.XDim()*m_Data.YDim();
    }
    else if (k==m_Data.ZDim()-1) { // border offset
      negZOffset = -(int)(m_Data.XDim()*m_Data.YDim()); posZOffset = 0;
    }
    else { // normal offset
      negZOffset = -(int)(m_Data.XDim()*m_Data.YDim()); posZOffset =  m_Data.XDim()*m_Data.YDim();
      scaleZ = 0.5;
    }
    
    if (j==0) { // border offset
      negYOffset = 0; posYOffset =  m_Data.XDim();
    }
    else if (j==m_Data.YDim()-1) { // border offset
      negYOffset = -(int)m_Data.XDim(); posYOffset = 0;
    }				
    else { // normal offset
      negYOffset = -(int)m_Data.XDim(); posYOffset =  m_Data.XDim();
      scaleY = 0.5;
    }
    
    if (i==0) { // border offset
      negXOffset = 0; posXOffset =  1;
    }
    else if (i==m_Data.XDim()-1) { // border offset
      negXOffset = -1; posXOffset = 0;
    }				
    else { // normal offset
      negXOffset = -1; posXOffset =  1;
      scaleX = 0.5;
    }
    
    inputindex = k*m_Data.XDim()*m_Data.YDim()+j*m_Data.XDim()+i;
    dx = (double)(m_Data(inputindex+negXOffset) - m_Data(inputindex+posXOffset)) * scaleX;
    dy = (double)(m_Data(inputindex+negYOffset) - m_Data(inputindex+posYOffset)) * scaleY;
    dz = (double)(m_Data(inputindex+negZOffset) - m_Data(inputindex+posZOffset)) * scaleZ;
    length = sqrt(dx*dx+dy*dy+dz*dz);
    if (length!=0) {
      nx = (float)(dx/length);
      ny = (float)(dy/length);
      nz = (float)(dz/length);
    }
  }

  static inline void isCached(bool* cachedTable, unsigned int i, unsigned int j, unsigned int k)
  {
    // bool ret = false;
    // ret = j!=0 || k!=0;

    // edge 0
    cachedTable[0] = j!=0 || k!=0;

    // edge 1
    cachedTable[1] = j!=0;

    // edge 2
    cachedTable[2] = j!=0;

    // edge 3
    cachedTable[3] = i!=0 || j!=0;

    // edge 4
    cachedTable[4] = k!=0;

    // edge 5
    cachedTable[5] = false;

    // edge 6
    cachedTable[6] = false;

    // edge 7
    cachedTable[7] = i!=0;

    // edge 8
    cachedTable[8] = k!=0 || i!=0;

    // edge 9
    cachedTable[9] = k!=0;

    // edge 10
    cachedTable[10] = i!=0;

    // edge 11
    cachedTable[11] = false;
  }

  inline unsigned int ContourExtractor::determineCase(/*unsigned char* data, */unsigned int* offsetTable, unsigned int index /*, GLfloat isovalue*/) const
  {
    unsigned int cubeCase = 0;

	// determine the marching cube case
/*
	if ((GLfloat)data[index+offsetTable[0]] < isovalue) cubeCase |= 1<<0;
	if ((GLfloat)data[index+offsetTable[1]] < isovalue) cubeCase |= 1<<1;
	if ((GLfloat)data[index+offsetTable[2]] < isovalue) cubeCase |= 1<<2;
	if ((GLfloat)data[index+offsetTable[3]] < isovalue) cubeCase |= 1<<3;
	if ((GLfloat)data[index+offsetTable[4]] < isovalue) cubeCase |= 1<<4;
	if ((GLfloat)data[index+offsetTable[5]] < isovalue) cubeCase |= 1<<5;
	if ((GLfloat)data[index+offsetTable[6]] < isovalue) cubeCase |= 1<<6;
	if ((GLfloat)data[index+offsetTable[7]] < isovalue) cubeCase |= 1<<7;
*/
	/*
	cubeCase |= ((GLfloat)data[index+offsetTable[0]] < isovalue)<<0;
	cubeCase |= ((GLfloat)data[index+offsetTable[1]] < isovalue)<<1;
	cubeCase |= ((GLfloat)data[index+offsetTable[2]] < isovalue)<<2;
	cubeCase |= ((GLfloat)data[index+offsetTable[3]] < isovalue)<<3;
	cubeCase |= ((GLfloat)data[index+offsetTable[4]] < isovalue)<<4;
	cubeCase |= ((GLfloat)data[index+offsetTable[5]] < isovalue)<<5;
	cubeCase |= ((GLfloat)data[index+offsetTable[6]] < isovalue)<<6;
	cubeCase |= ((GLfloat)data[index+offsetTable[7]] < isovalue)<<7;
	*/

    cubeCase |= (m_Buffers.m_VertClassifications[0][index+offsetTable[0]])<<0;
    cubeCase |= (m_Buffers.m_VertClassifications[0][index+offsetTable[1]])<<1;
    cubeCase |= (m_Buffers.m_VertClassifications[1][index+offsetTable[2]])<<2;
    cubeCase |= (m_Buffers.m_VertClassifications[1][index+offsetTable[3]])<<3;
    cubeCase |= (m_Buffers.m_VertClassifications[0][index+offsetTable[4]])<<4;
    cubeCase |= (m_Buffers.m_VertClassifications[0][index+offsetTable[5]])<<5;
    cubeCase |= (m_Buffers.m_VertClassifications[1][index+offsetTable[6]])<<6;
    cubeCase |= (m_Buffers.m_VertClassifications[1][index+offsetTable[7]])<<7;
    return cubeCase;
  }

  static inline void buildOffsetTable(unsigned int* table, unsigned int width, unsigned int height)
  {
    // determine the offset of each vertex of the cube give the width and height of the data
    unsigned int c;
    for (c=0; c<8; c++)
      {
	table[c] = verts[c][2] + width*verts[c][1] + width*height*verts[c][0];
      }
  }

  static inline void buildWithinSliceOffsetTable(unsigned int* table, unsigned int width)
  {
    // determine the offset of each vertex of the cube give the width and height of the data
    unsigned int c;
    for (c=0; c<8; c++)
      {
	table[c] = verts[c][2] + width*verts[c][1];
      }
  }

  static inline void buildEdgeCacheOffsetTable(unsigned int * table, unsigned int width)
  {
    // calculates where edges are cached in their respective table
    // is an offset from k*width*height + j*width + i
    
    // edge 0
    table[0] = 0;
    
    // edge 1
    table[1] = 1;
    
    // edge 2
    table[2] = 0;
    
    // edge 3
    table[3] = 0;
    
    // edge 4
    table[4] = width;
    
    // edge 5
    table[5] = width+1;
    
    // edge 6
    table[6] = width;
    
    // edge 7
    table[7] = width;
    
    // edge 8
    table[8] = 0;
    
    // edge 9
    table[9] = 1;
    
    // edge 10
    table[10] = 0;
    
    // edge 11
    table[11] = 1;
  }

#if 0
  static inline void computeVertFromOffset(float& vx, float& vy, float& vz, unsigned int offset, unsigned int width, unsigned int height)
  {
    vx = (float)(offset % (width));
    vy = (float)((offset/width) % height);
    vz = (float)(offset / (width*height));
  }
#endif

#if 0
  static inline void calcOffsets(unsigned int* offsets, unsigned int width, unsigned int height)
  {
    // calculate the table of offsets given the width and height
    offsets[0] = verts[0][0]*width*height + verts[0][1]*width + verts[0][2];
    offsets[1] = verts[1][0]*width*height + verts[1][1]*width + verts[1][2];
    offsets[2] = verts[2][0]*width*height + verts[2][1]*width + verts[2][2];
    offsets[3] = verts[3][0]*width*height + verts[3][1]*width + verts[3][2];
    offsets[4] = verts[4][0]*width*height + verts[4][1]*width + verts[4][2];
    offsets[5] = verts[5][0]*width*height + verts[5][1]*width + verts[5][2];
    offsets[6] = verts[6][0]*width*height + verts[6][1]*width + verts[6][2];
    offsets[7] = verts[7][0]*width*height + verts[7][1]*width + verts[7][2];
    offsets[8] = verts[8][0]*width*height + verts[8][1]*width + verts[8][2];
  }
#endif

  TriSurf ContourExtractor::extractContour(double isovalue,
					   double R, double G, double B)
  {
    TriSurf result;
    ContourGeometry contourGeometry;
    contourGeometry.allocateVertexBuffers(m_Data.XDim()*m_Data.YDim());
    contourGeometry.allocateTriangleBuffers(m_Data.XDim()*m_Data.YDim());
    contourGeometry.setUseColors(true);
    contourGeometry.setIsovalue(isovalue);
    unsigned int numCellsX = m_Data.XDim()-1;
    unsigned int numCellsY = m_Data.YDim()-1;
    unsigned int numCellsZ = m_Data.ZDim()-1;

    bool cacheAvailable[12];
    
    unsigned int offsetTable[8];
    unsigned int withinSliceOffsetTable[8];
    unsigned int cacheOffsetTable[12];
    buildOffsetTable(offsetTable, m_Data.XDim(), m_Data.YDim());
    buildWithinSliceOffsetTable(withinSliceOffsetTable, m_Data.XDim());
    buildEdgeCacheOffsetTable(cacheOffsetTable, m_Data.XDim());

    unsigned int cubeCase;

    unsigned int i,j,k, voxelOffset, offsetInSlice, newEdge, edgeCounter, triangleCounter, edgeID;

    // the new vertex
    float v1x, v1y, v1z, v2x, v2y, v2z, den1, den2;
    float n1x, n1y, n1z, n2x, n2y, n2z;

    // edges involved in a triangle
    unsigned int e1,e2,e3;
    unsigned int rowStart, inSliceRowStart;

    classifyVertices(0, m_Buffers.m_VertClassifications[0].get(), isovalue);

    for (k=0; k<numCellsZ; k++) { // for every slice
      classifyVertices(k+1, m_Buffers.m_VertClassifications[1].get(), isovalue);

      for (j=0; j<numCellsY; j++) { // for every row
	rowStart = m_Data.XDim()*m_Data.YDim()*k+m_Data.XDim()*j;
	inSliceRowStart = m_Data.XDim()*j;
	for (i=0; i<numCellsX; i++) { // for every cell
	  if (i==0 || i==1 || i==numCellsX-1) {
	    // isCached only changes at the start and end of a scan line
	    // as well as the second voxel of every scan line
	    isCached(cacheAvailable, i,j,k/*,m_Data.XDim(),m_Data.YDim()*/);
	  }
	  voxelOffset = rowStart+i;
	  offsetInSlice = inSliceRowStart+i;
	  cubeCase = determineCase(/*m_Data, */withinSliceOffsetTable, offsetInSlice/*, isovalue*/);
	  if (cubeCase!=0 && cubeCase!=255) {
	    // for each edge involved
	    for (edgeCounter=0; edgeCounter<cubeedges[cubeCase][0]; edgeCounter++) {
	      edgeID = cubeedges[cubeCase][edgeCounter+1];
	      // if the edge isnt cached yet, cache it
	      if (!cacheAvailable[edgeID]) {

		// add the edge and get its index
		/*computeVertFromOffset(v1x, v1y, v1z, offsetTable[edges[edgeID][0]],
		  m_Width, m_Height);*/
		/*computeVertFromOffset(v2x, v2y, v2z, offsetTable[edges[edgeID][1]],
		  m_Width, m_Height);*/
		v1x = (float)(i + verts[edges[edgeID][0]][2]);
		v1y = (float)(j + verts[edges[edgeID][0]][1]);
		v1z = (float)(k + verts[edges[edgeID][0]][0]);
		v2x = (float)(i + verts[edges[edgeID][1]][2]);
		v2y = (float)(j + verts[edges[edgeID][1]][1]);
		v2z = (float)(k + verts[edges[edgeID][1]][0]);
		den1 = m_Data(voxelOffset+offsetTable[edges[edgeID][0]]);
		den2 = m_Data(voxelOffset+offsetTable[edges[edgeID][1]]);
		getNormal(i + verts[edges[edgeID][0]][2], 
			  j + verts[edges[edgeID][0]][1], 
			  k + verts[edges[edgeID][0]][0], 
			  n1x, n1y, n1z);
		getNormal(i + verts[edges[edgeID][1]][2], 
			  j + verts[edges[edgeID][1]][1], 
			  k + verts[edges[edgeID][1]][0], 
			  n2x, n2y, n2z);
		// no color
		//contourGeometry->addEdge(newEdge, v1x, v1y, v1z, n1x, n1y, n1z,
		//	v2x, v2y, v2z, n2x, n2y, n2z, den1, den2);
		// color
		contourGeometry.addEdge(newEdge, v1x, v1y, v1z, n1x, n1y, n1z,
					 R,G,B, v2x, v2y, v2z, n2x, n2y, n2z, R,G,B, den1, den2);

		// The location in the cache is determined by two table lookups:
		// FromEdgeToChacheLookup[edgeID] finds which table the edge is cached in
		// offsetInSlice+cacheOffsetTable[edgeID] determines the location in that table
		m_Buffers.m_EdgeCaches[FromEdgeToChacheLookup[edgeID]][offsetInSlice+cacheOffsetTable[edgeID]] = newEdge;

	      }
	    }
	    // All appropriate edges are now cached
	    // Build the triangles, using indexes from the cache
	    // for each triangle
					
	    for (triangleCounter=0; triangleCounter<cubes[cubeCase][0]; triangleCounter++) {
	      e1 = cubes[cubeCase][triangleCounter*3+1];
	      e2 = cubes[cubeCase][triangleCounter*3+2];
	      e3 = cubes[cubeCase][triangleCounter*3+3];
						
	      contourGeometry.addTriangle(
					   m_Buffers.m_EdgeCaches[FromEdgeToChacheLookup[e1]][offsetInSlice+cacheOffsetTable[e1]], 
					   m_Buffers.m_EdgeCaches[FromEdgeToChacheLookup[e2]][offsetInSlice+cacheOffsetTable[e2]], 
					   m_Buffers.m_EdgeCaches[FromEdgeToChacheLookup[e3]][offsetInSlice+cacheOffsetTable[e3]]);
	    }

	  }

	}

      }
      // swap the edge buffers
      m_Buffers.swapEdgeBuffers();
    }

    int nextVert = 0, nextTri = 0;
    contourGeometry.addToGeometry(result, m_SaveMatrix, nextVert, nextTri);

    return result;
  }
}
