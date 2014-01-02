#ifndef __LBIE__FASTCONTOURING_H__
#define __LBIE__FASTCONTOURING_H__

#include <vector>
#include <VolMagick.h>

#include "Matrix.h"
#include "MarchingCubesBuffers.h"

namespace FastContouring
{
  //dead simple contour surface description
  struct TriSurf
  {
    std::vector<double> verts;
    std::vector<double> normals;
    std::vector<double> colors;
    std::vector<unsigned int> tris;
  };

  class ContourExtractor
  {
  public:
    ContourExtractor() {}
    ContourExtractor(const ContourExtractor& copy)
      : m_Data(copy.m_Data),
        m_Buffers(copy.m_Buffers),
        m_SaveMatrix(copy.m_SaveMatrix)
      {}
    ~ContourExtractor() {}
    
    ContourExtractor& operator=(const ContourExtractor& copy)
    {
      m_Data = copy.m_Data;
      m_Buffers = copy.m_Buffers;
      m_SaveMatrix = copy.m_SaveMatrix;
      return *this;
    }

    void setVolume(const VolMagick::Volume& vol);
    const VolMagick::Volume& getVolume() const { return m_Data; }
    
    TriSurf extractContour(double isovalue,
			   double R = 1.0, double G = 1.0, double B = 1.0);

  private:
    void classifyVertices(unsigned int k, unsigned int* cacheMemory, float isovalue) const;
    void getNormal(unsigned int i, unsigned int j, unsigned int k,
		   float& nx, float& ny, float& nz) const;
    unsigned int determineCase(/*unsigned char* data, */
			       unsigned int* offsetTable, unsigned int index
			       /*, GLfloat isovalue*/) const;

    VolMagick::Volume m_Data;
    // buffers used to speed up marching cubes
    MarchingCubesBuffers m_Buffers;
    Matrix m_SaveMatrix;
  };
			 
};

#endif
