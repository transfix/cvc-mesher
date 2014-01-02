#ifndef __LIBLBIE_H__
#define __LIBLBIE_H__

#include <iostream>
#include <VolMagick.h>

#include "octree.h"
#include "LBIE_geoframe.h"
#include <FastContouring.h>

namespace LBIE
{
  const float DEFAULT_ERR = 1.2501f;	    //-0.0001f    //99.99f		//0.0001f
  const float DEFAULT_ERR_IN = 0.0001f; 
  const float DEFAULT_IVAL = 0.5001f; //-0.0001f    //-0.5000f	//-0.0001f  //-1.0331		0.0001
  const float DEFAULT_IVAL_IN = 9.5001f;    //10.000f
  
  class Mesher
  {
  public:
    enum MeshType { SINGLE, TETRA, QUAD, HEXA, DOUBLE, TETRA2 };
    enum ImproveMethod { NO_IMPROVE, GEO_FLOW, EDGE_CONTRACT, JOE_LIU, MINIMAL_VOL, OPTIMIZATION };
    enum NormalType { BSPLINE_CONVOLUTION, CENTRAL_DIFFERENCE, BSPLINE_INTERPOLATION };
    enum ExtractionMethod { DUALLIB, FASTCONTOURING, LIBISOCONTOUR };

    Mesher(float isovalue = DEFAULT_IVAL,
	   float isovalue_in = DEFAULT_IVAL_IN,
	   float err = DEFAULT_ERR,
	   float err_in = DEFAULT_ERR_IN,
	   MeshType meshtype = SINGLE,
	   ImproveMethod improveMethod = NO_IMPROVE,
	   NormalType normtype = BSPLINE_CONVOLUTION,
	   ExtractionMethod extractionMethod = DUALLIB,
	   bool dual_contouring = false)
      : _isovalue(isovalue * -1.0), _isovalue_in(isovalue_in * -1.0),
        _err(err), _err_in(err_in), _meshType(meshtype),
        _improveMethod(improveMethod), _normalType(normtype),
        _extractionMethod(extractionMethod), _dual(dual_contouring),
        _octree_volume_updated(false)
   {
     _octree.set_isovalue(_isovalue);
     _octree.set_isovalue_in(_isovalue_in);
     _octree.setMeshType(int(_meshType));
     _octree.setNormalType(int(_normalType));
     if(_meshType == DOUBLE || _meshType == TETRA2)
       dual(true);
   }

   Mesher(const Mesher& mesher)
     : _isovalue(mesher._isovalue), _isovalue_in(mesher._isovalue_in),
       _err(mesher._err), _err_in(mesher._err_in), _meshType(mesher._meshType),
       _improveMethod(mesher._improveMethod), _normalType(mesher._normalType),
       _extractionMethod(mesher._extractionMethod), _dual(mesher._dual),
       _octree(mesher._octree), _geoframe(mesher._geoframe),
       _contourExtractor(mesher._contourExtractor),
       _octree_volume_updated(mesher._octree_volume_updated)
   {
     _octree.set_isovalue(_isovalue);
     _octree.set_isovalue_in(_isovalue_in);
     _octree.setMeshType(int(_meshType));
     _octree.setNormalType(int(_normalType));
     if(_meshType == DOUBLE || _meshType == TETRA2)
       _dual = true;
   }

   ~Mesher() {}

   Mesher& operator=(const Mesher& mesher)
   {
     _octree = mesher._octree;
     _contourExtractor = mesher._contourExtractor;
     isovalue(mesher.isovalue());
     isovalue_in(mesher.isovalue_in());
     err(mesher.err());
     err_in(mesher.err_in());
     meshType(mesher.meshType());
     improveMethod(mesher.improveMethod());
     normalType(mesher.normalType());
     extractionMethod(mesher.extractionMethod());
     dual(mesher.dual());
     mesh(mesher.mesh());
     _octree_volume_updated = mesher._octree_volume_updated;
     return *this;
   }

   void isovalue(float val) { _isovalue = val * -1.0; _octree.set_isovalue(_isovalue); }
   float isovalue() const { return _isovalue * -1.0; }
   void isovalue_in(float val) { _isovalue_in = val * -1.0; _octree.set_isovalue_in(_isovalue_in); }
   float isovalue_in() const { return _isovalue_in * -1.0; }

   void err(float val) { _err = val; }
   float err() const { return _err; }
   void err_in(float val) { _err_in = val; }
   float err_in() const { return _err_in; }

   void meshType(MeshType mt)
   {
     _meshType = mt;
     if(_meshType == DOUBLE || _meshType == TETRA2)
       _dual = true;
     _octree.setMeshType(int(_meshType));
   }
   MeshType meshType() const { return _meshType; }

   void improveMethod(ImproveMethod im) { _improveMethod = im; }
   ImproveMethod improveMethod() const { return _improveMethod; }

   void normalType(NormalType normtype)
   {
     _normalType = normtype;
     _octree.setNormalType(int(_normalType));
   }
   NormalType normalType() const { return _normalType; }

   void extractionMethod(ExtractionMethod et) { _extractionMethod = et; }
   ExtractionMethod extractionMethod() const { return _extractionMethod; }

   void dual(bool d) { if(_meshType != DOUBLE && _meshType != TETRA2) _dual = d; }
   bool dual() const { return _dual; }

   //mesh extraction... returns a reference to this class's geoframe
   void setVolume(const VolMagick::Volume& vol);
   geoframe& extractMesh(const VolMagick::Volume& vol);
   geoframe& extractMesh();

   //set mesh
   void mesh(const geoframe& g)
   { 
     _geoframe = g;
     meshType(MeshType(_geoframe.mesh_type));
   }
   geoframe& mesh() { return _geoframe; }
   const geoframe& mesh() const { return _geoframe; }

   //quality improve mesh
   geoframe& qualityImprove(unsigned int iterations = 1);

  protected:
   void set_octree_volume(const VolMagick::Volume& vol);
   void do_duallib_extract();
   void do_fastcontouring_extract();
   void do_libisocontour_extract();

   float _isovalue;
   float _isovalue_in;
   float _err;
   float _err_in;

   MeshType _meshType;
   ImproveMethod _improveMethod;
   NormalType _normalType;
   ExtractionMethod _extractionMethod;

   bool _dual; //if true, use isovalue/error interval

   Octree _octree; //where all the meshing happens
   geoframe _geoframe; //the mesh itself
   FastContouring::ContourExtractor _contourExtractor; //fast contour extractor for single tri meshes
   bool _octree_volume_updated; //flag for updating volume in octree before duallib mesh extraction
 };

}

#endif
