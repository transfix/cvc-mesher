#include <iostream>
#include <algorithm>

#include "LBIE_Mesher.h"

//libcontour headers
#include <contour.h>
#include <datasetreg3.h>

namespace LBIE
{
  static inline unsigned int upToPowerOfTwo(unsigned int value)
  {
    unsigned int c = 0;
    unsigned int v = value;
  
    // round down to nearest power of two
    while (v>1) {
      v = v>>1;
      c++;
    }
  
    // if that isn't exactly the original value
    if ((v<<c)!=value) {
      // return the next power of two
      return (v<<(c+1));
    }
    else {
      // return this power of two
      return (v<<c);
    }
  }

  void Mesher::setVolume(const VolMagick::Volume& vol)
  {
    using namespace std;

    cerr << "LBIE::Mesher::setVolume(): set volume: "; 
    VolMagick::Volume localvol(vol);
    _contourExtractor.setVolume(vol);
    cerr << "done" << endl;

    _octree_volume_updated = false;
  }

  void Mesher::set_octree_volume(const VolMagick::Volume &vol)
  {
    using namespace std;

    VolMagick::Volume localvol(vol);
    
    cerr << "LBIE::Mesher::set_octree_volume(): set volume: ";
    //LBIE requires volumes with dimension (2^n+1)^3
    unsigned int dim[3] = { localvol.XDim(), localvol.YDim(), localvol.ZDim() };
    unsigned int maxdim = *std::max_element(dim,dim+3);
    if((dim[0] != dim[1]) ||
       (dim[0] != dim[2]) ||
       upToPowerOfTwo(maxdim-1) != (maxdim-1))
      localvol.resize(VolMagick::Dimension(upToPowerOfTwo(maxdim)+1,
					   upToPowerOfTwo(maxdim)+1,
					   upToPowerOfTwo(maxdim)+1));
    _octree.setVolume(localvol);
    cerr << "done" << endl;

    _octree_volume_updated = true;
  }

  void Mesher::do_duallib_extract()
  {
    using namespace std;

    cerr << "LBIE::Mesher::extractMesh(): DUALLIB mesh extraction" << endl;
    
    _geoframe.setSpan(_octree.spans[0],_octree.spans[1],_octree.spans[2]);
    _geoframe.setMin(_octree.minext[0],_octree.minext[1],_octree.minext[2]);
    
    if(!_octree_volume_updated)
      set_octree_volume(_contourExtractor.getVolume());
    
    if(!_dual)
      {
	_octree.collapse();
	cerr << "LBIE::Mesher::extractMesh(): finished collapse" << endl;
	_octree.compute_qef();
	cerr << "LBIE::Mesher::extractMesh(): finished compute qef" << endl;
	_octree.traverse_qef(_err);
	cerr << "LBIE::Mesher::extractMesh(): finished traverse qef" << endl;
	_octree.mesh_extract(_geoframe,_err);
	cerr << "LBIE::Mesher::extractMesh(): finished mesh extraction" << endl;
      }
    else
      {
	_octree.collapse_interval();
	cerr << "LBIE::Mesher::extractMesh(): finished collapse interval" << endl;
	_octree.compute_qef_interval();
	cerr << "LBIE::Mesher::extractMesh(): finished compute qef interval" << endl;
	_octree.traverse_qef_interval(_err,_err_in);
	cerr << "LBIE::Mesher::extractMesh(): finished traverse qef interval" << endl;
	_octree.mesh_extract(_geoframe,_err);
	cerr << "LBIE::Mesher::extractMesh(): finished mesh extraction" << endl;
      }
    
    _octree.quality_improve(_geoframe,_improveMethod);
    _geoframe.mesh_type = geoframe::GEOTYPE(_meshType);
    
    //re-adjust vertices so they're back to their original location
    for(int i=0; i<_geoframe.numverts; i++)
      for(int j=0; j<3; j++)
	_geoframe.verts[i][j] =
	  _octree.minext[j] + _octree.spans[j] * _geoframe.verts[i][j];
  }

  void Mesher::do_fastcontouring_extract()
  {
    using namespace std;

    cerr << "LBIE::Mesher::extractMesh(): FASTCONTOURING mesh extraction" << endl;

    _geoframe.reset();
    meshType(SINGLE);
    _geoframe.mesh_type = geoframe::GEOTYPE(_meshType);

    _geoframe.setSpan(_octree.spans[0],_octree.spans[1],_octree.spans[2]);
    _geoframe.setMin(_octree.minext[0],_octree.minext[1],_octree.minext[2]);

    FastContouring::TriSurf result = _contourExtractor.extractContour(isovalue());

    //convert the result to a geoframe
    _geoframe.numverts = result.verts.size()/3;
    _geoframe.numtris = result.tris.size()/3;
    
    _geoframe.verts.resize(_geoframe.numverts);
    _geoframe.normals.resize(_geoframe.numverts);
    _geoframe.color.resize(_geoframe.numverts);
    _geoframe.triangles.resize(_geoframe.numtris);

    for(int i=0;i<_geoframe.numverts;i++)
      {
	for(int j=0; j<3; j++)
	  _geoframe.verts[i][j] = result.verts[i*3+j];
	for(int j=0; j<3; j++)
	  _geoframe.normals[i][j] = result.normals[i*3+j];
	for(int j=0; j<3; j++)
	  _geoframe.color[i][j] = result.colors[i*3+j];
      }
    for (int i=0;i<_geoframe.numtris;i++)
      for(int j=0; j<3; j++)
	_geoframe.triangles[i][j] = result.tris[i*3+j];
  }

  void Mesher::do_libisocontour_extract()
  {
    using namespace std;

    cerr << "LBIE::Mesher::extractMesh(): LIBISOCONTOUR mesh extraction" << endl;

    _geoframe.reset();
    meshType(SINGLE);
    _geoframe.mesh_type = geoframe::GEOTYPE(_meshType);

    _geoframe.setSpan(_octree.spans[0],_octree.spans[1],_octree.spans[2]);
    _geoframe.setMin(_octree.minext[0],_octree.minext[1],_octree.minext[2]);

    ConDataset* the_data;
    Contour3dData *contour3d;
    int dim[3];
    float span[3], orig[3];

    //just grab the volume from _contourExtractor
    VolMagick::Volume vol(_contourExtractor.getVolume());

    dim[0] = vol.XDim(); dim[1] = vol.YDim(); dim[2] = vol.ZDim();
    span[0] = vol.XSpan(); span[1] = vol.YSpan(); span[2] = vol.ZSpan();
    orig[0] = vol.XMin(); orig[1] = vol.YMin(); orig[2] = vol.ZMin();
    
    //convert it to a supported libcontour type and load it into libisocontour
    switch(vol.voxelType())
      {
      case VolMagick::UInt:
      case VolMagick::Double:
      case VolMagick::UInt64:
      case VolMagick::Float:
	vol.voxelType(VolMagick::Float);
	the_data = newDatasetReg(CONTOUR_FLOAT, CONTOUR_REG_3D, 1, 1, dim, *vol);
	break;
      case VolMagick::UChar:
	the_data = newDatasetReg(CONTOUR_UCHAR, CONTOUR_REG_3D, 1, 1, dim, *vol);
	break;
      case VolMagick::UShort:
	the_data = newDatasetReg(CONTOUR_USHORT, CONTOUR_REG_3D, 1, 1, dim, *vol);
	break;
      }

    ((Datareg3 *)the_data->data->getData(0))->setOrig(orig);
    ((Datareg3 *)the_data->data->getData(0))->setSpan(span);

    //actually extract the contour
    contour3d = getContour3d(the_data,0,0,isovalue(),NO_COLOR_VARIABLE);

    //convert the result to a geoframe
    _geoframe.numverts = contour3d->nvert;
    _geoframe.numtris = contour3d->ntri;
    
    _geoframe.verts.resize(_geoframe.numverts);
    _geoframe.normals.resize(_geoframe.numverts);
    _geoframe.color.resize(_geoframe.numverts);
    _geoframe.triangles.resize(_geoframe.numtris);

    for(int i=0;i<_geoframe.numverts;i++)
      {
	for(int j=0; j<3; j++)
	  _geoframe.verts[i][j] = contour3d->vert[i][j];
	for(int j=0; j<3; j++)
	  _geoframe.normals[i][j] = contour3d->vnorm[i][j];
	for(int j=0; j<3; j++)
	  _geoframe.color[i][j] = 1.0;
      }
    for (int i=0;i<_geoframe.numtris;i++)
      for(int j=0; j<3; j++)
	_geoframe.triangles[i][j] = contour3d->tri[i][j];

    delete the_data;
    delete contour3d;
  }

  geoframe& Mesher::extractMesh(const VolMagick::Volume& vol)
  {
    setVolume(vol);

    switch(extractionMethod())
      {
      case DUALLIB:
	do_duallib_extract();
	break;
      case FASTCONTOURING:
	do_fastcontouring_extract();
	break;
      case LIBISOCONTOUR:
	do_libisocontour_extract();
	break;
      }

    _geoframe.calculateExtents();
    return _geoframe;
  }

  geoframe& Mesher::extractMesh()
  {
    switch(extractionMethod())
      {
      case DUALLIB:
	do_duallib_extract();
	break;
      case FASTCONTOURING:
	do_fastcontouring_extract();
	break;
      case LIBISOCONTOUR:
	do_libisocontour_extract();
	break;
      }

    _geoframe.calculateExtents();
    return _geoframe;
  }

  geoframe& Mesher::qualityImprove(unsigned int iterations)
  {
    using namespace std;

    //make sure the mesh type is the same as the geoframe
    if(_geoframe.mesh_type != int(_meshType))
      meshType(MeshType(_geoframe.mesh_type));

    for(unsigned int i = 0; i < iterations; i++)
      {
	_octree.quality_improve(_geoframe,_improveMethod);
	cerr << 
	  "LBIE::Mesher::qualityImprove(): "
	  "Quality improvement iteration " << i << " finished." << endl;
      }
    return _geoframe;
  }
}
