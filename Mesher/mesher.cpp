#include <VolMagick.h>
#include <Exceptions.h>
#include <mesher.h>

#include <iostream>
#include <string>

// Change Log
// ----------
// 01/11/2014 - Joe R. - separating meshing and quality improve operations into their own functions
//                       to make them easy to call from outside.

namespace LBIE
{
  VOLMAGICK_DEF_EXCEPTION(cvc_mesher_exception);

  // -------
  // do_mesh
  // -------
  // Purpose: 
  //   Main entry point into LBIE meshing.
  // ---- Change History ----
  // 01/11/2014 -- Joe R. -- Creation.
  geoframe do_mesh(const VolMagick::Volume& vol,
		   float isovalue, float isovalue_in, float err, float err_in,
		   const std::string& meshtype, const std::string& improve_method, const std::string& normaltype, 
		   const std::string& extract_method, int improve_iterations, bool dual_contouring,
		   bool verbose)
  {
    using namespace std;

    Mesher mesher;

    if(verbose)
      {
	cout << "isovalue: " << isovalue << endl;
	cout << "inner_isovalue: " << isovalue_in << endl;
	cout << "error: " << err << endl;
	cout << "inner_error: " << err_in << endl;
	cout << "mesh_type: " << meshtype << endl;
	cout << "improvement_method: " << improve_method << endl;
	cout << "iterations: " << improve_iterations << endl;
	cout << "normal_type: " << normaltype << endl;
	cout << "dual: " << (dual_contouring ? "true" : "false") << endl;
      }

    mesher.isovalue(isovalue);
    mesher.isovalue_in(isovalue_in);
    mesher.err(err);
    mesher.err_in(err_in);
      
    Mesher::MeshType mt;
    if(meshtype == "single") mt = Mesher::SINGLE;
    else if(meshtype == "tetra") mt = Mesher::TETRA;
    else if(meshtype == "quad") mt = Mesher::QUAD;
    else if(meshtype == "hexa") mt = Mesher::HEXA;
    else if(meshtype == "double") mt = Mesher::DOUBLE;
    else if(meshtype == "tetra2") mt = Mesher::TETRA2;
    else throw cvc_mesher_exception("invalid mesh type '" + meshtype + "'");
    mesher.meshType(mt);

    Mesher::ExtractionMethod em;
    if(extract_method == "duallib") em = Mesher::DUALLIB;
    else if(extract_method == "fastcontouring") em = Mesher::FASTCONTOURING;
    else if(extract_method == "libisocontour") em = Mesher::LIBISOCONTOUR;
    else throw cvc_mesher_exception("invalid mesh extraction method '" + extract_method + "'");
    mesher.extractionMethod(em);
      
    Mesher::ImproveMethod im;
    if(improve_method == "no_improve") im = Mesher::NO_IMPROVE;
    else if(improve_method == "geo_flow") im = Mesher::GEO_FLOW;
    else if(improve_method == "edge_contract") im = Mesher::EDGE_CONTRACT;
    else if(improve_method == "joe_liu") im = Mesher::JOE_LIU;
    else if(improve_method == "minimal_vol") im = Mesher::MINIMAL_VOL;
    else if(improve_method == "optimization") im = Mesher::OPTIMIZATION;
    else throw cvc_mesher_exception("invalid quality improvement method '" + improve_method + "'");
    mesher.improveMethod(im);
  
    Mesher::NormalType nt;
    if(normaltype == "bspline_convolution") nt = Mesher::BSPLINE_CONVOLUTION;
    else if(normaltype == "central_difference") nt = Mesher::CENTRAL_DIFFERENCE;
    else if(normaltype == "bspline_interpolation") nt = Mesher::BSPLINE_INTERPOLATION;
    else throw cvc_mesher_exception("invalid normal type '" + normaltype + "'");
    mesher.normalType(nt);
    mesher.dual(dual_contouring);
    mesher.extractMesh(vol); //sets the internal geoframe to the extracted mesh
    mesher.qualityImprove(improve_iterations);
    return mesher.mesh();
  }

  // ---------------
  // quality_improve
  // ---------------
  // Purpose: 
  //   Main entry point into LBIE mesh quality improvement
  // ---- Change History ----
  // 01/11/2014 -- Joe R. -- Creation.
  geoframe quality_improve(const geoframe& g_frame, const std::string& improve_method, int improve_iterations,
			   bool verbose)
  {
    using namespace std;
    LBIE::Mesher mesher;
    
    if(verbose)
      {
	cout << "improvement_method: " << improve_method << endl;
	cout << "iterations: " << improve_iterations << endl;
      }

    LBIE::Mesher::ImproveMethod im;
    if(improve_method == "no_improve") im = LBIE::Mesher::NO_IMPROVE;
    else if(improve_method == "geo_flow") im = LBIE::Mesher::GEO_FLOW;
    else if(improve_method == "edge_contract") im = LBIE::Mesher::EDGE_CONTRACT;
    else if(improve_method == "joe_liu") im = LBIE::Mesher::JOE_LIU;
    else if(improve_method == "minimal_vol") im = LBIE::Mesher::MINIMAL_VOL;
    else if(improve_method == "optimization") im = LBIE::Mesher::OPTIMIZATION;
    else cvc_mesher_exception("Error: invalid quality improvement method '" + improve_method + "'");
    mesher.improveMethod(im);
    mesher.mesh(g_frame);
    mesher.qualityImprove(improve_iterations);
    return mesher.mesh(); //get the mesh
  }
}
