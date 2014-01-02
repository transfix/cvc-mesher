#include <iostream>
#include <string>
#include <cstdlib>
#include <VolMagick.h>
#include <LBIE_Mesher.h>
#include <boost/program_options.hpp>
#include <stdexcept>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
  LBIE::Mesher mesher;

  float isovalue, isovalue_in, err, err_in;
  string input_file, output_file;
  string operation, meshtype, improve_method, normaltype, extract_method;
  int improve_iterations;

  try{
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message")
      ("operation", po::value<string>(&operation), "Operation: mesh, quality_improve")
      ("input,i", po::value<string>(&input_file),
       "input file: for meshing, input a volume.  for quality_improve, input a raw geometry file")
      ("output,o", po::value<string>(&output_file), "output raw geometry file")
      ("isovalue", po::value<float>(&isovalue)->default_value(LBIE::DEFAULT_IVAL), "outer iso value")
      ("inner_isovalue", po::value<float>(&isovalue_in)->default_value(LBIE::DEFAULT_IVAL_IN), "inner iso value")
      ("error", po::value<float>(&err)->default_value(LBIE::DEFAULT_ERR), "error threshold")
      ("inner_error", po::value<float>(&err_in)->default_value(LBIE::DEFAULT_ERR_IN), "inner error threshold")
      ("mesh_type", po::value<string>(&meshtype)->default_value("single"), 
       "mesh type to extract: single, tetra, quad, hexa, double, tetra2")
      ("extraction_method", po::value<string>(&extract_method)->default_value("duallib"),
       "mesh extraction method: duallib, fastcontouring, libisocontour")
      ("improvement_method", po::value<string>(&improve_method)->default_value("geo_flow"),
       "mesh quality improvement method: no_improve, geo_flow, edge_contract, joe_liu, minimal_vol, optimization")
      ("iterations", po::value<int>(&improve_iterations)->default_value(1), "quality improvement iterations")
      ("normal_type", po::value<string>(&normaltype)->default_value("bspline_convolution"),
       "mesh normal type: bspline_convolution, central_difference, bspline_interpolation")
      ("dual", "enable dual contouring. Forced for double and tetra2 mesh types")
      ;

    po::positional_options_description p;
    p.add("operation", 1);
    p.add("input", 1);
    p.add("output", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc,argv).
	      options(desc).positional(p).run(), vm);
    po::notify(vm);

    if(vm.count("help"))
      {
	cerr << desc << endl;
	cerr << "Examples:" << endl;
	cerr << argv[0] << " mesh input.rawiv output.raw --isovalue -0.5001" << endl;
	cerr << argv[0] << " quality_improve input.raw output.raw" << endl;
	return EXIT_FAILURE;
      }

    if(!vm.count("operation"))
      {
	cerr << "Error: operation not specified!" << endl;
	cerr << desc << endl;
	return EXIT_FAILURE;
      }

    if(!vm.count("input") || !vm.count("output"))
      {
	cerr << "Error: input and output files must be specified!" << endl;
	cerr << desc << endl;
	return EXIT_FAILURE;
      }

    if(operation == "mesh")
      {
	cout << "isovalue: " << isovalue << endl;
	cout << "inner_isovalue: " << isovalue_in << endl;
	cout << "error: " << err << endl;
	cout << "inner_error: " << err_in << endl;
	cout << "mesh_type: " << meshtype << endl;
	cout << "improvement_method: " << improve_method << endl;
	cout << "iterations: " << improve_iterations << endl;
	cout << "normal_type: " << normaltype << endl;
	cout << "dual: " << (vm.count("dual") ? "true" : "false") << endl;

	mesher.isovalue(isovalue);
	mesher.isovalue_in(isovalue_in);
	mesher.err(err);
	mesher.err_in(err_in);
      
	LBIE::Mesher::MeshType mt;
	if(meshtype == "single") mt = LBIE::Mesher::SINGLE;
	else if(meshtype == "tetra") mt = LBIE::Mesher::TETRA;
	else if(meshtype == "quad") mt = LBIE::Mesher::QUAD;
	else if(meshtype == "hexa") mt = LBIE::Mesher::HEXA;
	else if(meshtype == "double") mt = LBIE::Mesher::DOUBLE;
	else if(meshtype == "tetra2") mt = LBIE::Mesher::TETRA2;
	else
	  {
	    cerr << "Error: invalid mesh type!" << endl;
	    return EXIT_FAILURE;
	  }
	mesher.meshType(mt);

	LBIE::Mesher::ExtractionMethod em;
	if(extract_method == "duallib") em = LBIE::Mesher::DUALLIB;
	else if(extract_method == "fastcontouring") em = LBIE::Mesher::FASTCONTOURING;
	else if(extract_method == "libisocontour") em = LBIE::Mesher::LIBISOCONTOUR;
	else
	  {
	    cerr << "Error: invalid mesh extraction method!" << endl;
	    return EXIT_FAILURE;
	  }
	mesher.extractionMethod(em);
      
	LBIE::Mesher::ImproveMethod im;
	if(improve_method == "no_improve") im = LBIE::Mesher::NO_IMPROVE;
	else if(improve_method == "geo_flow") im = LBIE::Mesher::GEO_FLOW;
	else if(improve_method == "edge_contract") im = LBIE::Mesher::EDGE_CONTRACT;
	else if(improve_method == "joe_liu") im = LBIE::Mesher::JOE_LIU;
	else if(improve_method == "minimal_vol") im = LBIE::Mesher::MINIMAL_VOL;
	else if(improve_method == "optimization") im = LBIE::Mesher::OPTIMIZATION;
	else
	  {
	    cerr << "Error: invalid quality improvement method!" << endl;
	    return EXIT_FAILURE;
	  }
	mesher.improveMethod(im);

	LBIE::Mesher::NormalType nt;
	if(normaltype == "bspline_convolution") nt = LBIE::Mesher::BSPLINE_CONVOLUTION;
	else if(normaltype == "central_difference") nt = LBIE::Mesher::CENTRAL_DIFFERENCE;
	else if(normaltype == "bspline_interpolation") nt = LBIE::Mesher::BSPLINE_INTERPOLATION;
	else
	  {
	    cerr << "Error: invalid normal type!" << endl;
	    return EXIT_FAILURE;
	  }
	mesher.normalType(nt);
      
	mesher.dual(bool(vm.count("dual")));

	VolMagick::Volume vol;
	VolMagick::readVolumeFile(vol,input_file);
      
	mesher.extractMesh(vol); //sets the internal geoframe to the extracted mesh
	mesher.qualityImprove(improve_iterations);

	LBIE::geoframe g_frame = mesher.mesh(); //get the mesh
	g_frame.write_raw(output_file.c_str()); //write it out using geoframe's I/O

	cout << "Wrote mesh: " << output_file << endl;
      }
    else if(operation == "quality_improve")
      {
	cout << "improvement_method: " << improve_method << endl;
	cout << "iterations: " << improve_iterations << endl;

	LBIE::Mesher::ImproveMethod im;
	if(improve_method == "no_improve") im = LBIE::Mesher::NO_IMPROVE;
	else if(improve_method == "geo_flow") im = LBIE::Mesher::GEO_FLOW;
	else if(improve_method == "edge_contract") im = LBIE::Mesher::EDGE_CONTRACT;
	else if(improve_method == "joe_liu") im = LBIE::Mesher::JOE_LIU;
	else if(improve_method == "minimal_vol") im = LBIE::Mesher::MINIMAL_VOL;
	else if(improve_method == "optimization") im = LBIE::Mesher::OPTIMIZATION;
	else
	  {
	    cerr << "Error: invalid quality improvement method!" << endl;
	    return EXIT_FAILURE;
	  }
	mesher.improveMethod(im);
      
	LBIE::geoframe g_frame;
	g_frame.read_raw(input_file.c_str());
	mesher.mesh(g_frame); //set the mesh to be improved
	mesher.qualityImprove(improve_iterations);
      
	g_frame = mesher.mesh(); //get the mesh
	g_frame.write_raw(output_file.c_str()); //write it out using geoframe's I/O

	cout << "Wrote mesh: " << output_file << endl;
      }
    else
      {
	cerr << "Error: unknown operation '" << operation << "'" << endl;
	return EXIT_FAILURE;
      }
  }
  catch(exception &e)
    {
      cerr << "Error: " << e.what() << endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
