#include <iostream>
#include <string>
#include <cstdlib>
#include <VolMagick.h>
#include <Exceptions.h>
#include <mesher.h>
#include <boost/program_options.hpp>
#include <stdexcept>

int main(int argc, char **argv)
{
  using namespace std;
  namespace po = boost::program_options;

  float isovalue, isovalue_in, err, err_in;
  string input_file, output_file;
  string operation, meshtype, improve_method, normaltype, extract_method;
  bool dual_contouring;
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

    dual_contouring = bool(vm.count("dual"));

    if(operation == "mesh")
      {
	VolMagick::Volume vol;
	VolMagick::readVolumeFile(vol,input_file);
	LBIE::geoframe g_frame = LBIE::do_mesh(vol,
					       isovalue, isovalue_in, err, err_in,
					       meshtype, improve_method, normaltype, extract_method,
					       improve_iterations, dual_contouring, true);
	g_frame.write_raw(output_file.c_str());
	cout << "Wrote mesh: " << output_file << endl;
      }
    else if(operation == "quality_improve")
      {
	LBIE::geoframe g_frame;
	g_frame.read_raw(input_file.c_str());
	g_frame = LBIE::quality_improve(g_frame, improve_method, improve_iterations, true);
	g_frame.write_raw(output_file.c_str()); //write it out using geoframe's I/O
      }
  }
  catch(exception &e)
    {
      cerr << "Error: " << e.what() << endl;
      return EXIT_FAILURE;
    }
}
