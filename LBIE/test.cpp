#include <iostream>
#include <stdio.h>

#include <math.h>

#include "octree.h"
#include "LBIE_Mesher.h"

using namespace LBIE;

Octree oc;

int main( int argc, char** argv )
{	
	// We use head65.rawiv as a generic test input file
//	const char* inputFile = "head65.rawiv";
	const char* inputFile = "/Workspace/transfix/head65.rawiv";
	
	// Files for storing the various different outputs 

	const char* output_single = "test_single.raw";
	const char* output_hexa = "test_hexa.raw";
	const char* output_double = "test_double.raw";
	const char* output_tetra = "test_tetra.raw";
	const char* output_t_4_h = "test_t_4_h.raw";
	const char* output_tetra2 = "test_tetra2.raw";

	
	/****
	Arguements for library tester
	LIBE_Mesher( input filename, 
				 output filename,
				 outer isovalue, 
				 inner isovalue, 
				 outer err tol, 
				 inner err tol, 
				 meshtype
				) 
	To test the library, we first create an output .raw file with some arguements.
	Then we compare this with the .raw file created by LBIE-Mesher program.
	If they match, then the library is similar to the program.
	****/


//	LBIE_Mesher( inputFile, output_single, 127.0, 177.0, 0.2, 0.0, 0 );
//	std::cout<< "Done!\n\n";

	LBIE_Mesher( inputFile, output_hexa,   0.5, 0.0, 0.2, 0.0, 1 );
	std::cout<< "Done!\n\n";
#if 0	
	LBIE_Mesher( inputFile, output_tetra,  0.5, 0.0, 0.2, 0.0, 3 );
	std::cout<< "Done!\n\n";
	
	LBIE_Mesher( inputFile, output_t_4_h,  0.5, 0.0, 0.2, 0.0, 4 );
	std::cout<< "Done!\n\n";
#endif
	/*
	LBIE_Mesher( inputFile, output_double, 0.5, -0.1,0.2, -0.3,2 );
	std::cout<< "Done!\n\n";
	*/
/*	
	LBIE_Mesher( inputFile, output_tetra2, 0.5, -0.1,0.2, -0.3,5 );
	std::cout<< "Done!\n\n";
	*/
	
	return 0;

}
