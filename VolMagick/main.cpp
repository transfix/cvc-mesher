/* $Id: main.cpp,v 1.2 2008/02/01 20:12:12 transfix Exp $ */

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>

#include "VolMagick.h"
#include "VolumeCache.h"
#include "endians.h"

#include <fstream>

using namespace std;

class VolMagickOpStatus : public VolMagick::VoxelOperationStatusMessenger
{
public:
  void start(const VolMagick::Voxels *vox, Operation op, VolMagick::uint64 numSteps) const
  {
    _numSteps = numSteps;
  }

  void step(const VolMagick::Voxels *vox, Operation op, VolMagick::uint64 curStep) const
  {
    const char *opStrings[] = { "CalculatingMinMax", "CalculatingMin", "CalculatingMax",
				"SubvolumeExtraction", "Fill", "Map", "Resize", "Composite",
				"BilateralFilter", "ContrastEnhancement"};

    fprintf(stderr,"%s: %5.2f %%\r",opStrings[op],(((float)curStep)/((float)((int)(_numSteps-1))))*100.0);
  }

  void end(const VolMagick::Voxels *vox, Operation op) const
  {
    printf("\n");
  }

private:
  mutable VolMagick::uint64 _numSteps;
};

int main(int argc, char **argv)
{
  if(argc < 2)
    {
      cerr << "Usage: " << argv[0] << " <volume file>" << endl;
      return 1;
    }

  try
    {
      VolMagick::VolumeFileInfo volinfo;
      VolMagick::VolumeCache volcache;
      VolMagickOpStatus status;

      VolMagick::setDefaultMessenger(&status);

#if 0
      VolMagick::Volume sphereVol;
      
      sphereVol.dimension(VolMagick::Dimension(128,128,128));//256,256,256));
      sphereVol.voxelType(VolMagick::UChar);

      double center_x = sphereVol.XDim()/2.0;
      double center_y = sphereVol.YDim()/2.0;      
      double center_z = sphereVol.ZDim()/2.0;
      double distance;

      for(unsigned int k=0; k<sphereVol.ZDim(); k++)
	for(unsigned int j=0; j<sphereVol.YDim(); j++)
	  for(unsigned int i=0; i<sphereVol.XDim(); i++)
	    {
	      distance = sqrt(double((i-center_x)*(i-center_x)+
				     (j-center_y)*(j-center_y)+
				     (k-center_z)*(k-center_z)));
	      //sphereVol(i,j,k, distance);
		
	      if((distance > 15.0) && (distance < 20.0))
		sphereVol(i,j,k, 20);//20.0+10*(distance-15.0)/(20.0-15.0));
	      if((distance >= 20.0) && (distance < 25.0))
		sphereVol(i,j,k, 30);//50.0+10*(distance-50.0)/(55.0-50.0));
	    }

      VolMagick::writeVolumeFile(sphereVol, argv[1]);

#endif

#if 0
      volinfo.read(argv[1]);
      cout << volinfo.filename() << ":" <<endl;
      cout << "Num Variables: " << volinfo.numVariables() << endl;
      cout << "Num Timesteps: " << volinfo.numTimesteps() << endl;
      cout << "Dimension: " << volinfo.XDim() << "x" << volinfo.YDim() << "x" << volinfo.ZDim() << endl;
      cout << "Bounding box: ";
      cout << "(" << volinfo.boundingBox().minx << "," << volinfo.boundingBox().miny << "," << volinfo.boundingBox().minz << ") ";
      cout << "(" << volinfo.boundingBox().maxx << "," << volinfo.boundingBox().maxy << "," << volinfo.boundingBox().maxz << ") ";
      cout << endl;
      double volmin = volinfo.min(), volmax = volinfo.max();
      cout << "Min voxel value: " << volmin << endl;
      cout << "Max voxel value: " << volmax << endl;
      cout << "Voxel type: " << volinfo.voxelTypeStr() << endl;
      cout << "Volume name: " << volinfo.name() << endl << endl;
#endif

      /*
      vol.dimension(VolMagick::Dimension(128,128,128));
      vol.voxelType(VolMagick::Float);
      for(unsigned int i=0; i<128; i++)
	for(unsigned int j=0; j<128; j++)
	  for(unsigned int k=0; k<128; k++)
	    vol(i,j,k,double(i)/128.0);
      */
     
#if 0 
      if(argc > 2)
	{
	  VolMagick::Volume vol;

	  readVolumeFile(vol,argv[1]);
	  vol.min(volinfo.min());
	  vol.max(volinfo.max());
	  //vol.map(0.0,255.0);
	  //vol.voxelType(VolMagick::UChar);
	  vol.resize(VolMagick::Dimension(512,512,512));
	  writeVolumeFile(vol,argv[2]);

	}
#endif

#if 0
      if(argc > 2)
	{
	  double realmin, realmax;
	  VolMagick::Volume vol;
	  std::vector<VolMagick::VolumeFileInfo> volinfos(argc-2);
	  volinfos[0].read(argv[1]);
	  realmin = volinfos[0].min();
	  realmax = volinfos[0].max();
	  for(unsigned int i=0; i<argc-2; i++)
	    {
	      volinfos[i].read(argv[i+1]);
	      if(realmin > volinfos[i].min()) realmin = volinfos[i].min();
	      if(realmax < volinfos[i].max()) realmax = volinfos[i].max();
	    }
	  
	  cout << "Realmin: " << realmin << endl;
	  cout << "Realmax: " << realmax << endl;

	  createVolumeFile(argv[argc-1],
			   volinfo.boundingBox(),
			   volinfo.dimension(),
			   std::vector<VolMagick::VoxelType>(1, VolMagick::UChar),
			   1, argc-2,
			   0.0, double(argc-3));

	  for(unsigned int i=0; i<argc-2; i++)
	    {
	      readVolumeFile(vol,argv[i+1]);
	      //so the mapping is correct across all timesteps, we must set the real min and max across time
	      vol.min(realmin);
	      vol.max(realmax);
	      vol.map(0.0,255.0);
	      vol.voxelType(VolMagick::UChar);
	      writeVolumeFile(vol,argv[argc-1],0,i);
	    }
	}
#endif

      //vector test
#if 0
      std::vector<VolMagick::Volume> vols(4);
      {
	char filename[256];

	cout << "Testing RawV!" << endl;
	vols[0].voxelType(VolMagick::UChar);
	vols[0].desc("red");
	vols[0].fill(130.0);
	vols[1].voxelType(VolMagick::UChar);
	vols[1].desc("green");
	vols[1].fill(20.0);
	vols[2].voxelType(VolMagick::UChar);
	vols[2].desc("blue");
	vols[2].fill(200.0);
	
	vols[3].voxelType(VolMagick::Float);
	vols[3].desc("alpha");
	for(unsigned int k=0; k<vols[3].ZDim(); k++)
	  for(unsigned int j=0; j<vols[3].YDim(); j++)
	    for(unsigned int i=0; i<vols[3].XDim(); i++)
	      vols[3](i,j,k,
		      sqrt(float((i - (vols[3].XDim()/2.0))*(i - (vols[3].XDim()/2.0)) + 
				 (j - (vols[3].YDim()/2.0))*(j - (vols[3].YDim()/2.0)) + 
				 (k - (vols[3].ZDim()/2.0))*(k - (vols[3].ZDim()/2.0)))));
	
	VolMagick::writeVolumeFile(vols,"test.rawv");
      }
#endif
      //cout << "std::streamoff size: " << sizeof(std::streamoff) << endl;

      //ifstream test("flkajlff");
    }
  catch(VolMagick::Exception &e)
    {
      cerr << e.what() << endl;
    }
  catch(std::exception &e)
    {
      cerr << e.what() << endl;
    }

  return 0;
}
