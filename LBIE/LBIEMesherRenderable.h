#ifndef __LBIEMESHERRENDERABLE_H__
#define __LBIEMESHERRENDERABLE_H__

//similar to MultiContour class
class LBIEMesherRenderable : public Renderable
{
 public:
  LBIEMesherRenderable();
  virtual ~LBIEMesherRenderable();

  //single variable data load
  void setData(unsigned char* data, 
	       unsigned int width, unsigned int height, unsigned int depth,
	       double aspectX, double aspectY, double aspectZ,
	       double subMinX, double subMinY, double subMinZ,
	       double subMaxX, double subMaxY, double subMaxZ,
	       double minX, double minY, double minZ,
	       double maxX, double maxY, double maxZ);

  //RGBA data load
  void setData(unsigned char* data, unsigned char* red, unsigned char* green, unsigned char* blue,
		unsigned int width, unsigned int height, unsigned int depth,
		double aspectX, double aspectY, double aspectZ,
		double subMinX, double subMinY, double subMinZ,
		double subMaxX, double subMaxY, double subMaxZ,
		double minX, double minY, double minZ,
		double maxX, double maxY, double maxZ);

  
};

#endif
