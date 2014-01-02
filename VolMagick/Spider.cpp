/* $Id: Spider.cpp,v 1.3 2008/08/18 17:00:10 transfix Exp $ */
/* http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ImportingExporting  */
/* Spider read/write routines ---------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>

#include <boost/scoped_array.hpp>
#include <boost/format.hpp>

#include "VolMagick.h"
#include "endians.h"

#warning TODO: modify this so I/O for spider files is out-of-core.

namespace VolMagick
{
  const char *read_err_strings[] = { "OK",
				     "Cannot open file",
				     "File is not a Spider volume",
				     "Problem when computing the size of the file",
				     "Filetype is neither volume ('V') or image ('I')" };

  typedef struct
  {
    /* 1 */
    float fNslice; // NUMBER OF SLICES (PLANES) IN VOLUME
    // (=1 FOR AN IMAGE)  FOR NEW LONG LABEL
    // FORMAT THE VALUE OF NSLICE STORED IN
    // THE FILE IS NEGATIVE.

    /* 2 */
    float fNrow; // NUMBER OF ROWS PER SLICE (Y)

    /* 3 */
    float fNrec; // TOTAL NUMBER OF RECORDS (SEE NOTE #3).

    /* 4 */
    float fNlabel; // AUXILIARY NUMBER TO COMPUTE TOTAL NUMBER OF RECS

    /* 5 */
    float fIform; // FILE TYPE SPECIFIER.
    // +3 FOR A 3-D FILE  (FLOAT)
    // +1 FOR A 2-D IMAGE (FLOAT)
    // -1 FOR A 2-D FOURIER TRANSFORM
    // -3 FOR A 3-D FOURIER TRANSFORM
    // -5 FOR A NEW 2-D FOURIER TRANSFORM
    // -7 FOR A NEW 3-D FOURIER TRANSFORM
    // +8 FOR A 2-D EIGHT BIT IMAGE FILE
    // +9 FOR A 2-D INT IMAGE FILE
    // 10 FOR A 3-D INT IMAGE FILE
    // 11 FOR A 2-D EIGHT BIT COLOR IMAGE FILE

    /* 6 */
    float fImami; // MAXIMUM/MINIMUM FLAG. IS SET AT 0 WHEN THE
    // FILE IS CREATED, AND AT 1 WHEN THE MAXIMUM AND
    // MINIMUM HAVE BEEN COMPUTED, AND HAVE BEEN STORED
    // INTO THIS LABEL RECORD (SEE FOLLOWING WORDS)

    /* 7 */
    float fFmax; // MAXIMUM VALUE

    /* 8 */
    float fFmin; // MINIMUM VALUE

    /* 9 */
    float fAv; // AVERAGE VALUE

    /* 10*/
    float fSig; // STANDARD DEVIATION. A VALUE OF -1. INDICATES
    // THAT SIG HAS NOT BEEN COMPUTED PREVIOUSLY.

    /* 11*/
    float fIhist; // FLAG INDICATING IF THE HISTOGRAM HAS BE
    // COMPUTED. NOT USED IN 3D FILES!

    /* 12*/
    float fNcol; // NUMBER OF PIXELS PER LINE (Columns X)

    /* 13*/
    float fLabrec; // NUMBER OF LABEL RECORDS IN FILE HEADER

    /* 14*/
    float fIangle; // FLAG THAT TILT ANGLES HAVE BEEN FILLED

    /* 15*/
    float fPhi; // EULER: ROTATIONAL ANGLE

    /* 16*/
    float fTheta; // EULER: TILT ANGLE

    /* 17*/
    float fPsi; // EULER: PSI  = TILT ANGLE

    /* 18*/
    float fXoff; // X TRANSLATION

    /* 19*/
    float fYoff; // Y TRANSLATION

    /* 20*/
    float fZoff; // Z TRANSLATION

    /* 21*/
    float fScale; // SCALE

    /* 22*/
    float fLabbyt; // TOTAL NUMBER OF BYTES IN LABEL

    /* 23*/
    float fLenbyt; // RECORD LENGTH IN BYTES
    char  fNada[24]; // this is a spider incongruence

    /* 30*/
    float fFlag; // THAT ANGLES ARE SET. 1 = ONE ADDITIONAL
    // ROTATION IS PRESENT, 2 = ADDITIONAL ROTATION
    // THAT PRECEEDS THE ROTATION THAT WAS STORED IN
    // 15 FOR DETAILS SEE MANUAL CHAPTER VOCEUL.MAN

    /* 31*/
    float fPhi1;

    /* 32*/
    float fTheta1;

    /* 33*/
    float fPsi1;

    /* 34*/
    float fPhi2;

    /* 35*/
    float fTheta2;

    /* 36*/
    float fPsi2;

    double fGeo_matrix[3][3]; // x9 = 72 bytes: Geometric info
    float fAngle1; // angle info

    float fr1;
    float fr2; // lift up cosine mask parameters

    /** Fraga 23/05/97  For Radon transforms **/
    float RTflag; // 1=RT, 2=FFT(RT)
    float Astart;
    float Aend;
    float Ainc;
    float Rsigma; // 4*7 = 28 bytes
    float Tstart;
    float Tend;
    float Tinc; // 4*3 = 12, 12+28 = 40B

    /** Sjors Scheres 17/12/04 **/
    float Weight; // For Maximum-Likelihood refinement
    float Flip; // 0=no flipping operation (false), 1=flipping (true)

    char fNada2[576]; // empty 700-76-40=624-40-8= 576 bytes

    /*212-214*/
    char szIDat[12]; // LOGICAL * 1 ARRAY DIMENSIONED 10, CONTAINING
    // THE DATE OF CREATION (10 CHARS)

    /*215-216*/
    char szITim[8]; // LOGICAL * 1 ARRAY DIMENSIONED 8, CONTAINING
    // THE TIME OF CREATION (8 CHARS)

    /*217-256*/
    char szITit[160]; // LOGICAL * 1 ARRAY DIMENSIONED 160
  } SpiderHeader;

  size_t FREAD(void *dest, size_t size, size_t nitems, FILE *&fp, int reverse)
  {
#if 0
    size_t retval;
    if (!reverse)
      retval = fread(dest, size, nitems, fp);
    else
      {
	char *ptr = (char *)dest;
	int end = 0;
	retval = 0;
	for (unsigned int n = 0; n < nitems; n++)
	  {
	    for (unsigned int i = size - 1; i >= 0; i--)
	      {
		if (fread(ptr + i, 1, 1, fp) != 1)
		  {
		    end = 1;
		    break;
		  }
	      }
	    if (end)
	      break;
	    else
	      retval++;
	    ptr += size;
	  }
      }
    return retval;
#endif

    size_t retval;
    retval = fread(dest, size, nitems, fp);
    if(reverse)
      {
	char *ptr = (char *)dest;
	switch(size)
	  {
	  case 2:
	    for(unsigned int n = 0; n < nitems; n++)
	      {
		SWAP_16(ptr);
		ptr += 2;
	      }
	    break;
	  case 4:
	    for(unsigned int n = 0; n < nitems; n++)
	      {
		SWAP_32(ptr);
		ptr += 4;
	      }
	    break;
	  case 8:
	    for(unsigned int n = 0; n < nitems; n++)
	      {
		SWAP_64(ptr);
		ptr += 8;
	      }
	    break;
	  default: break; //ignore weird sizes for now
	  }
      }

    return retval;
  }

  /* This function reads a Spider volume. Call it as

  int dim[3];
  float *data=NULL;
  readSpiderFile("myfile.vol",'V',dim,&data);
  readSpiderFile("myfile.img",'I',dim,&data);
    
  Code errors:
  0 - OK
  1 - Cannot open file
  2 - File is not a Spider volume
  3 - Problem when computing the size of the file
  4 - filetype is not 'V' or 'I'
  */
  int readSpiderFile(const char *filename, char filetype,
		     int dim[], float **data)
  {
    FILE* fp=NULL;
    union {
      unsigned char c[4];
      float         f;
    } file_type;
    int fileReversed=0, machineReversed=0, reversed=0;
    SpiderHeader header;
    unsigned long usfNcol, usfNrow, usfNslice, usfHeader, size, tmpSize,
      volSize, n, currentPosition;
    struct stat info;
    float *ptr;

    /* Set dimensions to 0, just in case it cannot be read */
    dim[0]=0;
    dim[1]=0;
    dim[2]=0;    

    /* Check that the filetype is correct */
    if (filetype!='V' and filetype!='I')
      return 4;

    /* Open file */
    if ((fp = fopen(filename, "rb")) == NULL)
      return 1;

    /* Check that the input file is really a Spider volume */
    currentPosition=ftell(fp);

    /* Check file type */
    fseek(fp, currentPosition+16, SEEK_SET);
    for (int i = 0; i < 4; i++)
      fread(&(file_type.c[i]), sizeof(unsigned char), 1, fp);
    fseek(fp,  currentPosition+0, SEEK_SET);
    switch (filetype)
      {
      case 'V':
	if (file_type.c[0]==  0 && file_type.c[1]== 0 &&
	    file_type.c[2]== 64 && file_type.c[3]==64)
	  {
	    fileReversed=0;
	  }
	else if (file_type.c[0]==64 && file_type.c[1]==64 &&
		 file_type.c[2]== 0 && file_type.c[3]== 0)
	  {
	    fileReversed=1;
	  }
	else
	  {
	    fclose(fp);
	    return 2;
	  }
	break;
      case 'I':
	if (file_type.c[0]==  0 && file_type.c[1]== 0 &&
	    file_type.c[2]==128 && file_type.c[3]==63)
	  {
	    fileReversed=0;
	  }
	else if (file_type.c[0]==63 && file_type.c[1]==128 &&
		 file_type.c[2]== 0 && file_type.c[3]==  0)
	  {
	    fileReversed=1;
	  }
	else
	  {
	    fclose(fp);
	    return 2;
	  }
	break;
      }

    /* Check machine type */
    file_type.f = 1;
    if (file_type.c[0]==63 && file_type.c[1]==128 && file_type.c[2]==0 &&
	file_type.c[3]==0)
      machineReversed=1;

    /* Read header */
    reversed=fileReversed ^ machineReversed;
    if (!reversed)
      fread(&header, sizeof(SpiderHeader), 1, fp);
    else
      {
	FREAD(&header,             sizeof(float),  36, fp, true);
	FREAD(&header.fGeo_matrix, sizeof(double),  9, fp, true);
	FREAD(&header.fAngle1,     sizeof(float),  13, fp, true);
	FREAD(&header.fNada2,      sizeof(char),  756, fp, true);
      }

    /* Compute file size, header size and volume dimensions */
    usfNcol = (unsigned long) header.fNcol;
    usfNrow = (unsigned long) header.fNrow;
    usfNslice = (unsigned long) header.fNslice;
    usfHeader = (unsigned long) header.fNcol *
      (unsigned long) header.fLabrec * sizeof(float);
    if (fstat(fileno(fp), &info))
      {
	fclose(fp);
	return 3;
      }

    /* Readjust the number of rows in "aberrant" images*/
    if (filetype=='I' || header.fIform == 1)
      if (usfNcol*usfNrow*sizeof(float) == info.st_size)
	{
	  usfNrow = (unsigned long)(--header.fNrow);
	  --header.fNrec;
	}

    /* Check that the file size is correct */
    switch (filetype)
      {
      case 'I':
	size = usfHeader + usfNcol * usfNrow * sizeof(float);
	if ((size != info.st_size) || (header.fIform != 1))
	  {
	    fclose(fp);
	    return 2;
	  }
	break;
      case 'V':
	size = usfHeader + usfNslice * usfNcol * usfNrow * sizeof(float);
	if ((size != info.st_size) || (header.fIform != 3))
	  {
	    fclose(fp);
	    return 2;
	  }
	break;
      }
    
    /* Read the extra filling header space */
    tmpSize = (unsigned long)(header.fNcol * header.fLabrec * 4);
    tmpSize -= sizeof(SpiderHeader);
    for (unsigned int i = 0; i < tmpSize / 4; i++)
      {
	float tmp;
	fread(&tmp, sizeof(float), 1, fp);
      }
    currentPosition=ftell(fp);

    /* Fill the dimensions */
    dim[0]=(int)header.fNcol;
    dim[1]=(int)header.fNrow;
    if (filetype=='V') dim[2]=(int)header.fNslice;
    else               dim[2]=1;
    volSize = (unsigned long) dim[0]*dim[1]*dim[2];

    /* Read the whole file */
    if (*data!=NULL)
      free(data);
    *data=(float *) malloc(volSize*sizeof(float));
    ptr=*data;
    for (n=0; n<volSize; n++)
      FREAD(ptr++, sizeof(float), 1, fp, reversed);

    /* Close file */
    fclose(fp);
    return 0;
  }

  /* This function writes a Spider volume. Call it as

  writeSpiderFile("myfile.vol",'V',dim,data);
  writeSpiderFile("myfile.img",'I',dim,data);
    
  Code errors:
  0 - OK
  1 - Cannot open file
  */
  int writeSpiderFile(const char *filename, char filetype,
		      int dim[], float *data)
  {
    FILE* fp;
    SpiderHeader header;
    unsigned long volSize, tmpSize;
    float tmp=0.0F;
    int headrec;

    /* Check that the filetype is correct */
    if (filetype!='V' and filetype!='I')
      return 4;

    /* Open file for writing */
    if ((fp = fopen(filename, "wb")) == NULL)
      return 1;

    /* Adjust header */
    header.fNcol = dim[0];
    header.fNrow = dim[1];
    header.fNslice=dim[2];
    if ((header.fNcol != 0) && (header.fNrow != 0))
      {
	header.fNlabel = (float)((int)(256 / header.fNcol + 1));
	header.fLabrec = (float) ceil((float)(256 / (float)header.fNcol));

	headrec = (int) 1024 / ((int)header.fNcol * 4);

	if ((1024 % (int)header.fNcol != 0))
	  {
	    header.fNrec = header.fNrow + 1;
	    headrec++;
	  }
	else
	  header.fNrec = header.fNrow;

	header.fLabbyt = header.fNcol * header.fLabrec * 4;
	header.fLenbyt = (float) header.fNcol * 4;
      }
    if (filetype=='V') header.fIform = 3;
    else if (filetype=='I') header.fIform = 1;
    header.fScale = 1;
    for (unsigned i = 0; i < 3; i++)
      for (unsigned j = 0; j < 3; j++)
	if (i == j)
	  header.fGeo_matrix[i][j] = (double)1.0;
	else
	  header.fGeo_matrix[i][j] = (double)0.0;
    header.fSig = -1;
    header.fImami = header.fFmax = header.fFmin = header.fAv =
      header.fIhist = header.fIangle = header.fPhi = header.fTheta =
      header.fPsi = header.fXoff = header.fYoff = header.fZoff =
      header.fPhi1 = header.fTheta1 = header.fPsi1 =
      header.fPhi2 = header.fTheta2 = header.fPsi2 = header.fAngle1 = 
      header.fr1 = header.fr2 = header.RTflag = header.Astart =
      header.Aend = header.Ainc = header.Rsigma = header.Tstart =
      header.Tend = header.Tinc = header.Weight = header.Flip = 0.0F;
    header.szIDat[0] = header.szITim[0] = header.szITit[0] =
      header.fNada2[0] = header.fNada[0] = 0;

    // Write the header and the necessary empty space
    fwrite(&header, sizeof(SpiderHeader), 1, fp);
    tmpSize = (header.fNcol * header.fLabrec * 4);
    tmpSize -= sizeof(SpiderHeader);
    for (unsigned int i = 0; i < tmpSize/4; i++)
      fwrite(&tmp, sizeof(float), 1, fp);

    // Now write the data
    fwrite(data, sizeof(float), (unsigned long) dim[0]*dim[1]*dim[2], fp);

    // Close file
    fclose(fp);
    return 0;
  }

  void VolumeFileInfo::readSpider(const std::string& file)
  {
    int dim[3];
    float *data=NULL;
    int retval = readSpiderFile(file.c_str(),'I',dim,&data);
    if(retval)
      retval = readSpiderFile(file.c_str(),'V',dim,&data);
    if(retval)
      {
	if(data) free(data);
	throw InvalidSpiderFile(read_err_strings[retval]);
      }
    
    /**** at this point, the volume header has no errors, so we may start modifiying this object ****/
    _numVariables = 1;
    _numTimesteps = 1;
    _tmin = _tmax = 0.0;
    
    _dimension = Dimension(dim);
    _boundingBox = BoundingBox(0.0, 0.0, 0.0,
			       double(dim[0]-1),double(dim[1]-1),double(dim[2]-1));

    _voxelTypes.clear();
    _voxelTypes.push_back(VolMagick::Float);
    _names.clear();
    _names.push_back("No Name");

    /* new volume, so min/max is now unset */
    _minIsSet.clear();
    _minIsSet.resize(_numVariables); for(int i=0; i<_minIsSet.size(); i++) _minIsSet[i].resize(_numTimesteps);
    _min.clear();
    _min.resize(_numVariables); for(int i=0; i<_min.size(); i++) _min[i].resize(_numTimesteps);
    _maxIsSet.clear();
    _maxIsSet.resize(_numVariables); for(int i=0; i<_maxIsSet.size(); i++) _maxIsSet[i].resize(_numTimesteps);
    _max.clear();
    _max.resize(_numVariables); for(int i=0; i<_max.size(); i++) _max[i].resize(_numTimesteps);
  }

  void readSpider(Volume& vol,
		  const std::string& filename, 
		  unsigned int var, unsigned int time,
		  uint64 off_x, uint64 off_y, uint64 off_z,
		  const Dimension& subvoldim)
  {
    int dim[3];
    float *data=NULL;

    if(var > 0)
      throw IndexOutOfBounds("Variable index out of bounds.");
    if(time > 0)
      throw IndexOutOfBounds("Timestep index out of bounds.");
    if(subvoldim.isNull())
      throw IndexOutOfBounds("Specified subvolume dimension is null.");

    int retval = readSpiderFile(filename.c_str(),'I',dim,&data);
    if(retval)
      retval = readSpiderFile(filename.c_str(),'V',dim,&data);
    if(retval)
      {
	if(data) free(data);
	throw InvalidSpiderFile(read_err_strings[retval]);
      }

    if((off_x + subvoldim[0] - 1 >= (unsigned int)dim[0]) ||
       (off_y + subvoldim[1] - 1 >= (unsigned int)dim[1]) ||
       (off_z + subvoldim[2] - 1 >= (unsigned int)dim[2]))
      throw IndexOutOfBounds("Subvolume specified is outside volume dimensions");

    vol.voxelType(Float);
    vol.dimension(subvoldim);
    vol.boundingBox(BoundingBox(double(off_x),
				double(off_y),
				double(off_z),
				double(off_x+subvoldim[0]-1),
				double(off_y+subvoldim[1]-1),
				double(off_z+subvoldim[2]-1)));

    for(uint64 k = 0; k < vol.ZDim(); k++)
      for(uint64 j = 0; j < vol.YDim(); j++)
	for(uint64 i = 0; i < vol.XDim(); i++)
	  vol(i,j,k, data[dim[0]*dim[1]*(k+off_z) + dim[0]*(j+off_y) + (i+off_x)]);
  }

  void createSpider(const std::string& filename,
		    const BoundingBox& boundingBox,
		    const Dimension& dimension,
		    const std::vector<VoxelType>& voxelTypes,
		    unsigned int numVariables, unsigned int numTimesteps,
		    double min_time, double max_time)
  {
    using namespace boost;
    int dim[3] = { dimension[0], dimension[1], dimension[2] };
    boost::scoped_array<float> data(new float[dimension.size()]);

    if(boundingBox.isNull())
      throw InvalidBoundingBox("Bounding box must not be null");
    if(dimension.isNull())
      throw InvalidBoundingBox("Dimension must not be null");
    if(numVariables > 1)
      throw InvalidSpiderFile(str(format("Spider format only supports 1 variable (%d requested)") % numVariables));
    if(numTimesteps > 1)
      throw InvalidSpiderFile(str(format("Spider format only supports 1 timestep (%d requested)") % numTimesteps));
    if(voxelTypes.size() > 1)
      throw InvalidSpiderFile("Spider format only supports 1 variable and 1 timestep. (too many voxel types specified)");
    if(voxelTypes[0] != Float)
      throw InvalidSpiderFile("Spider format only supports float valued voxels.");
    if(min_time != max_time)
      throw InvalidSpiderFile("Spider format does not support multiple timesteps. (min time and max time must be the same)");

    int retval = writeSpiderFile(filename.c_str(), 'V',dim,data.get());
    if(retval)
      throw InvalidSpiderFile(read_err_strings[retval]);
  }

  void writeSpider(const Volume& wvol, 
		   const std::string& filename,
		   unsigned int var, unsigned int time,
		   uint64 off_x, uint64 off_y, uint64 off_z)
  {
    Volume vol;

    if(var > 0)
      throw IndexOutOfBounds("Variable index out of bounds.");
    if(time > 0)
      throw IndexOutOfBounds("Timestep index out of bounds.");

    try
      {
	VolumeFileInfo volinfo(filename);
	readSpider(vol,filename,0,0,0,0,0,volinfo.dimension());
      }
    catch(ReadError &e)
      {
	createSpider(filename,
		     BoundingBox(0.0,0.0,0.0,1.0,1.0,1.0),
		     wvol.dimension(),
		     std::vector<VoxelType>(wvol.voxelType()),
		     0,0,0.0,0.0);
	readSpider(vol,filename,0,0,0,0,0,wvol.dimension());
      }

    for(uint64 k = 0; k < wvol.ZDim(); k++)
      for(uint64 j = 0; j < wvol.YDim(); j++)
	for(uint64 i = 0; i < wvol.XDim(); i++)
	  vol(i+off_x,j+off_y,k+off_z, wvol(i,j,k));

    int dim[] = { vol.XDim(), vol.YDim(), vol.ZDim() };
    vol.voxelType(Float); //make sure we have a float volume (though we shouldnt get here if we dont)
    int retval = writeSpiderFile(filename.c_str(),'V',dim,reinterpret_cast<float*>(*vol));
    if(retval)
      throw InvalidSpiderFile(read_err_strings[retval]);
  }

}
