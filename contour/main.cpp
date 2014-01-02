#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include "contour.h"
#include "datasetreg3.h"
//#include "cellqueue.h"

int main(){

  //char*file[]={{"pot-2eti-glucose-3fields.raw"}};
  char *file="/work/prok/data/rawiv/engine.rawiv";
  int fd;
  ConDataset* the_data; //ConDataset* the_data;
  u_char *buf, *start;
  struct stat st;
  //  the_data = loadDataset(CONTOUR_UCHAR, CONTOUR_REG_3D, 1, 1, file);
  int dim[3] = { 256,256,110 };
  
  fd = open(file,O_RDONLY);
  if(fd == -1)
    {
      fprintf(stderr,"Error opening file '%s'!\n",file);
      return -1;
    }
  fstat(fd, &st);
  start = buf = (u_char*)mmap(NULL,st.st_size,PROT_READ,MAP_PRIVATE,fd,0);
  if(buf == NULL)
    {
      fprintf(stderr,"mmap() failed!\n");
      return -1;
    }
  start+=68; /* skip the header */
  the_data = newDatasetReg(CONTOUR_UCHAR, CONTOUR_REG_3D, 1, 1, dim,start);

  float span[3] = {1.0,1.0,1.0}, orig[3] = {0.0,0.0,0.0};
  ((Datareg3 *)the_data->data->getData(0))->setOrig(orig);
  ((Datareg3 *)the_data->data->getData(0))->setSpan(span);

  int   isovar   = 0;
  int   timestep = 0;
  float isovalue = 91.636f;
  Contour3dData* isocontour;
  int i;

  /*for(i=0;i<4; i++)
    {
      isovalue += 1.0;
      isocontour = getContour3d(the_data,
				isovar   ,
				timestep ,
				isovalue, NO_COLOR_VARIABLE );
      printf("Number of vertices: %d\n", isocontour->nvert);
      delete isocontour;
    }
  for(i=0;i<4; i++)
    {
      isovalue -= 1.0;
      isocontour = getContour3d(the_data,
				isovar   ,
				timestep ,
				isovalue, NO_COLOR_VARIABLE );
      
      printf("Number of vertices: %d\n", isocontour->nvert);
      delete isocontour;
      }*/
  /*for (i=0; i < isocontour->nvert; i++)
    {
      printf("%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	     isocontour->vert[i][0],
	     isocontour->vert[i][1],
	     isocontour->vert[i][2],
	     isocontour->vnorm[i][0],
	     isocontour->vnorm[i][1],
	     isocontour->vnorm[i][2]);
    }

  printf("%d\n", isocontour->ntri);
  for (i=0; i < isocontour->ntri; i++)
    {
      printf("%d %d %d\n", 
	     isocontour->tri[i][0],
	     isocontour->tri[i][1],
	     isocontour->tri[i][2] );
    }
  */
  saveContour3d(the_data,0,0,isovalue,NO_COLOR_VARIABLE,"test.raw"); 
  
  munmap(buf,st.st_size);

  //delete isocontour;
  delete the_data;
  return 0;
}
