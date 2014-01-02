//------------------------------------------------------------------------
//
// conplotreg2.C - preprocess and extract contours from 2d scalar data
//
// Copyright (c) 1997 Dan Schikore
//------------------------------------------------------------------------

// $Id: conplotreg2.cpp,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#include <stdlib.h>
#ifndef MACOS_X
#include <malloc.h>
#endif
#include <memory.h>
#include <string.h>
#ifndef WIN32
#include <unistd.h>
#endif

#include "conplot.h"
#include "contour2d.h"
#include "range.h"
#include "segtree.h"
#include "conplotreg2.h"

extern int verbose;
extern void (*errorHandler)(char *, int);

//------------------------------------------------------------------------
//
// rect2d - a local static variable
//
//------------------------------------------------------------------------

static int rect2d[16][5] = {
   {0},
   {1, 3, 0},
   {1, 0, 1},
   {1, 3, 1},
   {1, 1, 2},
   {2, 3, 2, 1, 0},
   {1, 0, 2},
   {1, 3, 2},
   {1, 2, 3},
   {1, 2, 0},
   {2, 0, 3, 2, 1},
   {1, 2, 1},
   {1, 1, 3},
   {1, 1, 0},
   {1, 0, 3},
   {0}
};

//------------------------------------------------------------------------
//
// Conplotreg2() - create a contour plot for the given image.
//
//------------------------------------------------------------------------

Conplotreg2::Conplotreg2(Datasetreg2 *d) : Conplot(d)
{
   float min[3], max[3];
   int i;

   reg2 = d;

if (verbose > 1) {
   printf("***** Data Characteristics\n");
   printf("cells: %d\n", reg2->getNCells());
   printf("*****\n");
}

   contour2d = con2 = new Contour2d[reg2->nTime()];
   contour3d = NULL;
   data->getData(0)->getExtent(min, max);

if (verbose) {
printf("minextent: %f %f %f\n", min[0], min[1], min[2]);
printf("maxextent: %f %f %f\n", max[0], max[1], max[2]);
}
   for (i=0; i<reg2->nTime(); i++)
      con2[i].setExtent(min,max);

if (verbose > 1)
   printf("contour3d is %p, contour2d is %p\n", contour3d, con2);
}

//------------------------------------------------------------------------
//
// ~Conplotreg2() - destroy a plot
//
//------------------------------------------------------------------------

Conplotreg2::~Conplotreg2()
{
}

//------------------------------------------------------------------------
//
// InterpEdge() - interpolate ....
//
//------------------------------------------------------------------------

int
Conplotreg2::InterpEdge(int edge, float *val, float isovalue, int i, int j)
{
   float ival;
   float pt[2];

   switch (edge) {
      case 0:
         ival = (isovalue-val[1])/(val[0]-val[1]);
         pt[0] = (1.0f-ival) * curreg2->xCoord(i+1) +
                       ival  * curreg2->xCoord(i+0);
         pt[1] = curreg2->yCoord(j);
         break;
      case 1:
         ival = (isovalue-val[2])/(val[1]-val[2]);
         pt[0] = curreg2->xCoord(i+1);
         pt[1] = (1.0f-ival) * curreg2->yCoord(j+1) +
                       ival  * curreg2->yCoord(j+0);
         break;
      case 2:
         ival = (isovalue-val[3])/(val[2]-val[3]);
         pt[0] = (1.0f-ival) * curreg2->xCoord(i+0) +
                       ival  * curreg2->xCoord(i+1);
         pt[1] = curreg2->yCoord(j+1);
         break;
      case 3:
         ival = (isovalue-val[0])/(val[3]-val[0]);
         pt[0] = curreg2->xCoord(i);
         pt[1] = (1.0f-ival) * curreg2->yCoord(j+0) +
                       ival  * curreg2->yCoord(j+1);
         break;
   }

   return(curcon->AddVert(pt));
}

//------------------------------------------------------------------------
//
// EnqueueFaces() - enqueue adjacent faces for propagation of contour
//             code      = case table lookup code for current cell
//             i,j,k     = index of current cell
//             thisslice = queue of cells to be processed on current slice
//             nextslice = queue of cells to be processed on next slice
//             trackinz  = flag indicating whether to propagate in z-dir
//                         0 -> don't track in z, but queue cells in nextslice
//                         1 -> propagate contour in all 3 dim
//
//------------------------------------------------------------------------

inline void Conplotreg2::EnqueueFaces(int code, int c, CellQueue &thisslice)
{
}

//------------------------------------------------------------------------
//
// TrackContour() - compute and track a contour by table lookup and propagation
//                  through adjacent cells
//            isovalue  = surface value of interest
//            i,j,k     = index of seed cell
//            nextslice = queue of seeds for next slice (used only if tracking in 2d)
//
//------------------------------------------------------------------------

void Conplotreg2::TrackContour(float isovalue, int cell)
{
   float val[4];
   int i, j;
   u_int v1, v2;
   int code;
   int adj;
   int e;
   int nvert=0, nedge=0;			// to save isocontour components

   queue.Add(cell);

   curreg2 = (Datareg2*)data->getData(curtime);
   curcon = &con2[curtime];

   if (filePrefix)			// keep track of current nvert, nedge
      {
      nvert = curcon->getNVert();
      nedge = curcon->getNEdge();
      }

   while (queue.Get(cell) > 0)
      {
      curreg2->getCellValues(cell, val);
      curreg2->cell2index(cell, i, j);

      code = 0;
      if (val[0] < isovalue) code += 0x01;
      if (val[1] < isovalue) code += 0x02;
      if (val[2] < isovalue) code += 0x04;
      if (val[3] < isovalue) code += 0x08;

      for (e=0; e<rect2d[code][0]; e++)
	 {
         v1 = InterpEdge(rect2d[code][2*e+1], val, isovalue, i, j);
         v2 = InterpEdge(rect2d[code][2*e+2], val, isovalue, i, j);
         curcon->AddEdge(v1, v2);

         adj = curreg2->getCellAdj(cell, rect2d[code][2*e+1]);
         if (adj != -1 && !CellTouched(adj))
	    {
            TouchCell(adj);
            queue.Add(adj);
            }
         adj = curreg2->getCellAdj(cell, rect2d[code][2*e+2]);
         if (adj != -1 && !CellTouched(adj))
	    {
            TouchCell(adj);
            queue.Add(adj);
            }
         }

      EnqueueFaces(code, cell, queue);
      }

    if (filePrefix)			// write this isocontour component
	{
      	if (curcon->getNEdge() - nedge > 25)		// more than 25 edges
	    {
	    FILE	*fp;
	    int		v, e;
	    char	filename[200];

	    sprintf(filename, "%s%04d.ipoly", filePrefix, ncomponents);
	    if ( (fp = fopen(filename, "w")) )
		{
		fprintf(fp, "%d %d 0 0 0 0 0\n0 0 0\n", 
		    curcon->getNVert() - nvert, curcon->getNEdge() - nedge);

		for (v = nvert; v < curcon->getNVert(); v++)
		    fprintf(fp, "%g %g %g\n", curcon->vert[v][0],
					      curcon->vert[v][1], 0.0);
		
		fprintf(fp, "0 0\n");

		for (e = nedge; e < curcon->getNEdge(); e++)
		    fprintf(fp, "%d %d\n", curcon->edge[e][0], 
					   curcon->edge[e][1]);

		fclose(fp);
		ncomponents++;
		}
	    else
		{
                char    str[256];

                sprintf(str,
			"Conplotreg2::TrackContour: couldn't open file: %s",
                        filename);
                errorHandler(str, FALSE);
		}
	    }
	}
}
