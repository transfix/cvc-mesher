//------------------------------------------------------------------------
//
// conPlot3d.C - preprocess and extract contours from 3d scalar data
//
// Copyright (c) 1997 Dan Schikore - modified by Emilio Camahort, 1999
//
//------------------------------------------------------------------------

// $Id: conplot3d.cpp,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

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
#include "contour3d.h"
#include "range.h"
#include "segtree.h"
#include "conplot3d.h"

extern int verbose;
extern void (*errorHandler)(char *, int);

//------------------------------------------------------------------------
//
// unst3d, tetadj - static local variables
//
//------------------------------------------------------------------------

static int unst3d[16][7] = {
   {0},
   {1, 0, 2, 3},
   {1, 0, 4, 1},
   {2, 2, 3, 1, 3, 4, 1},
   {1, 2, 1, 5},
   {2, 3, 1, 5, 3, 0, 1},
   {2, 2, 0, 5, 0, 4, 5},
   {1, 3, 4, 5},
   {1, 3, 5, 4},
   {2, 2, 5, 0, 0, 5, 4},
   {2, 3, 5, 1, 3, 1, 0},
   {1, 2, 5, 1},
   {2, 2, 1, 3, 3, 1, 4},
   {1, 0, 1, 4},
   {1, 0, 3, 2},
   {0}
};

static int tetadj[16][5] = {
   {0},
   {3, 1, 2, 3},
   {3, 0, 2, 3},
   {4, 0, 1, 2, 3},
   {3, 0, 1, 3},
   {4, 0, 1, 2, 3},
   {4, 0, 1, 2, 3},
   {3, 0, 1, 2},
   {3, 0, 1, 2},
   {4, 0, 1, 2, 3},
   {4, 0, 1, 2, 3},
   {3, 0, 1, 3},
   {4, 0, 1, 2, 3},
   {3, 0, 2, 3},
   {3, 1, 2, 3},
   {0},
};

//------------------------------------------------------------------------
//
// Conplot() - create a contour plot for the given volume.
//
//------------------------------------------------------------------------

Conplot3d::Conplot3d(Datasetvol *d) : Conplot(d)
{
   float min[3], max[3];
   int i;

   vol = d;

if (verbose > 1) {
   printf("***** Data Characteristics\n");
   printf("cells: %d\n", vol->getNCells());
   printf("*****\n");
}

   contour2d = NULL;
   contour3d = con3 = new Contour3d[vol->nTime()];
   data->getData(0)->getExtent(min, max);

if (verbose) {
printf("minextent: %f %f %f\n", min[0], min[1], min[2]);
printf("maxextent: %f %f %f\n", max[0], max[1], max[2]);
}

   for (i=0; i<vol->nTime(); i++)
      con3[i].setExtent(min,max);

if (verbose > 1)
   printf("contour2d is %p, contour3d is %p\n", contour2d, con3);
}

//------------------------------------------------------------------------
//
// ~Conplot3d() - destroy a plot
//
//------------------------------------------------------------------------

Conplot3d::~Conplot3d()
{
}

//------------------------------------------------------------------------
//
// InterpEdge() -
//
//------------------------------------------------------------------------

int
Conplot3d::InterpEdge(int edge, float *val, u_int *v, float isovalue, int cell)
{
   float ival;
   float pt[3];
   float norm[3];

   switch (edge) {
      case 0:
         ival = (isovalue-val[1])/(val[0]-val[1]);
         pt[0] = curvol->getVert(v[1])[0]*(1.0f-ival) +
                 curvol->getVert(v[0])[0]*ival;
         pt[1] = curvol->getVert(v[1])[1]*(1.0f-ival) +
                 curvol->getVert(v[0])[1]*ival;
         pt[2] = curvol->getVert(v[1])[2]*(1.0f-ival) +
                 curvol->getVert(v[0])[2]*ival;
         norm[0] = curvol->getGrad(v[1])[0]*(1.0f-ival) +
                 curvol->getGrad(v[0])[0]*ival;
         norm[1] = curvol->getGrad(v[1])[1]*(1.0f-ival) +
                 curvol->getGrad(v[0])[1]*ival;
         norm[2] = curvol->getGrad(v[1])[2]*(1.0f-ival) +
                 curvol->getGrad(v[0])[2]*ival;
         break;
      case 1:
         ival = (isovalue-val[2])/(val[1]-val[2]);
         pt[0] = curvol->getVert(v[2])[0]*(1.0f-ival) +
                 curvol->getVert(v[1])[0]*ival;
         pt[1] = curvol->getVert(v[2])[1]*(1.0f-ival) +
                 curvol->getVert(v[1])[1]*ival;
         pt[2] = curvol->getVert(v[2])[2]*(1.0f-ival) +
                 curvol->getVert(v[1])[2]*ival;
         norm[0] = curvol->getGrad(v[2])[0]*(1.0f-ival) +
                 curvol->getGrad(v[1])[0]*ival;
         norm[1] = curvol->getGrad(v[2])[1]*(1.0f-ival) +
                 curvol->getGrad(v[1])[1]*ival;
         norm[2] = curvol->getGrad(v[2])[2]*(1.0f-ival) +
                 curvol->getGrad(v[1])[2]*ival;
         break;
      case 2:
         ival = (isovalue-val[0])/(val[2]-val[0]);
         pt[0] = curvol->getVert(v[0])[0]*(1.0f-ival) +
                 curvol->getVert(v[2])[0]*ival;
         pt[1] = curvol->getVert(v[0])[1]*(1.0f-ival) +
                 curvol->getVert(v[2])[1]*ival;
         pt[2] = curvol->getVert(v[0])[2]*(1.0f-ival) +
                 curvol->getVert(v[2])[2]*ival;
         norm[0] = curvol->getGrad(v[0])[0]*(1.0f-ival) +
                 curvol->getGrad(v[2])[0]*ival;
         norm[1] = curvol->getGrad(v[0])[1]*(1.0f-ival) +
                 curvol->getGrad(v[2])[1]*ival;
         norm[2] = curvol->getGrad(v[0])[2]*(1.0f-ival) +
                 curvol->getGrad(v[2])[2]*ival;
         break;
      case 3:
         ival = (isovalue-val[0])/(val[3]-val[0]);
         pt[0] = curvol->getVert(v[0])[0]*(1.0f-ival) +
                 curvol->getVert(v[3])[0]*ival;
         pt[1] = curvol->getVert(v[0])[1]*(1.0f-ival) +
                 curvol->getVert(v[3])[1]*ival;
         pt[2] = curvol->getVert(v[0])[2]*(1.0f-ival) +
                 curvol->getVert(v[3])[2]*ival;
         norm[0] = curvol->getGrad(v[0])[0]*(1.0f-ival) +
                 curvol->getGrad(v[3])[0]*ival;
         norm[1] = curvol->getGrad(v[0])[1]*(1.0f-ival) +
                 curvol->getGrad(v[3])[1]*ival;
         norm[2] = curvol->getGrad(v[0])[2]*(1.0f-ival) +
                 curvol->getGrad(v[3])[2]*ival;
         break;
      case 4:
         ival = (isovalue-val[1])/(val[3]-val[1]);
         pt[0] = curvol->getVert(v[1])[0]*(1.0f-ival) +
                 curvol->getVert(v[3])[0]*ival;
         pt[1] = curvol->getVert(v[1])[1]*(1.0f-ival) +
                 curvol->getVert(v[3])[1]*ival;
         pt[2] = curvol->getVert(v[1])[2]*(1.0f-ival) +
                 curvol->getVert(v[3])[2]*ival;
         norm[0] = curvol->getGrad(v[1])[0]*(1.0f-ival) +
                 curvol->getGrad(v[3])[0]*ival;
         norm[1] = curvol->getGrad(v[1])[1]*(1.0f-ival) +
                 curvol->getGrad(v[3])[1]*ival;
         norm[2] = curvol->getGrad(v[1])[2]*(1.0f-ival) +
                 curvol->getGrad(v[3])[2]*ival;
         break;
      case 5:
         ival = (isovalue-val[2])/(val[3]-val[2]);
         pt[0] = curvol->getVert(v[2])[0]*(1.0f-ival) +
                 curvol->getVert(v[3])[0]*ival;
         pt[1] = curvol->getVert(v[2])[1]*(1.0f-ival) +
                 curvol->getVert(v[3])[1]*ival;
         pt[2] = curvol->getVert(v[2])[2]*(1.0f-ival) +
                 curvol->getVert(v[3])[2]*ival;
         norm[0] = curvol->getGrad(v[2])[0]*(1.0f-ival) +
                 curvol->getGrad(v[3])[0]*ival;
         norm[1] = curvol->getGrad(v[2])[1]*(1.0f-ival) +
                 curvol->getGrad(v[3])[1]*ival;
         norm[2] = curvol->getGrad(v[2])[2]*(1.0f-ival) +
                 curvol->getGrad(v[3])[2]*ival;
         break;
   }

//printf("pt %f %f\n", pt[0], pt[1]);

   {
      float len = (float)sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
      if (len != 0.0) {
         norm[0]/=len;
         norm[1]/=len;
         norm[2]/=len;
      }
   }

   return(curcon->AddVert(pt, norm));
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

inline void Conplot3d::EnqueueFaces(int code, int c, CellQueue &thisslice)
{
   switch (code) {
      case 0:
         break;
      case 1:
         break;
      case 2:
         break;
   }
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

void Conplot3d::TrackContour(float isovalue, int cell)
{
   float val[4];
   u_int *verts;
   u_int v1, v2, v3;
   int code;
   int adj;
   int e, a;
   int nvert=0, ntri=0;			// to save isocontour components

   queue.Add(cell);

   curvol = (Datavol*)data->getData(curtime);
   curcon = &con3[curtime];

   if (filePrefix)			// keep track of current nvert, ntri
      {
      nvert = curcon->getNVert();
      ntri  = curcon->getNTri();
      }

   while (queue.Get(cell) > 0)
      {
      curvol->getCellValues(cell, val);
      verts = curvol->getCellVerts(cell);

      code = 0;
      if (val[0] < isovalue) code += 0x01;
      if (val[1] < isovalue) code += 0x02;
      if (val[2] < isovalue) code += 0x04;
      if (val[3] < isovalue) code += 0x08;

      for (e=0; e<unst3d[code][0]; e++)
	 {
         v1 = InterpEdge(unst3d[code][3*e+1], val, verts, isovalue, cell);
         v2 = InterpEdge(unst3d[code][3*e+2], val, verts, isovalue, cell);
         v3 = InterpEdge(unst3d[code][3*e+3], val, verts, isovalue, cell);
         curcon->AddTri(v1, v2, v3);

         for (a=0; a<tetadj[code][0]; a++)
	    {
            adj = curvol->getCellAdj(cell, tetadj[code][a+1]);
            if (adj != -1 && !CellTouched(adj))
	        {
                TouchCell(adj);
                queue.Add(adj);
		}
            }
         }
      }
   
    if (filePrefix)			// write this isocontour component
	{
	if (curcon->getNTri() - ntri > 25)	// more than 25 triangles
	    {
	    FILE	*fp;
	    int		v, t;
	    char	filename[200];

	    sprintf(filename, "%s%04d.ipoly", filePrefix, ncomponents);
	    if ( (fp = fopen(filename, "w")) )
		{
		fprintf(fp, "%d 0 %d 0 0 0 0\n0 0 0\n",
		    curcon->getNVert() - nvert, curcon->getNTri() - ntri);

		for (v = nvert; v < curcon->getNVert(); v++)
		    fprintf(fp, "%g %g %g\n", curcon->vert[v][0],
			    curcon->vert[v][1], curcon->vert[v][2]);

		fprintf(fp, "0 0\n");

		for (t = ntri; t < curcon->getNTri(); t++)
		    fprintf(fp, "3\n%d %d %d\n", curcon->tri[t][0],
			    curcon->tri[t][1], curcon->tri[t][2]);

		fclose(fp);
		ncomponents++;
		}
	    else
                {
                char    str[256];
   
                sprintf(str, "Conplot3d::TrackContour: couldn't open file: %s",
                        filename);
                errorHandler(str, FALSE);
                }
	    }
	}
}
