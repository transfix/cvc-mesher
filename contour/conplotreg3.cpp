//------------------------------------------------------------------------
//
// conPlotReg3.C - preprocess and extract contours from 3d scalar data
//
//------------------------------------------------------------------------

// $Id: conplotreg3.cpp,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

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
#include "conplotreg3.h"
#include "vtkmarchingcubescases.h"
#include "cubes.h"

#define UNIQUE

extern int verbose;
extern void (*errorHandler)(char *, int);

//------------------------------------------------------------------------
//
// edgeinfo - table of information about cube edges.
//            dir      = direction of edge (x=0, y=1, z=2)
//            di,dj,dk = the index of the lowermost point in reference
//                       to the cube
//            d1,d2    = the vertices of the edge
//
//------------------------------------------------------------------------

typedef struct {
   int dir;
   int di,dj,dk;
   int d1,d2;
} EdgeInfo;

static EdgeInfo edgeinfo[12] = {
   { 0, 0, 0, 0, 0, 1 },
   { 2, 1, 0, 0, 1, 2 },
   { 0, 0, 0, 1, 3, 2 },
   { 2, 0, 0, 0, 0, 3 },
   { 0, 0, 1, 0, 4, 5 },
   { 2, 1, 1, 0, 5, 6 },
   { 0, 0, 1, 1, 7, 6 },
   { 2, 0, 1, 0, 4, 7 },
   { 1, 0, 0, 0, 0, 4 },
   { 1, 1, 0, 0, 1, 5 },
   { 1, 0, 0, 1, 3, 7 },
   { 1, 1, 0, 1, 2, 6 }
};

//------------------------------------------------------------------------
//
// Conplot() - create a contour plot for the given volume.
//
//------------------------------------------------------------------------

Conplotreg3::Conplotreg3(Datasetreg3 *d) : Conplot(d)
{
   int		i;
   float	min[3], max[3];

   reg3 = d;

if (verbose > 1) {
   printf("***** Data Characteristics\n");
   printf("cells: %d\n", curreg3->getNCells());
   printf("*****\n");
}

   filePrefix = NULL;
   contour2d  = NULL;
   contour3d  = con3 = new Contour3d[reg3->nTime()];
   data->getData(0)->getExtent(min, max);

if (verbose) {
printf("minextent: %f %f %f\n", min[0], min[1], min[2]);
printf("maxextent: %f %f %f\n", max[0], max[1], max[2]);
}
    for (i = 0; i < reg3->nTime(); i++)
	con3[i].setExtent(min, max);

if (verbose > 1)
   printf("contour2d is %p, contour3d is %p\n", contour2d, con3);
}

//------------------------------------------------------------------------
//
// ~Conplotreg3() - destroy a plot
//
//------------------------------------------------------------------------

Conplotreg3::~Conplotreg3()
{
}

//------------------------------------------------------------------------
//
// interpRect3Dpts_x() - interpolate along an edge in the x direction
//
//------------------------------------------------------------------------
void
Conplotreg3::interpRect3Dpts_x(int i1, int j1, int k1, float *data, float *data2,
             float grad[3][8], int d1, int d2, float val, float *pt,
             float *norm, float *fval)
{
   double ival;

   ival = (val - data[d1])/(data[d2] - data[d1]);
   pt[0] = (float)(curreg3->orig[0] + curreg3->span[0]*(i1 + ival));
   pt[1] = (float)(curreg3->orig[1] + curreg3->span[1]*j1);
   pt[2] = (float)(curreg3->orig[2] + curreg3->span[2]*k1);

   norm[0] = (float)(grad[0][d1]*(1.0-ival) + grad[0][d2]*ival);
   norm[1] = (float)(grad[1][d1]*(1.0-ival) + grad[1][d2]*ival);
   norm[2] = (float)(grad[2][d1]*(1.0-ival) + grad[2][d2]*ival);

   *fval = (float)(data2[d1]*(1.0-ival) + data2[d2]*ival);
}


//------------------------------------------------------------------------
//
// interpRect3Dpts_y() - interpolate along an edge in the y direction
//
//------------------------------------------------------------------------
void
Conplotreg3::interpRect3Dpts_y(int i1, int j1, int k1, float *data, float *data2,
             float grad[3][8], int d1, int d2, float val, float *pt,
             float *norm, float *fval)
{
   double ival;

   ival = (val - data[d1])/(data[d2] - data[d1]);
   pt[0] = curreg3->orig[0] + curreg3->span[0]*i1;
   pt[1] = (float)(curreg3->orig[1] + curreg3->span[1]*(j1 + ival));
   pt[2] = curreg3->orig[2] + curreg3->span[2]*k1;

   norm[0] = (float)(grad[0][d1]*(1.0-ival) + grad[0][d2]*ival);
   norm[1] = (float)(grad[1][d1]*(1.0-ival) + grad[1][d2]*ival);
   norm[2] = (float)(grad[2][d1]*(1.0-ival) + grad[2][d2]*ival);

   *fval = (float)(data2[d1]*(1.0-ival) + data2[d2]*ival);
}


//------------------------------------------------------------------------
//
// interpRect3Dpts_z() - interpolate along an edge in the z direction
//
//------------------------------------------------------------------------
void
Conplotreg3::interpRect3Dpts_z(int i1, int j1, int k1, float *data, float *data2,
             float grad[3][8], int d1, int d2, float val, float *pt,
             float *norm, float *fval)
{
   double ival;

   ival = (val - data[d1])/(data[d2] - data[d1]);
   pt[0] = curreg3->orig[0] + curreg3->span[0]*i1;
   pt[1] = curreg3->orig[1] + curreg3->span[1]*j1;
   pt[2] = (float)(curreg3->orig[2] + curreg3->span[2]*(k1+ival));

   norm[0] = (float)(grad[0][d1]*(1.0-ival) + grad[0][d2]*ival);
   norm[1] = (float)(grad[1][d1]*(1.0-ival) + grad[1][d2]*ival);
   norm[2] = (float)(grad[2][d1]*(1.0-ival) + grad[2][d2]*ival);

   *fval = (float)(data2[d1]*(1.0-ival) + data2[d2]*ival);
}

//------------------------------------------------------------------------
//
// InterpEdge() - interpolate .....
//
//------------------------------------------------------------------------
int
Conplotreg3::InterpEdge(float *val, float *val2, float grad[3][8], float isovalue,
                        int i, int j, int k, int edge)
{
   float pt[3];
   float norm[3];
   float fval;
   EdgeInfo *ei = &edgeinfo[edge];
   int v;

					// do the interpolation

   switch (ei->dir) {
      case 0:
         interpRect3Dpts_x(i+ei->di, j+ei->dj, k+ei->dk, val, val2, grad,
                           ei->d1, ei->d2, isovalue, pt, norm, &fval);
         break;
      case 1:
         interpRect3Dpts_y(i+ei->di, j+ei->dj, k+ei->dk, val, val2, grad,
                           ei->d1, ei->d2, isovalue, pt, norm, &fval);
         break;
      case 2:
         interpRect3Dpts_z(i+ei->di, j+ei->dj, k+ei->dk, val, val2, grad,
                           ei->d1, ei->d2, isovalue, pt, norm, &fval);
         break;
   }
					// add to the coordinate list

   {
     double len = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
     norm[0] = (float)(norm[0]/len);
     norm[1] = (float)(norm[1]/len);
     norm[2] = (float)(norm[2]/len);
   }
#ifdef UNIQUE
   v = curcon->AddVertUnique(pt, norm, fval);
#else
   v = curcon->AddVert(pt, norm, fval);
#endif

   return(v);
}

//------------------------------------------------------------------------
//
// EnqueueFaces() - enqueue adjacent faces for propagation of contour
//             code      = case table lookup code for current cell
//             i,j,k     = index of current cell
//             queue = queue of cells to be processed on current slice
//             nextslice = queue of cells to be processed on next slice
//             trackinz  = flag indicating whether to propagate in z-dir
//                         0 -> don't track in z, but queue cells in nextslice
//                         1 -> propagate contour in all 3 dim
//
//------------------------------------------------------------------------
inline void
Conplotreg3::EnqueueFaces(int code, u_int i, u_int j, u_int k, CellQueue &queue)
{
   int f;
   int id;

   for (f=0; f<adjfaces[code][0]; f++) {
      switch (adjfaces[code][1+f]) {
         case 0:
            if (j >= 1) {
               id = curreg3->index2cell(i, j-1, k);
               if (!CellTouched(id)) {
                  TouchCell(id);
                  queue.Add(id);
               }
            }
            break;
         case 1:
            if (i >= 1) {
               id = curreg3->index2cell(i-1, j, k);
               if (!CellTouched(id)) {
                  TouchCell(id);
                  queue.Add(id);
               }
            }
            break;
         case 2:
            if (j <= curreg3->dim[1]-3) {
               id = curreg3->index2cell(i, j+1, k);
               if (!CellTouched(id)) {
                  TouchCell(id);
                  queue.Add(id);
               }
            }
            break;
         case 3:
            if (i <= curreg3->dim[0]-3) {
               id = curreg3->index2cell(i+1, j, k);
               if (!CellTouched(id)) {
                  TouchCell(id);
                  queue.Add(id);
               }
            }
            break;
         case 4:
            if (k <= curreg3->dim[2]-3) {
               id = curreg3->index2cell(i, j, k+1);
               if (!CellTouched(id)) {
                  TouchCell(id);
                  queue.Add(id);
               }
            }
            break;
         case 5:
            if (k >= 1) {
               id = curreg3->index2cell(i, j, k-1);
               if (!CellTouched(id)) {
                  TouchCell(id);
                  queue.Add(id);
               }
            }
            break;
      }
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

void Conplotreg3::TrackContour(float isovalue, int cell)
{
    int edge_v[12];
    float val[8];
    float val2[8];
    float grad[3][8];
    u_int v1, v2, v3;
    int i, j, k;
    int e, t;
    int code;
    int edge;
    int nvert=0, ntri=0;			// to save isocontour components

    queue.Add(cell);

    curreg3 = (Datareg3*)data->getData(curtime);
    curcon  = &con3[curtime];

    if (filePrefix)			// keep track of current nvert, ntri
	{
	nvert = curcon->getNVert();
	ntri  = curcon->getNTri();
	}

    if (curreg3->getColorFun() != -1)
	curcon->minmaxFun(data->getMinFun(curreg3->getColorFun()),
			  data->getMaxFun(curreg3->getColorFun()));

    while (queue.Get(cell) > 0)
	{
					// select contour for contour function

	// printf("track contour at cell:%d\n", cell);
	// by fan        curreg3->getCellValues(cell, val /*,0*/);
	// Emilio:	 curreg3->getCellValues(cell, val, funcontour);

	curreg3->getCellValues(cell, val, curreg3->getContourFun());

					// use color of color function

	if (curreg3->getColorFun() != -1)
	//if (curreg3->getNData() > 1)
	    // by fan   curreg3->getCellValues(cell, val2, 1);
	    // Emilio:	curreg3->getCellValues(cell, val2, funcolor);
	    curreg3->getCellValues(cell, val2, curreg3->getColorFun());

	curreg3->getCellGrad(cell, grad);
	curreg3->cell2index(cell, i, j, k);

	code = 0;
	if (val[0] < isovalue) code |= 0x01;
	if (val[1] < isovalue) code |= 0x02;
	if (val[2] < isovalue) code |= 0x04;
	if (val[3] < isovalue) code |= 0x08;
	if (val[4] < isovalue) code |= 0x10;
	if (val[5] < isovalue) code |= 0x20;
	if (val[6] < isovalue) code |= 0x40;
	if (val[7] < isovalue) code |= 0x80;

	for (e = 0; e < cubeedges[code][0]; e++)
	    {
	    edge = cubeedges[code][1+e];
            edge_v[edge] = InterpEdge(val, val2, grad, isovalue, i, j, k, edge);
	    }

        for (t=0; triCases[code].edges[t] != -1; )
	    {
            v1 = edge_v[triCases[code].edges[t++]];
            v2 = edge_v[triCases[code].edges[t++]];
            v3 = edge_v[triCases[code].edges[t++]];
            curcon->AddTri(v1, v2, v3);
	    }

        EnqueueFaces(code, i, j, k, queue);
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
   
                sprintf(str, "Conplotreg3::TrackContour: couldn't open file: %s",
                        filename);
                errorHandler(str, FALSE);
                }
	    }
	}
}
