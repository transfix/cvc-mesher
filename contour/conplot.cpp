//------------------------------------------------------------------------
//
// conplot.C - preprocess and extract contours from 3d scalar data
//
// Copyright (c) 1997 Dan Schikore
//------------------------------------------------------------------------

// $Id: conplot.cpp,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#include <stdlib.h>
#ifndef MACOS_X
#include <malloc.h>
#endif
#include <memory.h>
#include <string.h>
#ifndef WIN32
#include <unistd.h>
#include <sys/time.h>
#endif

#include "range.h"
#include "segtree.h"
#include "inttree.h"
#include "bucketsearch.h"

#include "conplot.h"

extern int verbose;

//#define NQUERY 10000
//#define TIME_SEARCH

//#define WRITE

//------------------------------------------------------------------------
//
// Conplot() - create a contour plot for the given volume.
//
//------------------------------------------------------------------------

//int  Conplot::funcontour;
//int  Conplot::funcolor;

Conplot::Conplot(Dataset *d)
{
   data	= d;
   contour2d = NULL;
   contour3d = NULL;
   filePrefix = NULL;

if (verbose) {
   printf("***** Data Characteristics\n");
   printf("cells: %d\n", data->getNCells());
   printf("*****\n");
}

   // initialize the bit array of 'touched' (visited) cells
   touched = (u_char *)malloc(sizeof(u_char) * (data->maxCellIndex()+7)>>3);

   int_cells = (u_int *)malloc(sizeof(u_int) * data->maxCellIndex());
if (verbose)
printf("initializing %d trees\n", data->nTime());
#ifdef USE_SEG_TREE
   tree = new SegTree[data->nTime()];
#elif defined USE_INT_TREE
   tree = new IntTree[data->nTime()];
#elif defined USE_BUCKETS
   tree = new BucketSearch[data->nTime()];
#endif
//   tree = new IntTree[data->nTime()];

   seeds = new SeedCells[data->nTime()];	// initialize seed data array

   curtime = 0;
}


//------------------------------------------------------------------------
//
// ~Conplot() - destroy a plot
//
//------------------------------------------------------------------------
Conplot::~Conplot()
{
}

void
Conplot::setTime(int t)
{
   curtime = t;
}

//------------------------------------------------------------------------
//
// ExtractAll() - extract an isosurface by propagation in 3d.  Data is
//                assumed to reside in memory.
//            isovalue  = surface value of interest
//
//------------------------------------------------------------------------
u_int
Conplot::ExtractAll(float isovalue)
{
    int n;
    int cur;
#ifdef TIME_SEARCH
    time_t start, finish;
    int t;
#endif

    if (isDone(curtime))
	return(Size(curtime));

#ifdef TIME_SEARCH
   start = clock();
   for (t=0; t<NQUERY; t++)
      n=tree[curtime].getCells(isovalue, int_cells);
   finish = clock();
   printf("%f seconds for %d queries\n", (finish-start)/(float)CLOCKS_PER_SEC, NQUERY);
   printf("%f seconds/query\n", (finish-start)/((float)(CLOCKS_PER_SEC)*NQUERY));
#endif

					// find the intersected seeds

    n = tree[curtime].getCells(isovalue, int_cells);
if (verbose)
    printf("%d intersected seeds\n", n);


					// flush the old surface
    Reset(curtime);
					// clear bit array of 'touched' cells
    ClearTouched();
//  memset(touched, 0, sizeof(u_char) * data->getNCells()>>3);

					// loop through the seeds in order
    for (cur = 0; cur < n; cur++)
	{
	if (!CellTouched(int_cells[cur]))
	    {
	    TouchCell(int_cells[cur]);
	    TrackContour(isovalue, int_cells[cur]);
	    }
	}

if (verbose)
    if (contour3d) printf("%d triangles\n", contour3d->getNTri());

    Done(curtime);

#ifdef WRITE
    if (contour3d) contour3d->write("output.tmesh");
#endif

    return(Size(curtime));
}
