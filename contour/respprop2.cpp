//------------------------------------------------------------------------
//
// respProp2.C - preprocessing of 2d volumes for seed set extraction
//
//------------------------------------------------------------------------

// $Id: respprop2.cpp,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#include <stdlib.h>
#ifndef MACOS_X
#include <malloc.h>
#endif
#include <memory.h>
#ifndef WIN32
#include <unistd.h>
#endif
#include "respprop2.h"
#include "datareg2.h"

#define DEBUGNo

extern int verbose;

//------------------------------------------------------------------------
//
//
//------------------------------------------------------------------------
void
respProp2::compSeeds(void)
{
   Datareg2 &reg2 = (Datareg2&)data;
   int i, j;
   int xdim, ydim;
   float val[4];
   Range prop, c_prop, done, resp, out;
   float min_x, min_y, max_x, max_y;
   int nseed;

   if (verbose)
      printf("***** Seed Creation\n");

   xdim = reg2.dim[0];
   ydim = reg2.dim[1];

   // proceed through the slices computing seeds
   nseed=0;

   // process the k'th slab
   for (i=0; i<xdim-1; i++)
      for (j=0; j<ydim-1; j++) {
         // load the voxel data
         reg2.getCellValues(i, j, val);

         min_x = MIN2(val[0], val[3]);
         max_x = MAX2(val[0], val[3]);
         min_y = MIN2(val[0], val[1]);
         max_y = MAX2(val[0], val[1]);

         // set the incoming values if on a border
         if (j==0) {
            prop.Set(min_y, max_y);
            c_prop.MakeEmpty();
         }

         if (i==0) {
            done.MakeEmpty();
            resp = Range(min_x, max_x);
         }
         else {
            done = Range(min_x, max_x);
            resp.MakeEmpty();
         }

         done += c_prop;

         resp = (prop + Range(MIN2(val[1],val[2]),MAX2(val[1],val[2]))) - done;

         if (j < ydim-2)
            out = Range(MIN2(val[2], val[3]), MAX2(val[2], val[3]));
         else
            out.MakeEmpty();

         // test for propagation of entire responsibility range
         if (resp.Empty() || (!out.Empty() &&
               out.MinAll() <= resp.MinAll() &&
               out.MaxAll() >= resp.MaxAll())) {

              prop = out - done;
              c_prop = out - prop;

           }
           else {
#ifdef DEBUG
printf("ADD SEED\n");
#endif
              // can't propagate all responsiblity, cell must be a seed
              seeds.AddSeed(reg2.index2cell(i,j), resp.MinAll(), resp.MaxAll());

#ifdef DEBUG
printf("seed at %d %d %d\n", i, j, k);
#endif

              nseed++;

              prop.MakeEmpty();
              c_prop = out;
           }
#ifdef DEBUG
sleep(1);
#endif
        }

   if (verbose)
      printf("computed %d seeds\n", nseed);
}
