//------------------------------------------------------------------------
//
// seedAll.C - preprocessing of 2d volumes for seed set extraction
//
//------------------------------------------------------------------------

// $Id: seedall.cpp,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#include <stdlib.h>
#ifndef MACOS_X
#include <malloc.h>
#endif
#include <memory.h>
#ifndef WIN32
#include <unistd.h>
#endif

#include "seedall.h"
#include "datareg2.h"

#define DEBUGNo

extern int verbose;

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------
void
seedAll::compSeeds(void)
{
   u_int c;
   float min, max;
   int nseed;

   if (verbose)
      printf("***** Seed Creation\n");

   // proceed through the slices computing seeds
   nseed=0;

   for (c=0; c<data.getNCells(); c++) {

         // load the voxel data
         data.getCellRange(c, min, max);

         seeds.AddSeed(c, min, max);

         nseed++;
   }

   if (verbose)
      printf("computed %d seeds\n", nseed);
}
