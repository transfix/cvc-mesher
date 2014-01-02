//------------------------------------------------------------
//
// seedCells.C - maintain the list of seed cells
//
// Copyright (c) 1997 Dan Schikore
//------------------------------------------------------------

// $Id: seedcells.cpp,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef MACOS_X
#include <malloc.h>
#endif
#include <stdlib.h>
#include "seedcells.h"

//------------------------------------------------------------
//
// SeedCells() - initialize the list
//
//------------------------------------------------------------
SeedCells::SeedCells()
{
   ncells = 0;
   cell_size = 10000;
   cells = (SeedCellP)malloc(sizeof(struct SeedCell) * cell_size);
}

//------------------------------------------------------------
//
// ~SeedCells() - free storage
//
//------------------------------------------------------------
SeedCells::~SeedCells()
{
   free(cells);
}

//------------------------------------------------------------
//
// AddSeed() - add a seed cell, increasing storage as necessary
//
//------------------------------------------------------------
int
SeedCells::AddSeed(u_int id, float min, float max)
{
   int n = ncells++;

   if (n >= cell_size) {
      cell_size *= 2;
      cells = (SeedCellP)realloc(cells, sizeof(struct SeedCell) * cell_size);
   }

   cells[n].cell_id = id;
   cells[n].min = min;
   cells[n].max = max;

   return(n);
}
