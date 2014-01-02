//----------------------------------------------------------------
//
// seedCells.h - maintain a list of seed cells
//
// Copyright (c) 1997 Dan Schikore
//----------------------------------------------------------------

// $Id: seedcells.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef SEED_CELLS_H
#define SEED_CELLS_H

#include <sys/types.h>

#ifdef WIN32
typedef unsigned int	u_int;
#endif

typedef struct SeedCell {
   float min, max;
   u_int cell_id;
} *SeedCellP;

class SeedCells {
   public:
     SeedCells();
     ~SeedCells();

     int    getNCells(void)    { return(ncells); }
     u_int  getCellID(int i)   { return(cells[i].cell_id); }
     float  getMin(int i)      { return(cells[i].min); }
     float  getMax(int i)      { return(cells[i].max); }
     void   Clear(void)        { ncells = 0; }
     SeedCell *getCellPointer(){ return(cells); }

     int AddSeed(u_int, float, float);
     void AddToRange(u_int i, float mn, float mx)
          {
             if (mn < cells[i].min)
                cells[i].min = mn;
             if (mx > cells[i].max)
                cells[i].max = mx;
          }

   private:
     int ncells;
     int cell_size;
     SeedCellP cells;
};

#endif
