//----------------------------------------------------------------------------
//
// cellSearch.h - cell search structure
//
//----------------------------------------------------------------------------

// $Id: cellsearch.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef CELL_SEARCH_H
#define CELL_SEARCH_H

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>

#ifdef WIN32
typedef unsigned int	u_int;
#endif

extern int verbose;

//----------------------------------------------------------------------------
//
// list of cells which cross a given segment
//
//----------------------------------------------------------------------------
class CellBucket {
   public:
      CellBucket();
      ~CellBucket();

      void insert(u_int cellid);

      int nCells(void)     { return(ncells); }
      u_int getCell(u_int i) { return(cells[i]); }
      void getCells(u_int *, u_int &);
      void traverseCells(void (*f)(u_int, void *), void*);
      void dump(char *str);
      u_int *getCells(void) { return(cells); }

   private:
      int ncells;
      int cellsize;
      u_int *cells;
};

inline
CellBucket::CellBucket()
{
   ncells = 0;
   cellsize = 0;
   cells = NULL;
}

inline
CellBucket::~CellBucket()
{
   if (cells != NULL)
      free(cells);
}

inline
void
CellBucket::insert(u_int cellid)
{
   int n = ncells++;

   if (n >= cellsize) {
      if (cellsize == 0) {
         cellsize = 5;
         cells    = (u_int *)malloc(sizeof(u_int)*cellsize);
      }
      else {
         cellsize *= 2;
         cells     = (u_int *)realloc(cells, sizeof(u_int)*cellsize);
      }
   }

   cells[n] = cellid;
}

inline
void
CellBucket::getCells(u_int *a, u_int &n)
{
   memcpy(&a[n], cells, sizeof(u_int)*ncells);
   n += ncells;
}

inline
void
CellBucket::traverseCells(void (*f)(u_int, void *), void *data)
{
   int i;

   for (i=0; i<ncells; i++)
      (*f)(cells[i], data);
}

inline
void
CellBucket::dump(char *str)
{
   int i;

   printf(str);
   for (i=0; i<ncells; i++) {
      printf("%d ", cells[i]);
   }
   printf("\n");
}

//----------------------------------------------------------------------------
//
// Abstract class for cell search structure
//
//----------------------------------------------------------------------------
class CellSearch {
   public:
      CellSearch() { if (verbose) printf("cellsearch constructor!!\n");}
      virtual ~CellSearch() { if (verbose) printf("cellsearch destructor\n"); }

      virtual void  Done(void) = 0;
      virtual void  Init(u_int, float *) = 0;
      virtual void  Dump(void) = 0;
      virtual void  Info(void) = 0;
      virtual void  Traverse(float, void (*f)(u_int, void*), void *) = 0;
      virtual u_int getCells(float, u_int *) = 0;
      virtual void  InsertSeg(u_int, float, float) = 0;

   protected:

   private:
};

#endif
