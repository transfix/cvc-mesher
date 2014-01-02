//----------------------------------------------------------------------------
//
// bucketSearch.h - segment tree data structure
//
//----------------------------------------------------------------------------

// $Id: bucketsearch.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef BUCKET_SEARCH_H
#define BUCKET_SEARCH_H

#include <sys/types.h>

#include "cellsearch.h"

//----------------------------------------------------------------------------
//
// Bucket search structure
//
//----------------------------------------------------------------------------
class BucketSearch : public CellSearch {
   public:
      BucketSearch(u_int n = 0, float *v = NULL);
      ~BucketSearch();

      void Init(u_int n, float *v);
      void InsertSeg(u_int cellid, float min, float max);
      void Dump(void);
      void Info(void);
      void Traverse(float, void (*f)(u_int, void*), void *);
      u_int getCells(float, u_int *);
      void Done(void);

   protected:
      u_int whichBucket(float f) { return u_int(f-minval); }

   private:
      int nbuckets;
      float minval, maxval;
      CellBucket *buckets;
};

#endif
