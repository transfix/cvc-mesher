
// $Id: rangesweep.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef RANGE_SWEEP_H
#define RANGE_SWEEP_H

#include "range.h"
// #include "spqueue.h"
#include "ipqueue.h"
#include "seedcells.h"
#include "conplot.h"
#include "data.h"

class RangeSweepRec {
   public:
//      operator <(RangeSweepRec&r2)  { return(this->cellid < r2.cellid); }
//      operator >(RangeSweepRec&r2)  { return(this->cellid > r2.cellid); }
//      operator ==(RangeSweepRec&r2) { return(this->cellid == r2.cellid); }

      int cellid;
      Range range;
};

class rangeSweep {
   public:
      rangeSweep(Data &d, SeedCells &s, Conplot &p) : data(d), seeds(s), plot(p) {}
      ~rangeSweep() {}

      void compSeeds(void);

   private:
      void PropagateRegion(int cellid, float min, float max);

//      SortedPriorityQueue<RangeSweepRec> queue;
      IndexedPriorityQueue<RangeSweepRec, double, int> queue;
      Data &data;
      SeedCells &seeds;
      Conplot   &plot;
};

#endif
