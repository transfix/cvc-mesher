
// $Id: seedchkr2.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef SEED_CHKR2_H
#define SEED_CHKR2_H

#include "range.h"
#include "seedcells.h"
#include "conplot.h"
#include "data.h"

class Datavol;
class Dataslc;

class seedChkr2 {
   public:
      seedChkr2(Data &d, SeedCells &s, Conplot &p) : data(d), seeds(s), plot(p) {}
      ~seedChkr2() {}

      void compSeeds(void);

   private:

      Data &data;
      SeedCells &seeds;
      Conplot   &plot;
};

#endif
