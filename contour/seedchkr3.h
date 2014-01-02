
// $Id: seedchkr3.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef SEED_CHKR3_H
#define SEED_CHKR3_H

#include "range.h"
#include "seedcells.h"
#include "conplot.h"
#include "data.h"

class Datavol;
class Dataslc;

class seedChkr3 {
   public:
      seedChkr3(Data &d, SeedCells &s, Conplot &p) : data(d), seeds(s), plot(p) {}
      ~seedChkr3() {}

      void compSeeds(void);

   private:

      Data &data;
      SeedCells &seeds;
      Conplot   &plot;
};

#endif
