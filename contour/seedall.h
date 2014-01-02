
// $Id: seedall.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef SEED_ALL_H
#define SEED_ALL_H

#include "range.h"
#include "seedcells.h"
#include "conplot.h"
#include "data.h"

class seedAll {
   public:
      seedAll(Data &d, SeedCells &s, Conplot &p) : data(d), seeds(s), plot(p) {}
      ~seedAll() {}

      void compSeeds(void);

   private:

      Data &data;
      SeedCells &seeds;
      Conplot   &plot;
};

#endif
