//------------------------------------------------------------------------
//
// conplot3d.h - class for preprocessing and extraction of surfaces from
//             3d data
//
// Copyright (c) 1997 Dan Schikore
//------------------------------------------------------------------------

// $Id: conplot3d.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef CONPLOT_3D_H
#define CONPLOT_3D_H

#include <sys/types.h>
#include "contour3d.h"
#include "dataset.h"
#include "segtree.h"
#include "seedcells.h"
#include "cellqueue.h"
#include "edgehash.h"
#include "range.h"
#include "datasetvol.h"

#include "conplot.h"

//------------------------------------------------------------------------
//
// conplot3d.h
//
//------------------------------------------------------------------------
class Conplot3d : public Conplot {
   public:
      Conplot3d(Datasetvol *d);
      ~Conplot3d();

   protected:
      // extract in 3d (from memory) or slice-by-slice (swap from disk)
      u_int ExtractAll(float isovalue);

      int InterpEdge(int, float *, u_int *, float, int);

      // track a contour from a seed cell
      void TrackContour(float, int);

      // enqueue faces for propagation of surface
      inline void EnqueueFaces(int, int, CellQueue &);

      void Reset(int t)   { con3[t].Reset();           }
      int  Size(int t)    { return(con3[t].getSize()); }
      int  isDone(int t)  { return(con3[t].isDone());  }
      void Done(int t)    { con3[t].Done(); }

   private:
      Datasetvol *vol;
      Datavol *curvol;
      Contour3d *con3, *curcon;
};

#endif
