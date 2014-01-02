//------------------------------------------------------------------------
//
// conPlot2d.h - class for preprocessing and extraction of surfaces from
//             3d data
//
// Copyright (c) 1997 Dan Schikore
//------------------------------------------------------------------------

// $Id: conplotreg3.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef CONPLOT_RECT_3D_H
#define CONPLOT_RECT_3D_H

#include <sys/types.h>
#include "contour3d.h"
#include "dataset.h"
#include "segtree.h"
#include "seedcells.h"
#include "cellqueue.h"
#include "edgehash.h"
#include "range.h"
#include "datasetreg3.h"
#include "conplot.h"

//------------------------------------------------------------------------
//
// conPlotreg3.h
//
//------------------------------------------------------------------------
class Conplotreg3 : public Conplot {
   public:
      Conplotreg3(Datasetreg3 *d);
      ~Conplotreg3();

   protected:
      // extract in 3d (from memory) or slice-by-slice (swap from disk)
      u_int ExtractAll(float isovalue);

      void interpRect3Dpts_x(int, int, int, float *, float *, float [3][8],
                       int, int, float, float *, float *, float *);
      void interpRect3Dpts_y(int, int, int, float *, float *, float [3][8],
                       int, int, float, float *, float *, float *);
      void interpRect3Dpts_z(int, int, int, float *, float *, float [3][8],
                       int, int, float, float *, float *, float *);

      int InterpEdge(float *, float *, float [3][8], float, int, int, int, int);

      // track a contour from a seed cell
      void TrackContour(float, int);

      // enqueue faces for propagation of surface
      inline void EnqueueFaces(int, u_int, u_int, u_int, CellQueue &);

      void Reset(int t)  { con3[t].Reset();           }
      int  Size(int t)   { return(con3[t].getSize()); }
      int  isDone(int t) { return(con3[t].isDone());  }
      void Done(int t)   { con3[t].Done(); }

    private:

	Datasetreg3	*reg3;
	Datareg3	*curreg3;

    public :					// modified by Emilio

	Contour3d *con3, *curcon;
};

#endif
