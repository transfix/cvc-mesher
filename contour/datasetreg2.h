//--------------------------------------------------------------------
//
// Datasetreg2 - representation for a 2D time-varying regular grid
//
// Copyright (c) 1997 Dan Schikore - updated by Emilio Camahort, 1999 
//
//--------------------------------------------------------------------

// $Id: datasetreg2.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef DATASET_REG2_H
#define DATASET_REG2_H

#include "dataset.h"
#include "datareg2.h"

extern int verbose;

//--------------------------------------------------------------------
//
// Datasetreg2 - a scalar time-varying dataset
//
//--------------------------------------------------------------------

class Datasetreg2 : public Dataset
{
   private:				// data member

      Datareg2	**reg2;

   public:				// constructors and destructors

      Datasetreg2(Data::DataType t, int ndata, int ntime, char *files[]);
      Datasetreg2(Data::DataType t, int ndata, int ntime, int *dim, 
		  u_char *data);
      ~Datasetreg2() {}

      					// member access methods

      float getMin(int t) const { return(reg2[t]->getMin()); }
      float getMax(int t) const { return(reg2[t]->getMax()); }
      float getMin()      const { return(min[0]); }
      float getMax()      const { return(max[0]); }

      // add by fan
      float getMinFun(int j)      const { return(min[j]); }
      float getMaxFun(int j)      const { return(max[j]); }


      Data	*getData(int i) { return(reg2[i]); }
      Datareg2	*getMesh(int i) { return(reg2[i]); }
};

//------------------------------------------------------------------------
//
// Datasetreg2() - usual constructor, reads data from one or more files
//
//------------------------------------------------------------------------

inline Datasetreg2::Datasetreg2(Data::DataType t, int nd, int nt, char *fn[])
		               : Dataset(t, nd, nt, fn)
{
   int i,j;

   meshtype = 4;
   reg2 = (Datareg2 **)malloc(sizeof(Datareg2 *)*nt);
   for (j = 0; j < nd; j++)
       {
       min[j] = 1e10;
       max[j] = -1e10;
       }
   ncells = 0;
   maxcellindex = 0;
   for (i=0; i<nt; i++)
      {
if (verbose)
printf("loading file: %s\n", fn[i]);
      reg2[i] = new Datareg2(t, nd, fn[i]);
      for (j = 0; j < nd; j++)
	  {
	  if (reg2[i]->getMin() < min[j])
	     min[j] = reg2[i]->getMin();
	  if (reg2[i]->getMax() > max[j])
	     max[j] = reg2[i]->getMax();
	  }
      if (reg2[i]->getNCells() > ncells)
         ncells = reg2[i]->getNCells();
      if (reg2[i]->maxCellIndex() > maxcellindex)
         maxcellindex = reg2[i]->maxCellIndex();
      }
}

//------------------------------------------------------------------------
//
// Datasetreg2() - alternative constructor for the libcontour library
//
//------------------------------------------------------------------------

inline Datasetreg2::Datasetreg2(Data::DataType t, int ndata, int ntime,
		    		int *dim, u_char *data)
				: Dataset(t, ndata, ntime, data)
{
    int	i;				// timestep index variable
    int j;				// a variable index
    int	size = 0;		// size of single timestep of data

    meshtype = 4;
    reg2 = (Datareg2 **)malloc(sizeof(Datareg2 *)*ntime);
    for (j = 0; j < ndata; j++)
	{
	min[j] = 1e10;
	max[j] = -1e10;
	}
    ncells = 0;
    maxcellindex = 0;

    switch (t)
	{
	case Data::UCHAR :  size = dim[0] * dim[1] * ndata * sizeof(u_char);
			    break;
	case Data::USHORT : size = dim[0] * dim[1] * ndata * sizeof(u_short);
			    break;
	case Data::FLOAT :  size = dim[0] * dim[1] * ndata * sizeof(float);
			    break;
	}

    for (i=0; i<ntime; i++)
	{
	reg2[i] = new Datareg2(t, ndata, dim, data + i*size);
	for (j = 0; j < ndata; j++)
	    {
	    if (reg2[i]->getMin() < min[j])
		min[j] = reg2[i]->getMin();
	    if (reg2[i]->getMax() > max[j])
		max[j] = reg2[i]->getMax();
	    }
	if (reg2[i]->getNCells() > ncells)
	    ncells = reg2[i]->getNCells();
	if (reg2[i]->maxCellIndex() > maxcellindex)
	    maxcellindex = reg2[i]->maxCellIndex();
      }
}

#endif
