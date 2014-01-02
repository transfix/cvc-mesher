//--------------------------------------------------------------------
//
// Datasetvol - representation for a time-varying volume
//
// Copyright (c) 1997 Dan Schikore - updated by Emilio Camahort, 1999
//
//--------------------------------------------------------------------

// $Id: datasetvol.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef DATASET_VOL_H
#define DATASET_VOL_H

#include "data.h"
#include "dataset.h"
#include "datavol.h"

//--------------------------------------------------------------------
//
// Datasetvol - a scalar time-varying dataset
//
//--------------------------------------------------------------------

class Datasetvol : public Dataset
{
   private:				// data members

      Datavol	**vol;

   public:				// constructors and destructors

      Datasetvol(Data::DataType t, int ndata, int ntime, char *files[]);
      Datasetvol(Data::DataType t, int ndata, int ntime, u_int nverts, u_int ncells,
	         double *verts, u_int *cells, int *celladj, u_char *data);
      ~Datasetvol() {}

      					// member access methods

      float getMin(int t) const { return(vol[t]->getMin()); }
      float getMax(int t) const { return(vol[t]->getMax()); }
      float getMin()      const { return(min[0]); }
      float getMax()      const { return(max[0]); }

      // add by fan
      float getMinFun(int j)      const { return(min[j]); }
      float getMaxFun(int j)      const { return(max[j]); }

      Data	*getData(int i) { return(vol[i]); }
      Datavol	*getVol(int i) { return(vol[i]); }
};

//------------------------------------------------------------------------
//
// Datasetvol() - usual constructor, reads data from one or more files
//
//------------------------------------------------------------------------

inline Datasetvol::Datasetvol(Data::DataType t, int nd, int nt, char *fn[])
			      : Dataset(t, nd, nt, fn)
{
   int i, j;

   meshtype = 3;
   vol = (Datavol **)malloc(sizeof(Datavol *)*nt);
   for (j = 0; j < nd; j++)
       {
       min[j] = 1e10;
       max[j] = -1e10;
       }
   ncells = 0;
   for (i=0; i<nt; i++)
      {
      vol[i] = new Datavol(t, nd, fn[i]);
      for (j = 0; j < nd; j++)
	  {
	  if (vol[i]->getMin() < min[j])
	     min[j] = vol[i]->getMin();
	  if (vol[i]->getMax() > max[j])
	     max[j] = vol[i]->getMax();
	  }
      if (vol[i]->getNCells() > ncells)
         ncells = vol[i]->getNCells();
      }
   maxcellindex=ncells;
}

//------------------------------------------------------------------------
//
// Datasetvol() - called by the constructors to initialize the data
//
//------------------------------------------------------------------------

inline Datasetvol::Datasetvol(Data::DataType t, int ndata, int ntime,
			      u_int nverts, u_int ncells, double *verts,
			      u_int *cells, int *celladj, u_char *data)
			      : Dataset(t, ndata, ntime, data)
{
    int		i;			// timestep index variable
    int		j;			// a variable index
    int		size = 0;	// size of single timestep of data

    meshtype = 3;
    vol = (Datavol **)malloc(sizeof(Datavol *)*ntime);
    for (j = 0; j < ndata; j++)
	{
	min[j] = 1e10;
	max[j] = -1e10;
	}
    // ncells = 0;	this was here to allow different ncells for different
    //			times, for now we don't allow it with this constructor
    Datasetvol::ncells = ncells;

    switch (t)
	{
	case Data::UCHAR :	size = nverts * ndata * sizeof(u_char);
				break;
	case Data::USHORT :	size = nverts * ndata * sizeof(u_short);
				break;
	case Data::FLOAT :	size = nverts * ndata * sizeof(float);
				break;
	}

    for (i=0; i<ntime; i++)
	{
	vol[i] = new Datavol(t, ndata, nverts, ncells, verts, cells, celladj,
							      data + i*size);
	for (j = 0; j < ndata; j++)
	    {
	    if (vol[i]->getMin() < min[j])
		min[j] = vol[i]->getMin();
	    if (vol[i]->getMax() > max[j])
		max[j] = vol[i]->getMax();
	    }
	if (vol[i]->getNCells() > ncells)
            ncells = vol[i]->getNCells();
      }
    maxcellindex=ncells;
}

#endif
