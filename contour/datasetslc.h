//--------------------------------------------------------------------
//
// Datasetslc - representation for a time-varying volume
//
// Copyright (c) 1997 Dan Schikore - modified by Emilio Camahort, 1999
//
//--------------------------------------------------------------------

// $Id: datasetslc.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef DATASET_SLC_H
#define DATASET_SLC_H

#include "dataset.h"
#include "dataslc.h"

//--------------------------------------------------------------------
//
// Datasetslc - a scalar time-varying dataset
//
//--------------------------------------------------------------------

class Datasetslc : public Dataset
{
   private:				// data member

      Dataslc	**slc;

   public:				// constructors and destructors

      Datasetslc(Data::DataType t, int ndata, int ntime, char *files[]);
      Datasetslc(Data::DataType t, int ndata, int ntime, u_int nverts, u_int ncells,
	         double *verts, u_int *cells, int *celladj, u_char *data);
      virtual ~Datasetslc() {}

      					// member access methods

      float getMin(int t) const { return(slc[t]->getMin()); }
      float getMax(int t) const { return(slc[t]->getMax()); }
      float getMin()      const { return(min[0]); }
      float getMax()      const { return(max[0]); }

      // add by fan
      float getMinFun(int j)      const { return(min[j]); }
      float getMaxFun(int j)      const { return(max[j]); }


      Data	*getData(int i) { return(slc[i]); }
      Dataslc	*getSlc(int i) { return(slc[i]); }
};

//------------------------------------------------------------------------
//
// Datasetslc() - usual constructor, reads data from one or more files
//
//------------------------------------------------------------------------

inline Datasetslc::Datasetslc(Data::DataType t, int nd, int nt, char *fn[])
		       : Dataset(t, nd, nt, fn)
{
   int		i, j;
   
   meshtype = 2;
   slc = (Dataslc **)malloc(sizeof(Dataslc *)*nt);
   for (j = 0; j < nd; j++)
       {
       min[j] = 1e10;
       max[j] = -1e10;
       }
   ncells = 0;
   for (i=0; i<nt; i++)
      {
if (verbose)
printf("loading file: %s\n", fn[i]);
      slc[i] = new Dataslc(t, nd, fn[i]);
      for (j = 0; j < nd; j++)
	  {
	  if (slc[i]->getMin() < min[j])
	     min[j] = slc[i]->getMin();
	  if (slc[i]->getMax() > max[j])
	     max[j] = slc[i]->getMax();
	  }
      if (slc[i]->getNCells() > ncells)
         ncells = slc[i]->getNCells();
if (verbose) {
printf("step %d: min : %f max : %f\n", i, min[0], max[0]);
printf("step %d: tmin : %f tmax : %f\n", i, slc[i]->getMin(), slc[i]->getMax());
}
      }
   maxcellindex=ncells;
   if (verbose)
	for (i = 0; i < nd; i++)
	    printf("variable[%d]: min=%f, max=%f\n",i, min[i],max[i]);
}

//------------------------------------------------------------------------
//
// Datasetslc() - alternative constructor for the libcontour library
//
//------------------------------------------------------------------------

inline Datasetslc::Datasetslc(Data::DataType t, int ndata, int ntime, 
			      u_int nverts, u_int ncells, double *verts,
			      u_int *cells, int *celladj, u_char *data)
			      : Dataset(t, ndata, ntime, data)
{
    int	i;				// timestep index variable
    int	j;				// a variable index
    int	size = 0;		// size of single timestep of data
   
    meshtype = 2;
    slc = (Dataslc **)malloc(sizeof(Dataslc *)*ntime);
    for (j = 0; j < ndata; j++)
	{
    	min[j] = 1e10;
    	max[j] = -1e10;
	}
    //   ncells = 0;	this was here to allow different ncells for different
    //			times, for now we don't allow it with this constructor
    Datasetslc::ncells = ncells;

    switch (t)
	{
	case Data::UCHAR : 	size = nverts * ndata * sizeof(u_char);
				break;
	case Data::USHORT :	size = nverts * ndata * sizeof(u_short);
				break;
	case Data::FLOAT :	size = nverts * ndata * sizeof(float);
				break;
	}

    for (i = 0; i < ntime; i++)
	{
	slc[i] = new Dataslc(t, ndata, nverts, ncells, verts, cells, celladj,
							      data + i*size);
	for (j = 0; j < ndata; j++)
	    {
	    if (slc[i]->getMin() < min[j])
		min[j] = slc[i]->getMin();
	    if (slc[i]->getMax() > max[j])
		max[j] = slc[i]->getMax();
	    }
	if (slc[i]->getNCells() > ncells)
            ncells = slc[i]->getNCells();
if (verbose) {
printf("step %d: min : %f max : %f\n", i, min[0], max[0]);
printf("step %d: tmin : %f tmax : %f\n", i, slc[i]->getMin(), slc[i]->getMax());
}
        }

   maxcellindex=ncells;
   if (verbose)
	for (i = 0; i < ndata; i++)
	    printf("variable[%d]: min=%f, max=%f\n",i, min[i],max[i]);
}

#endif
