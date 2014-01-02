//--------------------------------------------------------------------
//
// Dataset - representation for a time-varying dataset
//
// Copyright (c) 1997 Dan Schikore - updated by Emilio Camahort, 1999
//
//--------------------------------------------------------------------

// $Id: dataset.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef DATASET_H
#define DATASET_H

#include "data.h"

//--------------------------------------------------------------------
//
// Dataset - a scalar time-varying dataset
//
//--------------------------------------------------------------------
class Dataset
{
   private:			// data members

      Data::DataType	type;		// data type: uchar, ushort, float
      int		ntime;		// number of timesteps
      int		ndata;		// add by fan
      char		**filenames;	// data filenames

   protected:

      u_int	ncells;     		// number of cells
      int	meshtype;		// 2d unstr, reg2, 3d unstr, reg3
      int	maxcellindex;		// maximum number of cells
      float	*min, *max;		// min/max values for each variable

   public:			// constructors and destructors

      Dataset(Data::DataType t, int ndata, int ntime, char *files[]);
      Dataset(Data::DataType t, int ndata, int ntime, u_char *data);
      virtual ~Dataset() {}

				// member access methods

      Data::DataType	dataType(void) const	{ return(type); }
      int		meshType(void) const	{ return(meshtype); }
      int		nTime(void) const	{ return(ntime); }
      int		nData(void) const	{ return(ndata); }
      char		**fileNames(void) const	{ return(filenames); }

      virtual float getMin(int t) const = 0;	// min, max for "0" variable
      virtual float getMax(int t) const = 0;	// at time step "t"
      virtual float getMin() const = 0;		// min, max for "0" variable
      virtual float getMax() const = 0;		// over all times

      // add by fan
      virtual float getMinFun(int f) const = 0; // min, max for "j" variable
      virtual float getMaxFun(int f) const = 0;	// over all times (reg3 only)
      // end fan

      virtual Data *getData(int i) = 0;

      u_int         getNCells(void) { return(ncells); }
      int           maxCellIndex(void) { return(maxcellindex); }
};

//------------------------------------------------------------------------
//
// Dataset() - the usual constructor, initializes some data
//
//------------------------------------------------------------------------

inline Dataset::Dataset(Data::DataType t, int nd, int nt, char *fn[])
{
    type      = t;			// base type of the data
    ntime     = nt;			// number of time step
    ndata     = nd;			// number of data: add by fan
    filenames = fn;			// filenames containing the data
}

//------------------------------------------------------------------------
//
// Dataset() - alternative constructor for the libcontour library
//
//------------------------------------------------------------------------

inline Dataset::Dataset(Data::DataType t, int nd, int nt, u_char *data)
{
    type      = t;
    ntime     = nt;
    ndata     = nd;
    filenames = NULL;
}

#endif
