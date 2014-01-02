//------------------------------------------------------------------------
//
// cellQueue.C - queue of cell identifiers.  The circular queue dyanmically
//               resizes itself when full.  Elements in the queue are of
//               indicies of i,j,k (specialized for 3d structured grids)
//
// Copyright (c) 1997 Dan Schikore
//------------------------------------------------------------------------

// $Id: cellqueue.cpp,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#ifndef MACOS_X
#include <malloc.h>
#endif
#include "cellqueue.h"
