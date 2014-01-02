/*------------------------------------------------------------------------
 *
 * contour.h - contour library include file
 *
 * Copyright (c) 1999 Emilio Camahort
 *
 *----------------------------------------------------------------------*/

/* $Id: contour.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $ */

#ifndef _CONTOUR_H_
#define _CONTOUR_H_

#ifndef WIN32
#include <unistd.h>
#endif
typedef unsigned char	u_char;
typedef unsigned short	u_short;
typedef unsigned int	u_int;

/*------------------------------------------------------------------------
 * 
 *   constant definitions
 * 
 *----------------------------------------------------------------------*/

#define CONTOUR_UCHAR	0		/* supported types of raw data */
#define CONTOUR_USHORT	1
#define CONTOUR_FLOAT	2

#define CONTOUR_2D	2		/* supported types of meshes */
#define CONTOUR_3D	3
#define CONTOUR_REG_2D	4
#define CONTOUR_REG_3D	5

#define NO_COLOR_VARIABLE	-1	/* do not color an isocontour */

/*------------------------------------------------------------------------
 * 
 *   data structures
 * 
 *----------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   data structures for storing some dataset information
 *----------------------------------------------------------------------*/

typedef struct DatasetInfo
{
    int		datatype;		/* the dataset's base type */
    int		meshtype;		/* type of mesh (see defs above) */

    int		nvars;			/* number of variables */
    int		ntime;			/* number of timesteps */

    u_int	dim[3];			/* dataset's x, y, z dimensions */
    float	orig[3];		/* coordinates of dataset's origin */
    float	span[3];		/* how far does the dataset go */

    float	minext[3];		/* min corner of dataset */
    float	maxext[3];		/* max corner of dataset */

    float	*minvar;		/* min value for each variable */
    float	*maxvar;		/* max value for each variable */
} DatasetInfo;

/*------------------------------------------------------------------------
 *   data structures for storing seed cells
 *----------------------------------------------------------------------*/

typedef struct Seed		/* store one seed cell (see seedCells.h) */
{
    float	min;			/* cell minimum value */
    float	max;			/* cell maximum value */
    u_int	cell_id;		/* cell identifier */
} Seed;

typedef struct SeedData		/* structure to store seeds */
{
    int		nseeds;			/* number of seed cells */
    Seed	*seeds;			/* array of seed cells */
} SeedData;

/*------------------------------------------------------------------------
 *   data structure for storing signature functions
 *----------------------------------------------------------------------*/

typedef struct Signature	/* a single signature function */
{
    char	*name;			/* fuction name */
    u_int	nval;			/* number of function values */
    float	*fx;			/* x-values of function */
    float	*fy;			/* y-values of function */
} Signature;

/*------------------------------------------------------------------------
 *   a slice of a 3D regular grid
 *----------------------------------------------------------------------*/

typedef struct SliceData	/* structure to store a 2d scalar data slice */
{
    int		width;			/* slice width */
    int		height;			/* slice height */
    int		datatype;		/* type of data: uchar, ushort, float */

    u_char	*ucdata;		/* slice data for uchar datatype */
    u_short	*usdata;		/* slice data for ushort datatype */
    float	*fdata;			/* slice data for float datatype */
} SliceData;

/*------------------------------------------------------------------------
 *   store a 2D isocontour
 *----------------------------------------------------------------------*/

typedef struct Contour2dData	/* structure to store a 2d isocontour */
{
    int		nvert;			/* size of vertex array */
    int		nedge;			/* size of edge array */

    float	(*vert)[2];		/* polyline vertex array */
    u_int	(*edge)[2];		/* array of edge indices */
} Contour2dData;

/*------------------------------------------------------------------------
 *   store a 3D isocontour
 *----------------------------------------------------------------------*/

typedef struct Contour3dData	/* structure to store a 3d isocontour */
{
    int		nvert;			/* number of vertices in vertex array */
    int		ntri;			/* triangles in triangle array */

    float	(*vert)[3];		/* isosurface vertex array */
    float	(*vnorm)[3];		/* array of vertex normals */
    float	(*vfun);		/* color values at vertices */

    u_int	(*tri)[3];		/* triangle mesh array */

    int		colorvar;		/* color variable (-1 if no color) */
    float	fmin, fmax;		/* min and max color values */
} Contour3dData;

/*------------------------------------------------------------------------
 *   a structure to store a dataset of scalar data
 *----------------------------------------------------------------------*/

class Dataset;
class Conplot;

typedef struct ConDataset	// structure to store and entire dataset
{
    char	**vnames;		// variable names

    int		nsfun;			// number of signature functions
    Signature	***sfun;		// signature functions

    Dataset	*data;			// input data/dataset
    Conplot	*plot;			// isocontouring algorithm and data
} ConDataset;

/*------------------------------------------------------------------------
 *
 *  library routines
 *
 *----------------------------------------------------------------------*/

		/* set the library's level of verbosity, default is 1 */

void		setVerboseLevel(int level);

		/* set the library's error handler */

void		setErrorHandler(void (*handler)(char *, int));

		/* create a new dataset for unstructured data */

ConDataset	*newDatasetUnstr(int datatype, int meshtype, int nvars, 
			         int ntime, int nverts, int ncells, 
				 double *verts, u_int *cells,
				 int *celladj, u_char *data);

		/* create a new dataset for a regular grid of data */

ConDataset	*newDatasetReg(int datatype, int meshtype, int nvars, 
			       int ntime, int *dim, u_char *data);

		/* load from disk a new dataset */

ConDataset	*loadDataset(int datatype, int meshtype, int nvars, 
			    int ntime, char **files);

		/* get the dataset's basic information */

DatasetInfo	*getDatasetInfo(ConDataset *dataset);

		/* get variable names for multi-variate datasets */

char		**getVariableNames(ConDataset *dataset);

		/* extract seed cells for variable and timestep */

SeedData	*getSeedCells(ConDataset *dataset, int variable, int timestep);

		/* return the number of signature functions */

int		getNumberOfSignatures(ConDataset *dataset);

		/* compute signature functions for variable and timestep */

Signature	*getSignatureFunctions(ConDataset *dataset, int variable, 
			    int timestep);

		/* get signature values for isovalue, variable and timestep */

float		*getSignatureValues(ConDataset *dataset, int variable, 
			    int timestep, float isovalue);

		/* extract a dataset slice */

SliceData 	*getSlice(ConDataset *dataset, int variable, int timestep, 
			    char axis, u_int index);

		/* extract a 2d isocontour */

Contour2dData	*getContour2d(ConDataset *dataset, int variable, int timestep, 
			    float isovalue);

		/* extract a 3d isocontour */

Contour3dData	*getContour3d(ConDataset *dataset, int variable, int timestep, 
			    float isovalue, int colorvar);

		/* extract a 2d isocontour and write it to a file */

void		saveContour2d(ConDataset *dataset, int variable, int timestep, 
			    float isovalue, char *filename);

		/* extract a 3d isocontour and write it to a file */

void		saveContour3d(ConDataset *dataset, int variable, int timestep, 
			    float isovalue, int colorvar, char *filename);

		/* write isocontour connected components to ipoly files */

void    	writeIsoComponents(ConDataset *dataset, int variable,
			    int timestep, float isovalue, int colorvar,
			    char *fprefix);

		/* clear (delete) dataset from memory */

void		clearDataset(ConDataset *dataset);

#endif /* of _CONTOUR_H_ */
