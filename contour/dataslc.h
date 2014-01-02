//--------------------------------------------------------------------
//
// Dataslc - class for a triangular slice of scalar data
//
// Copyright (c) 1997 Dan Schikore
//--------------------------------------------------------------------

// $Id: dataslc.h,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#ifndef SLICE_H
#define SLICE_H

#include "data.h"

//--------------------------------------------------------------------
//
// Dataslc - a volume of scalar data.
//
//--------------------------------------------------------------------
class Dataslc : public Data
{
   private:				// data members

      double (*verts)[2];		// array of mesh vertices
      float (*vgrad)[3];		// function gradients at vertices
      u_int (*cells)[3];		// cells (triangles) in the mesh
      int   (*celladj)[3];		// indices to adjacent triangles

   public:				// constructors and destructors

      Dataslc(Data::DataType t, int ndata, char *rawfile=NULL);
      Dataslc(Data::DataType t, int ndata, u_int nverts, u_int ncells, 
	      double *verts, u_int *cells, int *celladj, u_char *data);
      inline ~Dataslc();

   					// signature function methods

      int getNFunctions(void)          { return(4); }
      float *compFunction(int, u_int &, float **);
      float *compFunction(int, u_int &, float ***,
                          float ***, float ***){return(NULL);} // add by fan
 
      char *fName(int);

   protected :				// signature functions

      float *compLength(u_int &, float **);
      float *compArea(u_int &, float **);
      float *compMaxArea(u_int &, float **);
      float *compGradient(u_int &, float **);

   public :

      // get data or gradient approximations (by differencing)
      u_int *getCell(int c) const
           { return(cells[c]); }
      int *getCellAdjs(int c) const
           { return(celladj[c]); }
      double *getVert(int v) const
           { return(verts[v]); }
      void getCellValues(int c, float *val)
           { val[0] = getValue(cells[c][0]);
             val[1] = getValue(cells[c][1]);
             val[2] = getValue(cells[c][2]);
           }
      u_int *getCellVerts(int c) const
           { return(cells[c]); }
      u_int   getCellVert(int c, int v) {return(cells[c][v]); }

      u_int getNCellVerts(void) { return(3); }
      u_int getNCellFaces(void) { return(3); }
      int getCellAdj(int c, int f) { return(celladj[c][f]); }

      void getCellGrad(int c, float grad[2]) {
         double u[3], v[3];
         u[0] = verts[cells[c][1]][0] - verts[cells[c][0]][0];
         u[1] = verts[cells[c][1]][1] - verts[cells[c][0]][1];
         u[2] = getValue(cells[c][1]) - getValue(cells[c][0]);
         v[0] = verts[cells[c][2]][0] - verts[cells[c][0]][0];
         v[1] = verts[cells[c][2]][1] - verts[cells[c][0]][1];
         v[2] = getValue(cells[c][2]) - getValue(cells[c][0]);
         grad[0] = (float)(u[1]*v[2] - u[2]*v[1]);
         grad[1] = (float)(u[2]*v[0] - u[0]*v[2]);
      }

      void getCellGrad3(int c, float grad[3]) {
         double u[3], v[3];
         u[0] = verts[cells[c][1]][0] - verts[cells[c][0]][0];
         u[1] = verts[cells[c][1]][1] - verts[cells[c][0]][1];
         u[2] = getValue(cells[c][1]) - getValue(cells[c][0]);
         v[0] = verts[cells[c][2]][0] - verts[cells[c][0]][0];
         v[1] = verts[cells[c][2]][1] - verts[cells[c][0]][1];
         v[2] = getValue(cells[c][2]) - getValue(cells[c][0]);
         grad[0] = (float)(u[1]*v[2] - u[2]*v[1]);
         grad[1] = (float)(u[2]*v[0] - u[0]*v[2]);
         grad[2] = (float)(u[0]*v[1] - u[1]*v[0]);
      }

      void normalToFace(int c, int f, float norm[2]) {
         double u[2];
         u[0] = verts[cells[c][f]][0] - verts[cells[c][f==2?0:f+1]][0];
         u[1] = verts[cells[c][f]][1] - verts[cells[c][f==2?0:f+1]][1];
         norm[0] = (float)-u[1];
         norm[1] =  (float)u[0];
      }

      void getCellRange(int c, float &min, float &max)
           {
              float t;
              max = min = getValue(cells[c][0]);
              if ((t=getValue(cells[c][1])) < min)
                 min = t;
              if (t > max)
                 max = t;
              if ((t=getValue(cells[c][2])) < min)
                 min = t;
              if (t > max)
                 max = t;
           }
      void getFaceRange(u_int c, u_int f, float &min, float &max)
           {
              float t;
              min = max = getValue(cells[c][f]);
              if ((t=getValue(cells[c][f==2?0:f+1])) < min)
                 min = t;
              if (t > max)
                 max = t;
           }

};

//------------------------------------------------------------------------
//
// ~Dataslc() - destroy a volume
//
//------------------------------------------------------------------------

inline Dataslc::~Dataslc()
{
    if (filename)
	{
	free(verts);
	free(cells);
	free(celladj);
	}
}

#endif
