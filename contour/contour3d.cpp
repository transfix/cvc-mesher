//  ______________________________________________________________________
//
//    NAME
//      Contour3d - Class for a contour surface
//
//      Copyright (c) 1998 Emilio Camahort, Dan Schikore
//
//    SYNOPSIS
//      #include <contour3d.h>
//  ______________________________________________________________________

// $Id: contour3d.cpp,v 1.1.1.1 2006/10/11 21:25:51 transfix Exp $

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef MACOS_X
#include <malloc.h>
#endif

#include "contour3d.h"

//#define WRITE_TMESH
//#define WRITE_IPOLY
#define WRITE_RAW

/******** Added by Joe R - 12/21/2005, used for fast searching of duplicate vertices  ********/
#ifdef USEDICT
struct vertex_dict_entry
{
  float x, y, z; // vertex coords
  int i;         // vertex index
};

static int VertexCompare(const void *a, const void *b)
{
  const vertex_dict_entry *v1 = (vertex_dict_entry*)a;
  const vertex_dict_entry *v2 = (vertex_dict_entry*)b;
  
  if(fabs(v1->x-v2->x) >= 0.00001) return v1->x < v2->x ? -1 : 1;
  if(fabs(v1->y-v2->y) >= 0.00001) return v1->y < v2->y ? -1 : 1;
  if(fabs(v1->z-v2->z) >= 0.00001) return v1->z < v2->z ? -1 : 1;

  return 0;
}

static dnode_t *VertexNodeAlloc(void *c)
{
  return dnode_create(NULL);
}

static void VertexNodeDelete(dnode_t *n, void *c)
{
  free((vertex_dict_entry*)dnode_get(n));
  free(n);
}
#endif
/********************************************************************************************/

//------------------------------------------------------------------------
//
// Contour3d() - basic constructor
//
//------------------------------------------------------------------------
Contour3d::Contour3d(int _vf)
{
   done = 0;

   nvert = 0;
   ntri = 0;
   vf = _vf;

   vsize =  500;
   tsize = 1000;

   vert  = (float (*)[3])malloc(sizeof(float[3]) * vsize);
   vnorm = (float (*)[3])malloc(sizeof(float[3]) * vsize);
   tri   = (u_int (*)[3])malloc(sizeof(u_int[3]) * tsize);

   vfun = (float *)malloc(sizeof(float)*vsize);

   color = (vf>1);

#ifdef USEDICT    /* added by Joe R. 12/21/2005 */
   dict_init(&vertex_dict,DICTCOUNT_T_MAX,VertexCompare);
   dict_set_allocator(&vertex_dict,VertexNodeAlloc,VertexNodeDelete,NULL);
#endif
}

//------------------------------------------------------------------------
//
// ~Contour3d() - free allocated memory
//
//------------------------------------------------------------------------
Contour3d::~Contour3d()
{
   free(vert);
   free(vnorm);
   free(tri);
#ifdef USEDICT
   dict_free(&vertex_dict); /* added by Joe R. 12/21/2005 - Deletes all vertex dictionary entries */
#endif
}

//------------------------------------------------------------------------
//
// AddVert() - add a vertex with the given (unit) normal
//
//------------------------------------------------------------------------
int
Contour3d::AddVert(float x, float y, float z, float nx, float ny, float nz, float f)
{
   int n = nvert++;

   if (nvert > vsize) {
      vsize<<=1;
      vert  = (float (*)[3])realloc(vert, sizeof(float[3]) * vsize);
      vnorm = (float (*)[3])realloc(vnorm, sizeof(float[3]) * vsize);
      vfun = (float *)realloc(vfun, sizeof(float)*vsize);
   }

   vert[n][0] = x;
   vert[n][1] = y;
   vert[n][2] = z;

   vnorm[n][0] = nx;
   vnorm[n][1] = ny;
   vnorm[n][2] = nz;

   vfun[n] = f;
// printf("f = %f\n", f);

#ifdef USEDICT   /* added by Joe R. 12/21/2005 */
   vertex_dict_entry *vde = (vertex_dict_entry*)malloc(sizeof(vertex_dict_entry));
   vde->x = x; vde->y = y; vde->z = z; vde->i = n;
   dict_alloc_insert(&vertex_dict, vde, vde);
#endif

   return(n);
}

//------------------------------------------------------------------------
//
// AddVertUnique() - add a vertex with the given (unit) normal
//
//------------------------------------------------------------------------
int
Contour3d::AddVertUnique(float x, float y, float z, float nx, float ny, float nz, float f)
{
#ifdef USEDICT  /* added by Joe R. 12/21/2005 */
   vertex_dict_entry vde;
   vde.x = x; vde.y = y; vde.z = z;
   dnode_t *n = dict_lookup(&vertex_dict,&vde);
   if(n != NULL)
     return ((vertex_dict_entry*)dnode_get(n))->i;
   else
     return(AddVert(x,y,z,nx,ny,nz,f));
#else
   int i;

   for (i=nvert-1; i>=0; i--) {
      if ((fabs(vert[i][0]-x) < 0.00001) &&
          (fabs(vert[i][1]-y) < 0.00001) &&
          (fabs(vert[i][2]-z) < 0.00001))
         return(i);
   }
#endif
   return(AddVert(x,y,z,nx,ny,nz,f));
}

//------------------------------------------------------------------------
//
// AddTri() - add a triangle indexed by it's 3 vertices
//
//------------------------------------------------------------------------
int
Contour3d::AddTri(u_int v1, u_int v2, u_int v3)
{
   int n = ntri++;

   if (ntri > tsize) {
      tsize<<=1;
      tri = (u_int (*)[3])realloc(tri, sizeof(u_int[3]) * tsize);
   }

   tri[n][0] = v1;
   tri[n][1] = v2;
   tri[n][2] = v3;

   return(n);
}

//------------------------------------------------------------------------
//
// Reset() - clear vertex and surface info
//
//------------------------------------------------------------------------
void
Contour3d::Reset(void)
{
   nvert = 0;
   ntri  = 0;
   done = 0;
#ifdef USEDICT
   dict_free(&vertex_dict); /* added by Joe R. 12/21/2005 - Deletes all vertex dictionary entries */
#endif
}


void
Contour3d::Done(void)
{
   done = 1;
}

//------------------------------------------------------------------------
//
// write() - write vertex and triangles to a file
//
//------------------------------------------------------------------------

int	Contour3d::write(char *filename)
{
#ifdef WRITE_IPOLY
   FILE *fp;
   int v, t;
                                                                                                                                            
   fp = fopen(filename, "w");
                                                                                                                                            
   // silent failure: changed by Emilio --> return 1 = ERROR
   if (fp == NULL)
      return 1;
                                                                                                                                            
   fprintf(fp, "%d 0 %d 0 0 0 0\n0 0 0\n", nvert, ntri);
                                                                                                                                            
   for (v=0; v<nvert; v++)
      fprintf(fp, "%g %g %g\n", vert[v][0], vert[v][1], vert[v][2]);
                                                                                                                                            
   fprintf(fp, "0 0\n");
                                                                                                                                            
   for (t=0; t<ntri; t++)
      fprintf(fp, "3\n%d %d %d\n", tri[t][0], tri[t][1], tri[t][2]);
                                                                                                                                            
   fclose(fp);
#elif defined WRITE_TMESH

   FILE *fp;
   int t;
   int nvf, nmat;
   u_int _tri[3];
   float vn[3];

   fp = fopen(filename, "w");

   // silent failure: changed by Emilio --> return 1 = ERROR
   if (fp == NULL)
      return 1;

   nvf = 1;
   nmat = 0;

   fwrite(&nvert, sizeof(int), 1, fp);
   fwrite(&ntri,  sizeof(int), 1, fp);
   fwrite(&nvf,   sizeof(int), 1, fp);
   fwrite(&nmat,  sizeof(int), 1, fp);

   fwrite(vert, sizeof(float[3]), nvert, fp);
   for (t=0; t<nvert; t++) {
      vn[0] = -vnorm[t][0];
      vn[1] = -vnorm[t][1];
      vn[2] = -vnorm[t][2];
     fwrite(vn, sizeof(float[3]), 1, fp);
   }

   // reverse the triangles
   for (t=0; t<ntri; t++) {
      _tri[0] = tri[t][1];
      _tri[1] = tri[t][0];
      _tri[2] = tri[t][2];
      fwrite(_tri, sizeof(int[3]), 1, fp);
   }

   // write the vertex function
   fwrite(vfun, sizeof(float), nvert, fp);

   fclose(fp);
#elif defined WRITE_RAW
   FILE *fp;
   int v, t;

   fp = fopen(filename, "w");

   if (fp == NULL)
      return 1;

   fprintf(fp, "%d %d\n", nvert, ntri);

   for (v=0; v<nvert; v++)
      fprintf(fp, "%7.3f %7.3f %7.3f\n", vert[v][0], vert[v][1], vert[v][2]);

   for (t=0; t<ntri; t++)
      fprintf(fp, "%d %d %d\n", tri[t][0], tri[t][1], tri[t][2]);

   fclose(fp);
#else /* write a simple poly file */
   FILE *fp;
   int t;

   fp = fopen(filename,"w");

   if(fp == NULL) return 1;

   fprintf(fp, "%d\n", ntri);
   for (t=0; t<ntri; t++) {
      fprintf(fp, "3\n%g %g %g\n",
              vert[tri[t][0]][0], vert[tri[t][0]][1], vert[tri[t][0]][2]);
      fprintf(fp, "%g %g %g\n",
              vert[tri[t][2]][0], vert[tri[t][2]][1], vert[tri[t][2]][2]);
      fprintf(fp, "%g %g %g\n",
              vert[tri[t][1]][0], vert[tri[t][1]][1], vert[tri[t][1]][2]);
   }

   fclose(fp);
#endif

   return 0;
}
