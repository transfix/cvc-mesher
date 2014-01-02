/******************************************************************************
				Copyright   

This code is developed within the Computational Visualization Center at The 
University of Texas at Austin.

This code has been made available to you under the auspices of a Lesser General 
Public License (LGPL) (http://www.ices.utexas.edu/cvc/software/license.html) 
and terms that you have agreed to.

Upon accepting the LGPL, we request you agree to acknowledge the use of use of 
the code that results in any published work, including scientific papers, 
films, and videotapes by citing the following references:

Y. Zhang, C. Bajaj
Adaptive and Quality Quadrilateral/Hexahedral Meshing from Volumetric Data
Computer Methods in Applied Mechanics and Engineering (CMAME), 
195(9-12):942-960, 2006.

Y. Zhang, C. Bajaj, B-S. Sohn
3D Finite Element Meshing from Imaging Data
The special issue of Computer Methods in Applied Mechanics and Engineering 
(CMAME) on Unstructured Mesh Generation, 194(48-49):5083-5106, 2005. 

Y. Zhang, C. Bajaj, G. Xu
Surface Smoothing and Quality Improvement of Quadrilateral/Hexahedral Meshes 
with Geometric Flow
Proceedings of 14th International Meshing Roundtable, pp. 449-468. San Diego, 
CA. September 11-14, 2005.  

If you desire to use this code for a profit venture, or if you do not wish to 
accept LGPL, but desire usage of this code, please contact Chandrajit Bajaj 
(bajaj@ices.utexas.edu) at the Computational Visualization Center at The 
University of Texas at Austin for a different license.
******************************************************************************/

#include "LBIE_geoframe.h"
#include <cstdio>
#include <cmath>
#include <cstring>

			 //#include <algorithm>
			 //#include <boost/multi_array.hpp>

#if 0
#include <qfiledialog.h>
#include <qfile.h>
#include <qfileinfo.h>
#include <qdir.h>
#include <qtextstream.h>
#endif

namespace LBIE
{

geoframe::geoframe()
  : avg_aspect(0.0), max_aspect(0.0), min_aspect(0.0),
    numverts(0), numtris(0), num_tris(0), numquads(0),
    numhexas(0), biggestDim(0.0), centerx(0.0), centery(0.0),
    centerz(0.0), max_x(0.0), min_x(0.0), max_y(0.0), min_y(0.0),
    max_z(0.0), min_z(0.0), mesh_type(SINGLE)
{
  for(int i=0; i<3; i++) span[i] = 1.0;
}

geoframe::geoframe(const geoframe& geofrm)
  : avg_aspect(geofrm.avg_aspect), max_aspect(geofrm.max_aspect),
    min_aspect(geofrm.min_aspect), numverts(geofrm.numverts),
    numtris(geofrm.numtris), num_tris(geofrm.num_tris), numquads(geofrm.numquads),
    numhexas(geofrm.numhexas), verts(geofrm.verts), normals(geofrm.normals),
    color(geofrm.color), curvatures(geofrm.curvatures), funcs(geofrm.funcs),
    triangles(geofrm.triangles), quads(geofrm.quads), bound_sign(geofrm.bound_sign),
    bound_tri(geofrm.bound_tri), vtx_idx_arr_extend(geofrm.vtx_idx_arr_extend),
    vtxnew_sign(geofrm.vtxnew_sign), bound_edge(geofrm.bound_edge),
    refine_edge(geofrm.refine_edge), refine_edgevtx(geofrm.refine_edgevtx),
    biggestDim(geofrm.biggestDim), centerx(geofrm.centerx), centery(geofrm.centery),
    centerz(geofrm.centerz), max_x(geofrm.max_x), min_x(geofrm.min_x),
    max_y(geofrm.max_y), min_y(geofrm.min_y), max_z(geofrm.max_z), min_z(geofrm.min_z),
    mesh_type(SINGLE)
{
  for(int i=0; i<3; i++) span[i] = geofrm.span[i];
}

geoframe::~geoframe()
{
}

geoframe& geoframe::operator=(const geoframe& geofrm)
{
  avg_aspect = geofrm.avg_aspect;
  max_aspect = geofrm.max_aspect;
  min_aspect = geofrm.min_aspect;
  numverts = geofrm.numverts;
  numtris = geofrm.numtris;
  num_tris = geofrm.num_tris;
  numquads = geofrm.numquads;
  numhexas = geofrm.numhexas;
  verts = geofrm.verts;
  normals = geofrm.normals;
  color = geofrm.color;
  curvatures = geofrm.curvatures;
  funcs = geofrm.funcs;
  triangles = geofrm.triangles;
  quads = geofrm.quads;
  bound_sign = geofrm.bound_sign;
  bound_tri = geofrm.bound_tri;
  vtx_idx_arr_extend = geofrm.vtx_idx_arr_extend;
  vtxnew_sign = geofrm.vtxnew_sign;
  bound_edge = geofrm.bound_edge;
  refine_edge = geofrm.refine_edge;
  refine_edgevtx = geofrm.refine_edgevtx;
  biggestDim = geofrm.biggestDim;
  centerx = geofrm.centerx;
  centery = geofrm.centery;
  centerz = geofrm.centerz;
  max_x = geofrm.max_x;
  min_x = geofrm.min_x;
  max_y = geofrm.max_y;
  min_y = geofrm.min_y;
  max_z = geofrm.max_z;
  min_z = geofrm.min_z;
  for(int i=0; i<3; i++) span[i] = geofrm.span[i];
  mesh_type = geofrm.mesh_type;
  return *this;
}

/*
void geoframe::loadgeoframe(const char * name, int num)
{
  char buffer[256];
  
  sprintf(buffer, "%s%05d.raw", name, num );
  //sprintf(buffer, "%s", name);
  QFile file(buffer);

  file.open(IO_ReadOnly);
  QTextStream stream(&file);

  stream >> numverts >> numtris;
  verts = new double[numverts*3];
  triangles = new unsigned int[numtris*3];
  int c;
  //int min = 1<<30;
  for (c=0; c<numverts; c++) {
    stream >> verts[c*3+0];
    stream >> verts[c*3+1];
    stream >> verts[c*3+2];
  }
  for (c=0; c<numtris; c++) {
    stream >> triangles[c*3+0];
    //min = obj.faceArray[c*3+0]<min?obj.faceArray[c*3+0]:min;
    stream >> triangles[c*3+1];
    //min = obj.faceArray[c*3+0]<min?obj.faceArray[c*3+1]:min;
    stream >> triangles[c*3+2];
    //min = obj.faceArray[c*3+0]<min?obj.faceArray[c*3+2]:min;
  }
  calculatenormals();
  calculateExtents();
}
*/

void cross(float* dest, const float* v1, const float* v2)
{
	dest[0] = v1[1]*v2[2] - v1[2]*v2[1];
	dest[1] = v1[2]*v2[0] - v1[0]*v2[2];
	dest[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

void geoframe::calculateAspectRatio()
{
  avg_aspect = max_aspect = 0.0;
  //min_aspect = 

  printf("calculate aspect ratio\n");
  int vert1,vert2,vert3;
  float tmp_aspect = 0.0;
  for (int c=0; c<numtris; c++) {
    vert1 = triangles[c][0];
    vert2 = triangles[c][1];
    vert3 = triangles[c][2];
    tmp_aspect = get_aspect_ratio(vert1,vert2,vert3);
    if(tmp_aspect > max_aspect)
      max_aspect = tmp_aspect;
    avg_aspect += tmp_aspect;
  }
  avg_aspect /= (float)numtris;
  printf("max aspect ratio:%f avg aspect ratio:%f\n",max_aspect,avg_aspect);
}


void geoframe::calculateTriangleNormal(float* norm, unsigned int c)
{
	
  float v1[3], v2[3];
  int vert;
  vert = triangles[c][0];
  v1[0] = v2[0] = -verts[vert][0];
  v1[1] = v2[1] = -verts[vert][1];
  v1[2] = v2[2] = -verts[vert][2];

  vert = triangles[c][1];
  v1[0] += verts[vert][0];
  v1[1] += verts[vert][1];
  v1[2] += verts[vert][2];

  vert = triangles[c][2];
  v2[0] += verts[vert][0];
  v2[1] += verts[vert][1];
  v2[2] += verts[vert][2];

  cross(norm, v1, v2);
  
}

void geoframe::calculatenormals()
{
	
	
	int c, vert;
	float normal[3];
	float len;
	
	// for each triangle
	for (c=0; c<numtris; c++) {
		calculateTriangleNormal(normal, c);
		normals[c][0] = -normal[0];
		normals[c][1] = -normal[1];
		normals[c][2] = -normal[2];
	}
	
	// normalize the vectors
	for (vert=0; vert<numtris; vert++) {
		len = (float) sqrt(
			normals[vert][0] * normals[vert][0] +
			normals[vert][1] * normals[vert][1] +
			normals[vert][2] * normals[vert][2]);
		if(len == 0.0)
		{
			normals[vert][0] = 1.0;
			normals[vert][1] = 0.0;
			normals[vert][2] = 0.0;
			printf("error\n"); continue;
		}
		normals[vert][0]/=len;
		normals[vert][1]/=len;
		normals[vert][2]/=len;
	}
	
}

void geoframe::calculateExtents()
{
  int c;
  max_x=0.f, min_x=0.f;
  max_y=0.f, min_y=0.f;
  max_z=0.f, min_z=0.f;
  float value;

  for (c=0; c<numverts; c++) {
    if (c==0) {
      max_x = min_x = verts[c][0];
      max_y = min_y = verts[c][1];
      max_z = min_z = verts[c][2];
    }
    else {
      value = verts[c][0];
      max_x = (value>max_x?value:max_x);
      min_x = (value<min_x?value:min_x);

      value = verts[c][1];
      max_y = (value>max_y?value:max_y);
      min_y = (value<min_y?value:min_y);

      value = verts[c][2];
      max_z = (value>max_z?value:max_z);
      min_z = (value<min_z?value:min_z);
    }
  }

    biggestDim = (max_y-min_y>max_x-min_x?max_y-min_y:max_x-min_x);
    biggestDim = (max_z-min_z>biggestDim?max_z-min_z:biggestDim);
    centerx = (max_x+min_x)/2.0;
    centery = (max_y+min_y)/2.0;
    centerz = (max_z+min_z)/2.0;
}

/*void geoframe::display()
{
  int vert;
  glBegin(GL_TRIANGLES);
  for (int c=0; c<numtris; c++) {
    vert = triangles[c*3+0];
    glNormal3d( 
      normals[vert*3+0],
      normals[vert*3+1],
      normals[vert*3+2]);
    glVertex3d(
      verts[vert*3+0],
      verts[vert*3+1],
      verts[vert*3+2]);
    vert = triangles[c*3+1];
    glNormal3d( 
      normals[vert*3+0],
      normals[vert*3+1],
      normals[vert*3+2]);
    glVertex3d(
      verts[vert*3+0],
      verts[vert*3+1],
      verts[vert*3+2]);
    vert = triangles[c*3+2];
    glNormal3d( 
      normals[vert*3+0],
      normals[vert*3+1],
      normals[vert*3+2]);
    glVertex3d(
      verts[vert*3+0],
      verts[vert*3+1],
      verts[vert*3+2]);
  }
  glEnd();
}

void geoframe::read_raw(const char * rawiv_fname) {
	
	float temp0, temp1, temp2;
	int nv, ntri, t0, t1, t2;
	int fscanf_return = 0;

	FILE* vol_fp = fopen(rawiv_fname,"r");
	
	if (vol_fp==NULL) {
		printf("wrong name : %s\n",rawiv_fname);
		return -1;
	}

	fscanf_return = fscanf(vol_fp,"%d %d\n", &nv, &ntri);
	numverts = nv;
	numtris = ntri;

	verts  = (float (*)[3])malloc(sizeof(float[3]) * nv);
	triangles   = (int (*)[3])malloc(sizeof(int[3]) * ntri);
	int i;
	for (i=0;i<nv;i++) {
		fscanf_return = fscanf(vol_fp,"%f %f %f\n", &temp0, &temp1, &temp2);
		verts[i][0] = temp0;
		verts[i][1] = temp1;
		verts[i][2] = temp2;
	}
	for (i=0;i<ntri;i++) {
		fscanf_return = fscanf(vol_fp,"%d %d %d\n", &t0, &t1, &t2);
		triangles[i][0] = t0;
		triangles[i][1] = t1;
		triangles[i][2] = t2;
	}
	fclose(vol_fp);

}*/

int geoframe::read_raw(const char * rawiv_fname) {
	
  	float temp0, temp1, temp2, r, g, b;
	int nv, t0, t1, t2, t3, flag_mesh;
	int i, j, k, sign, vv0, vv1, vv2, bool_0;
	//int **neighbor;
	typedef std::vector<std::vector<unsigned int> > neighbor_array_type;
	float nx, ny, nz;

	reset();

	FILE* vol_fp = fopen(rawiv_fname, "r");
	if (vol_fp==NULL) {
		printf("wrong name : %s\n",rawiv_fname);
		return -1;
	}

	flag_mesh = -1;
	if(strstr(rawiv_fname, "_tri") != NULL) flag_mesh = 0;
	if(strstr(rawiv_fname, "_tet") != NULL) flag_mesh = 1;
	if(strstr(rawiv_fname, "_quad") != NULL) flag_mesh = 2;
	if(strstr(rawiv_fname, "_hex") != NULL) flag_mesh = 3;
	if(strstr(rawiv_fname, "_nurbs") != NULL) flag_mesh = 4;

	nx = 0.0f;	ny = 0.0f;	nz = 0.0f;

	switch (flag_mesh) {

	// read triangular mesh
	case 0:
	  {
		int ntri;
		fscanf(vol_fp,"%d %d\n", &nv, &ntri);
		numverts = nv;
		numtris = ntri;

		//verts  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		verts.resize(nv);
		vtx_idx_arr_extend.resize(nv,-1);
		//normals  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		normals.resize(nv);
		//color  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		color.resize(nv);
		//triangles   = (unsigned int (*)[3])malloc(sizeof(unsigned int[3]) * ntri);
		triangles.resize(ntri);
		//bound_sign = (unsigned int (*))realloc(bound_sign,sizeof(unsigned int)*vsize);
		bound_sign.resize(nv);
		//bound_tri = (unsigned int (*))realloc(bound_tri,sizeof(unsigned int) * ntri);
		bound_tri.resize(ntri);

		normals[0][0] = 9999.0f;
		for (i=0;i<nv;i++) {
			if(strstr(rawiv_fname, ".rawnc") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f %f %f %f\n", &temp0, &temp1, &temp2, &nx, &ny, &nz, &r, &g, &b);
			else if(strstr(rawiv_fname, ".rawn") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f\n", &temp0, &temp1, &temp2, &nx, &ny, &nz);
			else if(strstr(rawiv_fname, ".rawc") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f\n", &temp0, &temp1, &temp2, &r, &g, &b);
			else if(strstr(rawiv_fname, ".raw") != NULL)
				fscanf(vol_fp,"%f %f %f\n", &temp0, &temp1, &temp2);
			else {
				printf("wrong name : %s\n",rawiv_fname);
				return flag_mesh;
			}

			if(strstr(rawiv_fname, ".rawc") == NULL && strstr(rawiv_fname, ".rawnc") == NULL) {
				//r = 1.0f;	g = 0.3f; b = 0.3f;
				//r = 0.3f;	g = 1.0f; b = 0.3f;
				//r = 1.0f;	g = 0.0f; b = 1.0f;
				//r = 0.6f;	g = 0.6f; b = 1.0f;
				//if(sign == 0) {
					r = 1.0f;	g = 0.7f; b = 0.0f;
					//r = 1.0f;	g = 0.4f; b = 0.4f;
					//r = 0.0f;	g = 0.4f; b = 0.0f;
				//}
			}

			verts[i][0] = temp0;	verts[i][1] = temp1;	verts[i][2] = temp2;
			if(strstr(rawiv_fname, ".rawn") != NULL || strstr(rawiv_fname, ".rawnc") != NULL) {
				normals[i][0] = -nx;	normals[i][1] = -ny;	normals[i][2] = -nz;
			}
			color[i][0] = r;	color[i][1] = g;	color[i][2] = b;
			bound_sign[i] = 1;
		}
		for (i=0;i<ntri;i++) {
			fscanf(vol_fp,"%d %d %d\n", &t0, &t1, &t2);
			triangles[i][0] = t0;
			triangles[i][1] = t1;
			triangles[i][2] = t2;
			bound_tri[i] = 1;
		}
		break;
	  }

	// read tetra mesh
	case 1:	
	  {
	        unsigned int ntet, itet;
		fscanf(vol_fp,"%d %d\n", &nv, &ntet);
		numverts = nv;
		numtris = ntet*4;

		//verts  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		verts.resize(nv);
		vtx_idx_arr_extend.resize(nv,-1);
		//normals  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		normals.resize(nv);
		//color  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		color.resize(nv);
		//triangles   = (unsigned int (*)[3])malloc(sizeof(unsigned int[3]) * ntet*4);
		triangles.resize(numtris);
		//bound_sign = (unsigned int (*))realloc(bound_sign,sizeof(unsigned int)*vsize);
		bound_sign.resize(nv);
		//bound_tri = (unsigned int (*))realloc(bound_tri,sizeof(unsigned int) * ntet*4);
		bound_tri.resize(numtris);

#if 0
		neighbor = 0;
		neighbor = new int*[nv];
		for(i = 0; i < nv; i++)	neighbor[i] = 0;
		for(i = 0; i < nv; i++)	neighbor[i] = new int[50];
		for (i = 0; i < nv; i++) {
			for (j = 0; j < 50; j++) neighbor[i][j] = -1;
		}
#endif
		neighbor_array_type neighbor(nv);

		normals[0][0] = 9999.0f;
		for (i = 0; i < nv; i++) {
			if(strstr(rawiv_fname, ".rawnc") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f %f %f %f %d\n", &temp0, &temp1, &temp2, &nx, &ny, &nz, &r, &g, &b, &sign);
			else if(strstr(rawiv_fname, ".rawn") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f %d\n", &temp0, &temp1, &temp2, &nx, &ny, &nz, &sign);
			else if(strstr(rawiv_fname, ".rawc") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f %d\n", &temp0, &temp1, &temp2, &r, &g, &b, &sign);
			else if(strstr(rawiv_fname, ".raw") != NULL)
				fscanf(vol_fp,"%f %f %f %d\n", &temp0, &temp1, &temp2, &sign);
			else {
				printf("wrong name : %s\n",rawiv_fname);
				return flag_mesh;
			}

			if(strstr(rawiv_fname, ".rawc") == NULL && strstr(rawiv_fname, ".rawnc") == NULL) {
				//r = 1.0f;	g = 0.3f; b = 0.3f;
				//r = 0.3f;	g = 1.0f; b = 0.3f;
				//r = 1.0f;	g = 0.0f; b = 1.0f;
				//r = 0.6f;	g = 0.6f; b = 1.0f;
				//if(sign == 0) {
					r = 1.0f;	g = 0.7f; b = 0.0f;
					//r = 0.0f;	g = 0.4f; b = 0.0f;
				//}
			}

			verts[i][0] = temp0+32.0f;
			verts[i][1] = temp1+32.0f;
			verts[i][2] = temp2+32.0f;
			if(strstr(rawiv_fname, ".rawn") != NULL || strstr(rawiv_fname, ".rawnc") != NULL) {
				normals[i][0] = nx;	normals[i][1] = ny;	normals[i][2] = nz;
			}
			color[i][0] = r;	color[i][1] = g;	color[i][2] = b;
			bound_sign[i] = sign;
		}
		for (i=0;i<ntet;i++) {
			fscanf(vol_fp,"%d %d %d %d\n", &t0, &t1, &t2, &t3);
			triangles[4*i][0] = t0;		triangles[4*i][1] = t2;		triangles[4*i][2] = t1;
			triangles[4*i+1][0] = t1;	triangles[4*i+1][1] = t2;	triangles[4*i+1][2] = t3;
			triangles[4*i+2][0] = t0;	triangles[4*i+2][1] = t3;	triangles[4*i+2][2] = t2;
			triangles[4*i+3][0] = t0;	triangles[4*i+3][1] = t1;	triangles[4*i+3][2] = t3;

			// find which face is on the boundary
#if 0
			for(j = 0; j < 50; j++) {
				if(neighbor[t0][j] == -1) {neighbor[t0][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t1][j] == -1) {neighbor[t1][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t2][j] == -1) {neighbor[t2][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t3][j] == -1) {neighbor[t3][j] = i;	break;}
			}
#endif
			neighbor[t0].push_back(i);
			neighbor[t1].push_back(i);
			neighbor[t2].push_back(i);
			neighbor[t3].push_back(i);
		}

		for (i = 0; i < ntet; i++) {
			for(j = 0; j < 4; j++) {
				vv0 = triangles[4*i+j][0];	vv1 = triangles[4*i+j][1];	vv2 = triangles[4*i+j][2];	
				
				bound_tri[4*i+j] = 0;

				//				for(k = 0; k < 50; k++) {
				//					itet = neighbor[vv0][k];
				//					if(itet == -1) break;
				for(neighbor_array_type::value_type::iterator itet = neighbor[vv0].begin();
				    itet != neighbor[vv0].end();
				    itet++)
				  {
				    t0 = triangles[4*(*itet)][0];	t1 = triangles[4*(*itet)][1];
				    t2 = triangles[4*(*itet)][2];	t3 = triangles[4*(*itet)+1][2];

					bool_0 = 0;
					if(t0 == vv0 || t1 == vv0 || t2 == vv0 || t3 == vv0) bool_0++;
					if(t0 == vv1 || t1 == vv1 || t2 == vv1 || t3 == vv1) bool_0++;
					if(t0 == vv2 || t1 == vv2 || t2 == vv2 || t3 == vv2) bool_0++;

					if(bool_0 == 3) bound_tri[4*i+j]++; 
				}
			}
		}

#if 0
		for(i = 0; i < nv; i++)	delete [] neighbor[i];
		delete [] neighbor;
#endif

		break;
	  }

    // read quad mesh
	case 2:
	  {
		int nquad;
		fscanf(vol_fp,"%d %d\n", &nv, &nquad);
		numverts = nv;
		numquads = nquad;

		//verts  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		verts.resize(nv);
		vtx_idx_arr_extend.resize(nv,-1);
		//normals  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		normals.resize(nv);
		//color  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		color.resize(nv);
		//quads   = (unsigned int (*)[4])malloc(sizeof(unsigned int[4]) * nquad);
		quads.resize(nquad);
		//bound_sign = (unsigned int (*))realloc(bound_sign,sizeof(unsigned int)*nv);
		bound_sign.resize(nv);

		normals[0][0] = 9999.0f;
		for (i = 0; i < nv; i++) {
			if(strstr(rawiv_fname, ".rawnc") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f %f %f %f\n", &temp0, &temp1, &temp2, &nx, &ny, &nz, &r, &g, &b);
			else if(strstr(rawiv_fname, ".rawn") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f\n", &temp0, &temp1, &temp2, &nx, &ny, &nz);
			else if(strstr(rawiv_fname, ".rawc") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f\n", &temp0, &temp1, &temp2, &r, &g, &b);
			else if(strstr(rawiv_fname, ".raw") != NULL)
				fscanf(vol_fp,"%f %f %f\n", &temp0, &temp1, &temp2);
			else {
				printf("wrong name : %s\n",rawiv_fname);
				return flag_mesh;
			}

			if(strstr(rawiv_fname, ".rawc") == NULL && strstr(rawiv_fname, ".rawnc") == NULL) {
				//r = 1.0f;	g = 0.59f;	b = 0.59f;
				r = 1.0f;	g = 0.7f;	b = 0.0f;
				//r = 1.0f;	g = 0.3f;	b = 0.3f;
			}

			verts[i][0] = temp0;	verts[i][1] = temp1;	verts[i][2] = temp2;
			if(strstr(rawiv_fname, ".rawn") != NULL || strstr(rawiv_fname, ".rawnc") != NULL) {
				normals[i][0] = nx;	normals[i][1] = ny;	normals[i][2] = nz;
			}
			color[i][0] = r;		color[i][1] = g;		color[i][2] = b;
			bound_sign[i] = 1;
		}
		for (i = 0;i < nquad;i++) {
			fscanf(vol_fp,"%d %d %d %d\n", &t0, &t1, &t2, &t3);
			quads[i][0] = t0;		quads[i][1] = t3;
			quads[i][2] = t2;		quads[i][3] = t1;
			//bound_quad[i] = 1;
		}
		break;
	  }

	// read hex mesh
	case 3:	
	  {
		int nhex, t4, t5, t6, t7;
		int sign0, sign1, sign2, sign3, sign4, sign5;
		fscanf(vol_fp,"%d %d\n", &nv, &nhex);
		numverts = nv;
		numquads = nhex*6;
		numhexas = nhex;

		//verts  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		verts.resize(nv);
		vtx_idx_arr_extend.resize(nv,-1);
		//normals  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		normals.resize(nv);
		//color  = (float (*)[3])malloc(sizeof(float[3]) * nv);
		color.resize(nv);
		//quads   = (unsigned int (*)[4])malloc(sizeof(unsigned int[4]) * nhex*6);
		quads.resize(numquads);
		//bound_sign = (unsigned int (*))realloc(bound_sign,sizeof(unsigned int) * nv);
		bound_sign.resize(nv);
		//bound_quad = (unsigned int (*))realloc(bound_quad,sizeof(unsigned int) * nhex*6);

#if 0
		neighbor = 0;
		neighbor = new int*[nv];
		for(i = 0; i < nv; i++)	neighbor[i] = 0;
		for(i = 0; i < nv; i++)	neighbor[i] = new int[50];
		for (i = 0; i < nv; i++) {
			for (j = 0; j < 50; j++) neighbor[i][j] = -1;
		}
#endif
		neighbor_array_type neighbor(nv);
		
		//float minz = 10.0f, nn;
		normals[0][0] = 9999.0f;
		for (i = 0; i < nv; i++) {
			if(strstr(rawiv_fname, ".rawnc") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f %f %f %f %d\n", &temp0, &temp1, &temp2, &nx, &ny, &nz, &r, &g, &b, &sign);
			else if(strstr(rawiv_fname, ".rawn") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f %d\n", &temp0, &temp1, &temp2, &nx, &ny, &nz, &sign);
			else if(strstr(rawiv_fname, ".rawc") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f %d\n", &temp0, &temp1, &temp2, &r, &g, &b, &sign);
			else if(strstr(rawiv_fname, ".raw") != NULL)
				fscanf(vol_fp,"%f %f %f %d\n", &temp0, &temp1, &temp2, &sign);
				//fscanf(vol_fp,"%f %f %f\n", &temp0, &temp1, &temp2);
			else {
				printf("wrong name : %s\n",rawiv_fname);
				return flag_mesh;
			}

			//sign = 1;

			if(strstr(rawiv_fname, ".rawc") == NULL && strstr(rawiv_fname, ".rawnc") == NULL) {
				//r = 1.0f;	g = 0.59f;	b = 0.59f;
				//r = 1.0f;	g = 0.7f;	b = 0.0f;
				r = 0.6f;	g = 0.6f;	b = 1.0f;
				//if(i < 7522) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				//else {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			}

			//nn = (float)nz/sqrt(nx*nx + ny*ny + nz*nz);
			//if(temp2 < -9.0f && fabs(nn) > 0.9f) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			//else if(temp1 < 0.0f && temp2 < 5.6f && fabs(nn) > 0.9f) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			//else if(temp2 > 12.0f && fabs(nn) > 0.9f) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			//else {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			//if(temp2 < minz) minz = temp2;

			//nn = (float)nz/sqrt(nx*nx + ny*ny + nz*nz);
			//if(temp2 < -150.0f && fabs(nn) > 0.9f) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			//else if(temp2 > 122.0f && fabs(nn) > 0.95f) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			//else {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			//if(temp2 < minz) minz = temp2;

			if(sign > 1) {r = 1.0f;	g = 0.59f;	b = 0.59f; sign = 1;}
			else {
				
				r = 1.0f;	g = 0.7f;	b = 0.0f;
				//r = 0.6f;	g = 0.6f;	b = 1.0f;
				/*
				// template no < 7
				//r = 1.0f;	g = 0.3f;	b = 0.3f;
				if(i < 27) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
				else if(i < 66) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
				//else {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				else if(i < 96) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				else {r = 0.3f;	g = 0.7f;	b = 1.0f;}
				
				// template no == 7
				if(i < 27) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
				else if(i < 39) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
				else if(i < 78) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
				else if(i < 123) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				else {r = 0.3f;	g = 0.7f;	b = 1.0f;}
				// template no == 7c
				if(i < 27) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
				else if(i < 39) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
				else if(i < 78) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
				else if(i < 114) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				else {r = 0.3f;	g = 0.7f;	b = 1.0f;}
				// template no == 8, 8c
				if(i < 3*17) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
				else if(i < 7*17+6) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
				else if(i < 7*17+4*9+6) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				else {r = 0.3f;	g = 0.7f;	b = 1.0f;}
				
				// template no == 9, 9c
				if(i < 3*17) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
				else if(i < 7*17+12) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
				else if(i < 7*17+4*9+12) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				else if(i < 7*17+4*9*2+12) {r = 1.0f;	g = 0.6f;	b = 0.6f;}
				else if(i < 7*17+4*9*3+12) {r = 0.3f;	g = 0.7f;	b = 1.0f;}
				else {r = 1.0f;	g = 0.0f;	b = 0.7f;}
				
				// template no == 10, 10c
				if(i < 27) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
				else if(i < 66) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
				else if(i < 96) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				else if(i < 123) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
				else {r = 0.3f;	g = 0.7f;	b = 1.0f;}

				// template no == 11
				if(i < 3*9) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
				else if(i < 4*9+3) {r = 0.6f;	g = 0.6f;	b = 1.0f;}
				else if(i < 5*9+6) {r = 1.0f;	g = 0.6f;	b = 0.6f;}
				else if(i < 9*9+9) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
				else if(i < 9*9+5*9+9) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				else if(i < 9*9+5*9*2+9) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
				else {r = 0.3f;	g = 0.7f;	b = 1.0f;}

				// template no == 11c
				if(i < 3*9) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
				else if(i < 4*9+3) {r = 0.6f;	g = 0.6f;	b = 1.0f;}
				else if(i < 5*9+6) {r = 1.0f;	g = 0.6f;	b = 0.6f;}
				else if(i < 9*9+9) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
				else if(i < 9*9+4*9+9) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				else if(i < 9*9+4*9*2+9) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
				else {r = 0.3f;	g = 0.7f;	b = 1.0f;}
				
				// template no == 12, 12c
				if(i < 3*17) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
				else if(i < 7*17+12) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
				else if(i < 7*17+4*9+12) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				else if(i < 7*17+4*9*2+15) {r = 1.0f;	g = 0.6f;	b = 0.6f;}
				else if(i < 7*17+4*9*3+15) {r = 0.3f;	g = 0.7f;	b = 1.0f;}
				else if(i < 7*17+4*9*4+15) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
				else {r = 0.6f;	g = 0.6f;	b = 1.0f;}
				*/
			}

			verts[i][0] = temp0;	verts[i][1] = temp1;	verts[i][2] = temp2;
			if(strstr(rawiv_fname, ".rawn") != NULL || strstr(rawiv_fname, ".rawnc") != NULL) {
				normals[i][0] = nx;	normals[i][1] = ny;	normals[i][2] = nz;
			}
			color[i][0] = r;		color[i][1] = g;		color[i][2] = b;
			bound_sign[i] = sign;
		}
		for (i = 0; i < nhex; i++) {
			//fscanf(vol_fp,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7,
			//	&sign0, &sign1, &sign2, &sign3, &sign4, &sign5);
			//fscanf(vol_fp,"%d %d %d %d %d %d %d %d %d\n", &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &sign);
			//fscanf(vol_fp,"%d %d %d %d %d %d %d %d\n", &t0, &t1, &t3, &t2, &t4, &t5, &t7, &t6);
			fscanf(vol_fp,"%d %d %d %d %d %d %d %d\n", &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
			//if(i > 6146) {
			//	t0 += 7522;	t1 += 7522;	t2 += 7522;	t3 += 7522;
			//	t4 += 7522;	t5 += 7522;	t6 += 7522;	t7 += 7522;
			//}

			quads[6*i][0] = t0;		quads[6*i][1] = t3; 	quads[6*i][2] = t2;		quads[6*i][3] = t1;
			quads[6*i+1][0] = t4;	quads[6*i+1][1] = t5; 	quads[6*i+1][2] = t6;	quads[6*i+1][3] = t7;
			quads[6*i+2][0] = t0;	quads[6*i+2][1] = t4; 	quads[6*i+2][2] = t7;	quads[6*i+2][3] = t3;
			quads[6*i+3][0] = t1;	quads[6*i+3][1] = t2; 	quads[6*i+3][2] = t6;	quads[6*i+3][3] = t5;
			quads[6*i+4][0] = t0;	quads[6*i+4][1] = t1; 	quads[6*i+4][2] = t5;	quads[6*i+4][3] = t4;
			quads[6*i+5][0] = t2;	quads[6*i+5][1] = t3; 	quads[6*i+5][2] = t7;	quads[6*i+5][3] = t6;

			//sign0 = 1;	sign1 = 1;	sign2 = 1;	sign3 = 1;	sign4 = 1;	sign5 = 1;
			//bound_quad[6*i] = sign0;	bound_quad[6*i+1] = sign1;	bound_quad[6*i+2] = sign2;
			//bound_quad[6*i+3] = sign3;	bound_quad[6*i+4] = sign4;	bound_quad[6*i+5] = sign5;
			/*
			for(j = 0; j < 6; j++) {
				if(bound_sign[quads[6*i+j][0]]+bound_sign[quads[6*i+j][1]]+bound_sign[quads[6*i+j][2]]+bound_sign[quads[6*i+j][3]]==4)
					bound_quad[6*i+j] = 1;
				else
					bound_quad[6*i+j] = 0;
			}
			*/
#if 0
			for(j = 0; j < 50; j++) {
				if(neighbor[t0][j] == -1) {neighbor[t0][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t1][j] == -1) {neighbor[t1][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t2][j] == -1) {neighbor[t2][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t3][j] == -1) {neighbor[t3][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t4][j] == -1) {neighbor[t4][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t5][j] == -1) {neighbor[t5][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t6][j] == -1) {neighbor[t6][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t7][j] == -1) {neighbor[t7][j] = i;	break;}
			}
#endif
			neighbor[t0].push_back(i);
			neighbor[t1].push_back(i);
			neighbor[t2].push_back(i);
			neighbor[t3].push_back(i);
			neighbor[t4].push_back(i);
			neighbor[t5].push_back(i);
			neighbor[t6].push_back(i);
			neighbor[t7].push_back(i);
		}

		// find which face is on the boundary
#if 0
		int vv3, ihex;
		for (i = 0; i < nhex; i++) {
			for(j = 0; j < 6; j++) {
				vv0 = quads[6*i+j][0];	vv1 = quads[6*i+j][1];
				vv2 = quads[6*i+j][2];	vv3 = quads[6*i+j][3];
				//bound_quad[6*i+j] = 0;
				//				for(k = 0; k < 50; k++) {
				//					ihex = neighbor[vv0][k];
				//					if(ihex == -1) break;
				for(neighbor_array_type::value_type::iterator ihex = neighbor[vv0].begin();
				    ihex != neighbor[vv0].end();
				    ihex++)
				  {
					t0 = quads[6*(*ihex)][0];		t1 = quads[6*(*ihex)][1];
					t2 = quads[6*(*ihex)][2];		t3 = quads[6*(*ihex)][3];
					t4 = quads[6*(*ihex)+1][0];	t5 = quads[6*(*ihex)+1][1];
					t6 = quads[6*(*ihex)+1][2];	t7 = quads[6*(*ihex)+1][3];

					int bool_0 = 0;
					if( t0 == vv0 || t1 == vv0 || t2 == vv0 || t3 == vv0 ||
						t4 == vv0 || t5 == vv0 || t6 == vv0 || t7 == vv0 ) bool_0++;
					if( t0 == vv1 || t1 == vv1 || t2 == vv1 || t3 == vv1 ||
						t4 == vv1 || t5 == vv1 || t6 == vv1 || t7 == vv1 ) bool_0++;
					if( t0 == vv2 || t1 == vv2 || t2 == vv2 || t3 == vv2 ||
						t4 == vv2 || t5 == vv2 || t6 == vv2 || t7 == vv2 ) bool_0++;
					if( t0 == vv3 || t1 == vv3 || t2 == vv3 || t3 == vv3 ||
						t4 == vv3 || t5 == vv3 || t6 == vv3 || t7 == vv3 ) bool_0++;

					//if(bool_0 == 4) bound_quad[6*i+j]++;
				}
			}
		}
#endif

#if 0
		for(i = 0; i < nv; i++)	delete [] neighbor[i];
		delete [] neighbor;
#endif

		break;
	  }
	// read nurbs mesh
	case 4:	
	  {
	    int nhex, ntri, t4, t5, t6, t7;
		//int sign0, sign1, sign2, sign3, sign4, sign5;
	  int nv1, npatch/*, *patch*/;
	  npatch = 0;
		fscanf(vol_fp,"%d %d %d %d %d\n", &nv, &nhex, &nv1, &ntri, &npatch);
		numverts = nv+nv1;
		numquads = nhex*6;
		numhexas = nhex;
		numtris = ntri;

		//verts  = (float (*)[3])malloc(sizeof(float[3]) * numverts);
		verts.resize(numverts);
		vtx_idx_arr_extend.resize(nv,-1);
		//normals  = (float (*)[3])malloc(sizeof(float[3]) * numverts);
		normals.resize(numverts);
		//color  = (float (*)[3])malloc(sizeof(float[3]) * numverts);
		color.resize(numverts);
		//quads   = (unsigned int (*)[4])malloc(sizeof(unsigned int[4]) * nhex*6);
		quads.resize(numquads);
		//triangles   = (unsigned int (*)[3])malloc(sizeof(unsigned int[3]) * ntri);
		triangles.resize(numtris);
		//bound_sign = (unsigned int (*))realloc(bound_sign,sizeof(unsigned int) * numverts);
		bound_sign.resize(numverts);
		//bound_quad = (unsigned int (*))realloc(bound_quad,sizeof(unsigned int) * nhex*6);
		//bound_tri = (unsigned int (*))realloc(bound_tri,sizeof(unsigned int) * ntri);
		bound_tri.resize(ntri);
		//patch = (int(*))malloc(sizeof(int)*npatch);
		std::vector<int> patch(npatch);

		for(i = 0; i < npatch; i++) {
			fscanf(vol_fp,"%d\n", &t0);
			patch[i] = t0;
		}

#if 0
		neighbor = 0;
		neighbor = new int*[nv];
		for(i = 0; i < nv; i++)	neighbor[i] = 0;
		for(i = 0; i < nv; i++)	neighbor[i] = new int[50];
		for (i = 0; i < nv; i++) {
			for (j = 0; j < 50; j++) neighbor[i][j] = -1;
		}
#endif
		neighbor_array_type neighbor(nv);

		normals[0][0] = 9999.0f;
		for (i = 0; i < nv; i++) {
			if(strstr(rawiv_fname, ".rawnc") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f %f %f %f %d\n", &temp0, &temp1, &temp2, &nx, &ny, &nz, &r, &g, &b, &sign);
			else if(strstr(rawiv_fname, ".rawn") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f %d\n", &temp0, &temp1, &temp2, &nx, &ny, &nz, &sign);
			else if(strstr(rawiv_fname, ".rawc") != NULL)
				fscanf(vol_fp,"%f %f %f %f %f %f %d\n", &temp0, &temp1, &temp2, &r, &g, &b, &sign);
			else if(strstr(rawiv_fname, ".raw") != NULL)
				fscanf(vol_fp,"%f %f %f %d\n", &temp0, &temp1, &temp2, &sign);
			else {
				printf("wrong name : %s\n",rawiv_fname);
				return flag_mesh;
			}

			//if(temp2 == 0.0f) nx = 1.0f;

			if(strstr(rawiv_fname, ".rawc") == NULL && strstr(rawiv_fname, ".rawnc") == NULL) {
				//if(i >= patch[0] && i < patch[5]) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
				//else if(i < patch[6]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				//else if(i < patch[7]) {r = 0.6f;	g = 0.6f;	b = 1.0f;}
				//else if(i < patch[8]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
				//else {r = 0.6f;	g = 0.6f;	b = 1.0f;}
				r = 0.6f;	g = 0.6f;	b = 1.0f;
			}

			//nn = (float)nz/sqrt(nx*nx + ny*ny + nz*nz);
			//if(temp2 < -9.0f && fabs(nn) > 0.9f) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			//else if(temp1 < 0.0f && temp2 < 5.6f && fabs(nn) > 0.9f) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			//else if(temp2 > 12.0f && fabs(nn) > 0.9f) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			//else {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			//if(temp2 < minz) minz = temp2;

			//nn = (float)nz/sqrt(nx*nx + ny*ny + nz*nz);
			//if(temp2 < -150.0f && fabs(nn) > 0.9f) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			//else if(temp2 > 122.0f && fabs(nn) > 0.95f) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			//else {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			//if(temp2 < minz) minz = temp2;

			//if(sign > 1) {r = 1.0f;	g = 0.59f;	b = 0.59f; sign = 1;}
			//else {r = 0.6f;	g = 0.6f;	b = 1.0f;}

			verts[i][0] = temp0;	verts[i][1] = temp1;	verts[i][2] = temp2;
			if(strstr(rawiv_fname, ".rawn") != NULL || strstr(rawiv_fname, ".rawnc") != NULL) {
				normals[i][0] = nx;	normals[i][1] = ny;	normals[i][2] = nz;
				//normals[i][0] = nx;	normals[i][1] = ny;	normals[i][2] = nz;
			}
			color[i][0] = r;		color[i][1] = g;		color[i][2] = b;
			bound_sign[i] = sign;
		}
		for (i = 0; i < nhex; i++) {
			//fscanf(vol_fp,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7,
			//	&sign0, &sign1, &sign2, &sign3, &sign4, &sign5);
			//fscanf(vol_fp,"%d %d %d %d %d %d %d %d %d\n", &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &sign);
			//fscanf(vol_fp,"%d %d %d %d %d %d %d %d\n", &t0, &t1, &t3, &t2, &t4, &t5, &t7, &t6);
			fscanf(vol_fp,"%d %d %d %d %d %d %d %d\n", &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
			//if(i > 6146) {
			//	t0 += 7522;	t1 += 7522;	t2 += 7522;	t3 += 7522;
			//	t4 += 7522;	t5 += 7522;	t6 += 7522;	t7 += 7522;
			//}

			quads[6*i][0] = t0;		quads[6*i][1] = t3; 	quads[6*i][2] = t2;		quads[6*i][3] = t1;
			quads[6*i+1][0] = t4;	quads[6*i+1][1] = t5; 	quads[6*i+1][2] = t6;	quads[6*i+1][3] = t7;
			quads[6*i+2][0] = t0;	quads[6*i+2][1] = t4; 	quads[6*i+2][2] = t7;	quads[6*i+2][3] = t3;
			quads[6*i+3][0] = t1;	quads[6*i+3][1] = t2; 	quads[6*i+3][2] = t6;	quads[6*i+3][3] = t5;
			quads[6*i+4][0] = t0;	quads[6*i+4][1] = t1; 	quads[6*i+4][2] = t5;	quads[6*i+4][3] = t4;
			quads[6*i+5][0] = t2;	quads[6*i+5][1] = t3; 	quads[6*i+5][2] = t7;	quads[6*i+5][3] = t6;

			r = 0.7f;	g = 0.7f;	b = 1.0f;

			if(npatch == 0) {r = 1.0f;	g = 0.3f;	b = 0.3f;}
			else {
			// abdominal_Yuri
			if(i >= patch[0] && i < patch[4]) {r = 1.0f;	g = 0.3f;	b = 0.3f;}
			else if(i < patch[5]) {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			else if(i < patch[6]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[7]) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			else if(i < patch[8]) {r = 0.3f;	g = 0.7f;	b = 0.3f;}
			else if(i < patch[9]) {r = 0.7f;	g = 1.0f;	b = 0.3f;}
			else if(i < patch[10]) {r = 0.0f;	g = 0.7f;	b = 1.0f;}
			else if(i < patch[11]) {r = 1.0f;	g = 1.0f;	b = 0.4f;}
			else {r = 1.0f;	g = 0.4f;	b = 1.0f;}
		}

			r = 200.0f/255.0f;	g = 0.0f;	b = 0.0f;

			// propeller
			if(i >= patch[0] && i < patch[1]) {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			else if(i < patch[2]) {r = 1.0f;	g = 0.3f;	b = 0.3f;}
			else if(i < patch[3]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[4]) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			else if(i < patch[5]) {r = 0.0f;	g = 0.7f;	b = 1.0f;}
			else if(i < patch[6]) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
			else if(i < patch[7]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else if(i < patch[8]) {r = 0.3f;	g = 0.4f;	b = 1.0f;}
			else if(i < patch[9]) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
/*			else if(i < patch[10]) {r = 0.3f;	g = 0.7f;	b = 0.3f;}
			else if(i < patch[11]) {r = 1.0f;	g = 1.0f;	b = 0.4f;}
			else if(i < patch[12]) {r = 203.0f/255.0f;	g = 81.0f/255.0f;	b = 94.0f/255.0f;}
			else if(i < patch[13]) {r = 203.0f/255.0f;	g = 149.0f/255.0f;	b = 149.0f/255.0f;}
			else if(i < patch[14]) {r = 61.0f/255.0f;	g = 61.0f/255.0f;	b = 203.0f/255.0f;}
			else if(i < patch[15]) {r = 0.0f/255.0f;	g = 203.0f/255.0f;	b = 0.0f/255.0f;}
			else if(i < patch[16]) {r = 162.0f/255.0f;	g = 0.0f/255.0f;	b = 54.0f/255.0f;}
			else if(i < patch[17]) {r = 162.0f/255.0f;	g = 54.0f/255.0f;	b = 13.0f/255.0f;}
/*			else {r = 0.7f;	g = 1.0f;	b = 0.3f;}

			/*
			// abdominal_full
			if(i >= patch[0] && i < patch[1]) {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			else if(i < patch[2]) {r = 1.0f;	g = 0.3f;	b = 0.3f;}
			else if(i < patch[3]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[5]) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			else if(i < patch[6]) {r = 0.0f;	g = 0.7f;	b = 1.0f;}
			else if(i < patch[8]) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
			else if(i < patch[10]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else if(i < patch[12]) {r = 0.3f;	g = 0.4f;	b = 1.0f;}
			else if(i < patch[13]) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			else if(i < patch[14]) {r = 0.3f;	g = 0.7f;	b = 0.3f;}
			else if(i < patch[15]) {r = 1.0f;	g = 1.0f;	b = 0.4f;}
			else if(i < patch[16]) {r = 203.0f/255.0f;	g = 81.0f/255.0f;	b = 94.0f/255.0f;}
			else if(i < patch[17]) {r = 203.0f/255.0f;	g = 149.0f/255.0f;	b = 149.0f/255.0f;}
			else if(i < patch[18]) {r = 61.0f/255.0f;	g = 61.0f/255.0f;	b = 203.0f/255.0f;}
			else if(i < patch[19]) {r = 0.0f/255.0f;	g = 203.0f/255.0f;	b = 0.0f/255.0f;}
			else if(i < patch[20]) {r = 162.0f/255.0f;	g = 0.0f/255.0f;	b = 54.0f/255.0f;}
			else if(i < patch[21]) {r = 162.0f/255.0f;	g = 54.0f/255.0f;	b = 13.0f/255.0f;}
			else if(i < patch[22]) {r = 9.0f/255.0f;	g = 94.0f/255.0f;	b = 92.0f/255.0f;}
			else if(i < patch[24]) {r = 203.0f/255.0f;	g = 203.0f/255.0f;	b = 3.0f/255.0f;}
			else if(i < patch[25]) {r = 162.0f/255.0f;	g = 0.0f/255.0f;	b = 56.0f/255.0f;}
			else if(i < patch[27]) {r = 162.0f/255.0f;	g = 81.0f/255.0f;	b = 203.0f/255.0f;}
			else if(i < patch[28]) {r = 203.0f/255.0f;	g = 203.0f/255.0f;	b = 108.0f/255.0f;}
			else if(i < patch[29]) {r = 108.0f/255.0f;	g = 0.0f/255.0f;	b = 102.0f/255.0f;}
			else if(i < patch[32]) {r = 0.0f/255.0f;	g = 200.0f/255.0f;	b = 203.0f/255.0f;}
			else if(i < patch[33]) {r = 203.0f/255.0f;	g = 203.0f/255.0f;	b = 0.0f;}
			else {r = 0.7f;	g = 1.0f;	b = 0.3f;}
			
			//r = 1.0f;	g = 0.3f; b = 0.3f;
			//r = 0.3f;	g = 1.0f; b = 0.3f;
			// abdominal_3
			if(i >= patch[0] && i < patch[1]) {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			else if(i < patch[2]) {r = 1.0f;	g = 0.3f;	b = 0.3f;}
			else if(i < patch[3]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[4]) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			else if(i < patch[5]) {r = 0.0f;	g = 0.7f;	b = 1.0f;}
			else if(i < patch[6]) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
			else if(i < patch[7]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else if(i < patch[8]) {r = 0.3f;	g = 0.4f;	b = 1.0f;}
			else if(i < patch[9]) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			else if(i < patch[10]) {r = 0.3f;	g = 0.7f;	b = 0.3f;}
			else if(i < patch[11]) {r = 1.0f;	g = 1.0f;	b = 0.4f;}
			else {r = 0.7f;	g = 1.0f;	b = 0.3f;}

			// abdominal
			if(i >= patch[0] && i < patch[1]) {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			else if(i < patch[2]) {r = 1.0f;	g = 0.3f;	b = 0.3f;}
			else if(i < patch[3]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[4]) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			else if(i < patch[5]) {r = 0.0f;	g = 0.7f;	b = 1.0f;}
			else if(i < patch[6]) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
			else if(i < patch[7]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else if(i < patch[8]) {r = 0.3f;	g = 0.4f;	b = 1.0f;}
			else if(i < patch[9]) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			else {r = 0.3f;	g = 0.7f;	b = 0.3f;}

			// thoracic aorta
			if(i >= patch[0] && i < patch[2]) {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			else if(i < patch[4]) {r = 1.0f;	g = 0.3f;	b = 0.3f;}
			else if(i < patch[5]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[6]) {r = 1.0f;	g = 0.59f;	b = 0.59f;}
			else if(i < patch[7]) {r = 0.0f;	g = 0.7f;	b = 1.0f;}
			else if(i < patch[8]) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
			else if(i < patch[9]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else if(i < patch[10]) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			else {r = 0.3f;	g = 0.7f;	b = 0.3f;}
			
			// coronary
			if(i >= patch[0] && i < patch[1]) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			else if(i < patch[2]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[3]) {r = 1.0f;	g = 0.3f;	b = 0.3f;}
			else {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			
			// template 1-4
			if(i >= patch[0] && i < patch[1]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[2]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			
			// template 5, 6
			if(i >= patch[0] && i < patch[1]) {r = 0.3f;	g = 0.7f;	b = 1.0f;}
			else if(i < patch[2]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[3]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			
			// template 7
			if(i >= patch[0] && i < patch[1]) {r = 0.3f;	g = 0.7f;	b = 1.0f;}
			else if(i < patch[2]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[3]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else if(i < patch[4]) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
			else {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			
			// template 8
			if(i >= patch[0] && i < patch[1]) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			else if(i < patch[2]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else if(i < patch[3]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else {r = 0.3f;	g = 0.7f;	b = 1.0f;}
			
			// template 9
			if(i >= patch[0] && i < patch[1]) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			else if(i < patch[2]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else if(i < patch[3]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[4]) {r = 1.0f;	g = 0.6f;	b = 0.6f;}
			else if(i < patch[5]) {r = 0.3f;	g = 0.7f;	b = 1.0f;}
			else {r = 1.0f;	g = 0.0f;	b = 0.7f;}
			
			// template 10
			if(i >= patch[0] && i < patch[1]) {r = 0.3f;	g = 0.7f;	b = 1.0f;}
			else if(i < patch[2]) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
			else if(i < patch[3]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[4]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			// template 11
			if(i >= patch[0] && i < patch[1]) {r = 0.3f;	g = 0.7f;	b = 1.0f;}
			else if(i < patch[2]) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
			else if(i < patch[3]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[4]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else if(i < patch[5]) {r = 1.0f;	g = 0.6f;	b = 0.6f;}
			else if(i < patch[6]) {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			else {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			
			// template 12
			if(i >= patch[0] && i < patch[1]) {r = 1.0f;	g = 0.4f;	b = 1.0f;}
			else if(i < patch[2]) {r = 0.3f;	g = 1.0f;	b = 0.3f;}
			else if(i < patch[3]) {r = 1.0f;	g = 0.7f;	b = 0.0f;}
			else if(i < patch[4]) {r = 1.0f;	g = 0.6f;	b = 0.6f;}
			else if(i < patch[5]) {r = 0.3f;	g = 0.7f;	b = 1.0f;}
			else if(i < patch[6]) {r = 1.0f;	g = 0.0f;	b = 0.7f;}
			else {r = 0.6f;	g = 0.6f;	b = 1.0f;}
			*/

			/*for(j = 0; j < 4; j++) {
				color[quads[6*i][j]][0] = r;
				color[quads[6*i][j]][1] = g;
				color[quads[6*i][j]][2] = b;
				color[quads[6*i+1][j]][0] = r;
				color[quads[6*i+1][j]][1] = g;
				color[quads[6*i+1][j]][2] = b;
			}*/

			//sign0 = 1;	sign1 = 1;	sign2 = 1;	sign3 = 1;	sign4 = 1;	sign5 = 1;
			//bound_quad[6*i] = sign0;	bound_quad[6*i+1] = sign1;	bound_quad[6*i+2] = sign2;
			//bound_quad[6*i+3] = sign3;	bound_quad[6*i+4] = sign4;	bound_quad[6*i+5] = sign5;
			/*
			for(j = 0; j < 6; j++) {
				if(bound_sign[quads[6*i+j][0]]+bound_sign[quads[6*i+j][1]]+bound_sign[quads[6*i+j][2]]+bound_sign[quads[6*i+j][3]]==4)
					bound_quad[6*i+j] = 1;
				else
					bound_quad[6*i+j] = 0;
			}
			*/
#if 0
			for(j = 0; j < 50; j++) {
				if(neighbor[t0][j] == -1) {neighbor[t0][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t1][j] == -1) {neighbor[t1][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t2][j] == -1) {neighbor[t2][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t3][j] == -1) {neighbor[t3][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t4][j] == -1) {neighbor[t4][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t5][j] == -1) {neighbor[t5][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t6][j] == -1) {neighbor[t6][j] = i;	break;}
			}
			for(j = 0; j < 50; j++) {
				if(neighbor[t7][j] == -1) {neighbor[t7][j] = i;	break;}
			}
#endif
			neighbor[t0].push_back(i);
			neighbor[t1].push_back(i);
			neighbor[t2].push_back(i);
			neighbor[t3].push_back(i);
			neighbor[t4].push_back(i);
			neighbor[t5].push_back(i);
			neighbor[t6].push_back(i);
			neighbor[t7].push_back(i);
		}

		// find which face is on the boundary
#if 0
		//int vv3, ihex;
		for (i = 0; i < nhex; i++) {
			for(j = 0; j < 6; j++) {
				vv0 = quads[6*i+j][0];	vv1 = quads[6*i+j][1];
				vv2 = quads[6*i+j][2];	vv3 = quads[6*i+j][3];
				//bound_quad[6*i+j] = 0;
				//				for(k = 0; k < 50; k++) {
				//					ihex = neighbor[vv0][k];
				//					if(ihex == -1) break;
				for(neighbor_array_type::value_type::iterator ihex = neighbor[vv0].begin();
				    ihex != neighbor[vv0].end();
				    ihex++)
				  {
					t0 = quads[6*ihex][0];		t1 = quads[6*ihex][1];
					t2 = quads[6*ihex][2];		t3 = quads[6*ihex][3];
					t4 = quads[6*ihex+1][0];	t5 = quads[6*ihex+1][1];
					t6 = quads[6*ihex+1][2];	t7 = quads[6*ihex+1][3];

					int bool_0 = 0;
					if( t0 == vv0 || t1 == vv0 || t2 == vv0 || t3 == vv0 ||
						t4 == vv0 || t5 == vv0 || t6 == vv0 || t7 == vv0 ) bool_0++;
					if( t0 == vv1 || t1 == vv1 || t2 == vv1 || t3 == vv1 ||
						t4 == vv1 || t5 == vv1 || t6 == vv1 || t7 == vv1 ) bool_0++;
					if( t0 == vv2 || t1 == vv2 || t2 == vv2 || t3 == vv2 ||
						t4 == vv2 || t5 == vv2 || t6 == vv2 || t7 == vv2 ) bool_0++;
					if( t0 == vv3 || t1 == vv3 || t2 == vv3 || t3 == vv3 ||
						t4 == vv3 || t5 == vv3 || t6 == vv3 || t7 == vv3 ) bool_0++;

					//if(bool_0 == 4) bound_quad[6*i+j]++;
				}
			}
		}
#endif

#if 0
		for(i = 0; i < nv; i++)	delete [] neighbor[i];
		delete [] neighbor;
#endif

		for (i=0;i<nv1;i++) {
			fscanf(vol_fp,"%f %f %f\n", &temp0, &temp1, &temp2);

			if(strstr(rawiv_fname, ".rawc") == NULL && strstr(rawiv_fname, ".rawnc") == NULL) {
				//r = 1.0f;	g = 0.3f; b = 0.3f;
				//r = 0.3f;	g = 1.0f; b = 0.3f;
				//r = 1.0f;	g = 0.0f; b = 1.0f;
				//r = 0.6f;	g = 0.6f; b = 1.0f;
				//if(sign == 0) {
					r = 1.0f;	g = 0.7f; b = 0.0f;
					//r = 1.0f;	g = 0.4f; b = 0.4f;
					//r = 0.0f;	g = 0.4f; b = 0.0f;
				//}
			}

			verts[i+nv][0] = temp0;	verts[i+nv][1] = temp1;	verts[i+nv][2] = temp2;
			//color[i+nv][0] = r;	color[i+nv][1] = g;	color[i+nv][2] = b;
			bound_sign[i+nv] = 1;
		}
		for (i=0;i<ntri;i++) {
			fscanf(vol_fp,"%d %d %d\n", &t0, &t1, &t2);
			triangles[i][0] = t0+nv;
			triangles[i][1] = t1+nv;
			triangles[i][2] = t2+nv;
			bound_tri[i] = 1;
		}

		break;
	  }
	}

	fclose(vol_fp);
	mesh_type = GEOTYPE(flag_mesh);
	return flag_mesh;
}


/*void geoframe::duplicate()
{
	int nv = numverts;
	int nquad = numquads;

	for (i = 0; i < nv; i++) {
                
    	    verts[i+nv][0] = verts[i][0];    verts[i+nv][1] = verts[i][1] + 10.0;    verts[i+nv][2] = verts[i][2];
            normals[i+nv][0] = normals[i][0];    normals[i+nv][1] = normals[i][1];    normals[i+nv][2] = normals[i][2];
            bound_sign[i+nv] = 1;
        }
        for (i = 0;i < nquad;i++) {
                quads[i+nquad][0] = quads[i][0]+nv;	quads[i+nquad][1] = quads[i][1]+nv;
		quads[i+nquad][2] = quads[i][2]+nv;	quads[i+nquad][3] = quads[i][3]+nv;
        }
	nv *= 2;
	nquad *= 2;
}*/

void geoframe::write_raw( const char* filename, int meshtype ){

  //    updateBySpan();
    if ( meshtype == SINGLE || meshtype == DOUBLE ) {
                saveTriangle( filename );
        }
        else if ( meshtype == TETRA || meshtype == TETRA2 ){
                saveTetra( filename );
        }
        else if ( meshtype == HEXA ){
                saveHexa( filename );
        }
        else if ( meshtype == QUAD ){
                saveQuad( filename );
        }
}

void geoframe::write_raw( const char* filename)
{
  //    updateBySpan();
    if ( mesh_type == SINGLE || mesh_type == DOUBLE ) {
                saveTriangle( filename );
        }
        else if ( mesh_type == TETRA || mesh_type == TETRA2 ){
                saveTetra( filename );
        }
        else if ( mesh_type == HEXA ){
                saveHexa( filename );
        }
        else if ( mesh_type == QUAD ){
                saveQuad( filename );
        }
}

void geoframe::updateBySpan()
{
	printf("%f %f %f\n",min_x,min_y,min_z);
	int nv= getNVert();
	for (int i=0;i<nv;i++) {
		verts[i][0] = min_x + span[0] * verts[i][0];
		verts[i][1] = min_y + span[1] * verts[i][1];
		verts[i][2] = min_z + span[2] * verts[i][2];
        }
}


void geoframe::saveTriangle( const char* filename ){
	
	FILE* fp;
	
	fp=fopen(filename,"w");
	int v0, v1, v2, nv, ntri, i;
	
	float r, a, b, c, p, sum_area, center;
	sum_area = 0.0;
	center = 64.0f/2.0f;
	
	nv= getNVert();
	ntri= getNTri();

	printf("number of vertices: %d\n",nv);
	printf("number of triangles: %d\n",ntri);
	fprintf(fp,"%d %d\n", nv, ntri);
	for (i=0;i<nv;i++) {
		r = (verts[i][0]-center)*(verts[i][0]-center) +
					(verts[i][1]-center)*(verts[i][1]-center) +
					(verts[i][2]-center)*(verts[i][2]-center);
		fprintf(fp,"%f %f %f\n",verts[i][0],verts[i][1],verts[i][2]);
	}
	for (i=0;i<ntri;i++) {
		v0 = triangles[i][0];
		v1 = triangles[i][1];
		v2 = triangles[i][2];
		r = (verts[v0][0]-center)*(verts[v0][0]-center) +
					(verts[v0][1]-center)*(verts[v0][1]-center) +
					(verts[v0][2]-center)*(verts[v0][2]-center);
		if(sqrt(r) < 17.0) {
			a = sqrt((verts[v1][0]-verts[v0][0])*(verts[v1][0]-verts[v0][0]) +
					 (verts[v1][1]-verts[v0][1])*(verts[v1][1]-verts[v0][1]) +
					 (verts[v1][2]-verts[v0][2])*(verts[v1][2]-verts[v0][2]));
			b = sqrt((verts[v2][0]-verts[v1][0])*(verts[v2][0]-verts[v1][0]) +
					 (verts[v2][1]-verts[v1][1])*(verts[v2][1]-verts[v1][1]) +
					 (verts[v2][2]-verts[v1][2])*(verts[v2][2]-verts[v1][2]));
			c = sqrt((verts[v0][0]-verts[v2][0])*(verts[v0][0]-verts[v2][0]) +
				     (verts[v0][1]-verts[v2][1])*(verts[v0][1]-verts[v2][1]) +
					 (verts[v0][2]-verts[v2][2])*(verts[v0][2]-verts[v2][2]));
			p = (a + b + c)/2.0;
			sum_area += sqrt(p*(p - a)*(p - b)*(p - c));
		}
		fprintf(fp,"%d %d %d\n",triangles[i][0],triangles[i][1],triangles[i][2]);
	}
	//fprintf(fp, "area_14.5 = %f %f\n", sum_area, 4.0*3.1415926*10*10);
	printf("finished\n");
	fclose(fp);
}

void geoframe::saveTetra( const char* filename ){

	FILE* fp;
	fp=fopen(filename,"w");
	int nv, ntri, i;
	float center = (129.0f-1.0f)/2.0f;

	nv=getNVert();
	ntri=getNTri();

	fprintf(fp,"%d %d\n",nv,ntri/4);
	for (i=0;i<nv;i++) {
		//fprintf(fp,"%f %f %f\n",g_frame->verts[i][0],g_frame->verts[i][1],g_frame->verts[i][2]);
		//float r =	(g_frame->verts[i][0]-center)*(g_frame->verts[i][0]-center)*0.6798984*0.6798984 +
		//			(g_frame->verts[i][1]-center)*(g_frame->verts[i][1]-center)*0.6798984*0.6798984 +
		//			(g_frame->verts[i][2]-center)*(g_frame->verts[i][2]-center)*0.6779922*0.6779922;
		//fprintf(fp,"%f %f %f\n",(g_frame->verts[i][0]-center)*0.6798984*2,(g_frame->verts[i][1]-center)*0.6798984*2,
		//	(g_frame->verts[i][2]-center)*0.6779922*2);
		fprintf(fp,"%f %f %f\n",(verts[i][0]-center),(verts[i][1]-center),
			(verts[i][2]-center));
	}
	for (i=0;i<ntri/4;i++) {
		fprintf(fp,"%d %d %d %d\n",triangles[4*i][0],triangles[4*i][1],
			triangles[4*i][2], triangles[4*i+1][2]);
	}
	fclose(fp);
}

void geoframe::saveHexa( const char* filename ){

	FILE* fp;
	fp = fopen(filename,"w");
	int nv, nquad, i;

	nv = getNVert();
	nquad = getNQuad();

	fprintf(fp,"%d %d\n", nv, nquad/6);
	for (i = 0; i < nv; i++) {
		//fprintf(fp,"%f %f %f %f %f %f %d\n", verts[i][0], verts[i][1], verts[i][2],
		//	normals[i][0], normals[i][1], normals[i][2], bound_sign[i]);
		fprintf(fp,"%f %f %f %d\n", verts[i][0], verts[i][1], verts[i][2], bound_sign[i]);
	}
	for (i = 0;i < nquad/6; i++) {
		fprintf(fp,"%d %d %d %d %d %d %d %d\n", quads[6*i][0], quads[6*i][1],
							quads[6*i][2], quads[6*i][3], quads[6*i+1][1],
							quads[6*i+1][0], quads[6*i+1][3], quads[6*i+1][2]);
	}
	fclose(fp);

}

void geoframe::saveQuad( const char* filename ){

	FILE* fp;
	fp = fopen(filename,"w");
	int nv, nquad, i;

	nv=getNVert();
	nquad=getNQuad();

	//fprintf(fp,"OFF\n");	
	// Output 4-node quad element
	fprintf(fp,"%d %d\n", nv, nquad);
	for (i = 0; i < nv; i++) {
		fprintf(fp,"%f %f %f\n", verts[i][0], verts[i][1], verts[i][2]);
	}
	for (i = 0; i < nquad; i++) {
		fprintf(fp,"%d %d %d %d\n", quads[i][0], quads[i][1], quads[i][2], quads[i][3]);
	}
	
	fclose(fp);
}

}
