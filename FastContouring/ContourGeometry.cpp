/*
  Copyright 2002-2003 The University of Texas at Austin
  
	Authors: Anthony Thane <thanea@ices.utexas.edu>
	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of Volume Rover.

  Volume Rover is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Volume Rover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with iotree; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

// ContourGeometry.cpp: implementation of the ContourGeometry class.
//
//////////////////////////////////////////////////////////////////////

#include "FastContouring.h"
#include "ContourGeometry.h"
#include <Matrix.h>
#include <Vector.h>
#include <math.h>
//#include "../ipoly/src/ipoly.h"
//#include "../ipoly/src/ipolyutil.h"

namespace FastContouring
{

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ContourGeometry::ContourGeometry()
{
	setDefaults();
}

ContourGeometry::~ContourGeometry()
{
	destroyTriangleBuffers();
	destroyVertexBuffers();
	destroyQuadBuffers();
}

#if 0
void ContourGeometry::addToGeometry(Geometry* geometry, const Matrix& matrix, int& nextVert, int& nextTri)
{
	if (!m_Triangles || !m_Vertex1 || ! m_Normal1) 
		return;

	int nextVertPrivate = nextVert;

	Matrix inverseTranspose = matrix.inverseTranspose();

	doInterpolation(m_Isovalue);
	
	unsigned int c;
	Vector vector;
	Vector point, normal;
	for (c=0; c<m_NumVertices; c++) {
		point[0] = m_Vertex1[c*4+0];
		point[1] = m_Vertex1[c*4+1];
		point[2] = m_Vertex1[c*4+2];
		point[3] = 1.0;
		normal[0] = m_Normal1[c*3+0];
		normal[1] = m_Normal1[c*3+1];
		normal[2] = m_Normal1[c*3+2];
		normal[3] = 0.0;
		// set the data up for the right coord space
		point = matrix*point;
		normal = inverseTranspose*normal;

		geometry->m_TriVerts[nextVertPrivate*3+0] = point[0];
		geometry->m_TriVerts[nextVertPrivate*3+1] = point[1];
		geometry->m_TriVerts[nextVertPrivate*3+2] = point[2];
		geometry->m_TriVertNormals[nextVertPrivate*3+0] = normal[0];
		geometry->m_TriVertNormals[nextVertPrivate*3+1] = normal[1];
		geometry->m_TriVertNormals[nextVertPrivate*3+2] = normal[2];
		if (useColors()) {
			geometry->m_TriVertColors[nextVertPrivate*3+0] = m_Color1[c*3+0];
			geometry->m_TriVertColors[nextVertPrivate*3+1] = m_Color1[c*3+1];
			geometry->m_TriVertColors[nextVertPrivate*3+2] = m_Color1[c*3+2];
		}
		nextVertPrivate++;
	}
	for (c=0; c<m_NumTriangles; c++) {
		geometry->m_Tris[nextTri*3+0] = m_Triangles[c*3+0]+nextVert;
		geometry->m_Tris[nextTri*3+1] = m_Triangles[c*3+1]+nextVert;
		geometry->m_Tris[nextTri*3+2] = m_Triangles[c*3+2]+nextVert;
		nextTri++;
	}
	nextVert=nextVertPrivate;
	
}
#endif
	
void ContourGeometry::addToGeometry(TriSurf& geometry, const Matrix& matrix, int& nextVert, int& nextTri)
{
  if (!m_Triangles || !m_Vertex1 || ! m_Normal1) 
    return;

  int nextVertPrivate = nextVert;

  Matrix inverseTranspose = matrix.inverseTranspose();

  doInterpolation(m_Isovalue);
	
  unsigned int c;
  Vector vector;
  Vector point, normal;
  for (c=0; c<m_NumVertices; c++) {
    point[0] = m_Vertex1[c*4+0];
    point[1] = m_Vertex1[c*4+1];
    point[2] = m_Vertex1[c*4+2];
    point[3] = 1.0;
    normal[0] = m_Normal1[c*3+0];
    normal[1] = m_Normal1[c*3+1];
    normal[2] = m_Normal1[c*3+2];
    normal[3] = 0.0;
    // set the data up for the right coord space
    point = matrix*point;
    normal = inverseTranspose*normal;

    geometry.verts.push_back(point[0]);
    geometry.verts.push_back(point[1]);
    geometry.verts.push_back(point[2]);
    //geometry->m_TriVerts[nextVertPrivate*3+0] = point[0];
    //geometry->m_TriVerts[nextVertPrivate*3+1] = point[1];
    //geometry->m_TriVerts[nextVertPrivate*3+2] = point[2];
    geometry.normals.push_back(normal[0]);
    geometry.normals.push_back(normal[1]);
    geometry.normals.push_back(normal[2]);
    //geometry->m_TriVertNormals[nextVertPrivate*3+0] = normal[0];
    //    geometry->m_TriVertNormals[nextVertPrivate*3+1] = normal[1];
    //    geometry->m_TriVertNormals[nextVertPrivate*3+2] = normal[2];
    if (useColors()) {
      //geometry->m_TriVertColors[nextVertPrivate*3+0] = m_Color1[c*3+0];
      //geometry->m_TriVertColors[nextVertPrivate*3+1] = m_Color1[c*3+1];
      //geometry->m_TriVertColors[nextVertPrivate*3+2] = m_Color1[c*3+2];
      geometry.colors.push_back(m_Color1[c*3+0]);
      geometry.colors.push_back(m_Color1[c*3+1]);
      geometry.colors.push_back(m_Color1[c*3+2]);
    }
    nextVertPrivate++;
  }
  for (c=0; c<m_NumTriangles; c++) {
    geometry.tris.push_back(m_Triangles[c*3+0]+nextVert);
    geometry.tris.push_back(m_Triangles[c*3+1]+nextVert);
    geometry.tris.push_back(m_Triangles[c*3+2]+nextVert);
    //geometry->m_Tris[nextTri*3+0] = m_Triangles[c*3+0]+nextVert;
    //geometry->m_Tris[nextTri*3+1] = m_Triangles[c*3+1]+nextVert;
    //geometry->m_Tris[nextTri*3+2] = m_Triangles[c*3+2]+nextVert;
    nextTri++;
  }
  nextVert=nextVertPrivate;
	
}

bool ContourGeometry::addQuadVertex(float vx, float vy, float vz, float nx, float ny, float nz, float cx, float cy, float cz)
{
	if (m_NumVertices==m_NumVerticesAllocated) {
		if (!doubleVertexBuffers()) {
			return false;
		}
	}

	m_Vertex1[m_NumVertices*3+0] = vx;
	m_Vertex1[m_NumVertices*3+1] = vy;
	m_Vertex1[m_NumVertices*3+2] = vz;
	m_Normal1[m_NumVertices*3+0] = nx;
	m_Normal1[m_NumVertices*3+1] = ny;
	m_Normal1[m_NumVertices*3+2] = nz;
	m_Color1[m_NumVertices*3+0] = cx;
	m_Color1[m_NumVertices*3+1] = cy;
	m_Color1[m_NumVertices*3+2] = cz;

	m_NumVertices++;
	return true;
}

bool ContourGeometry::addQuad(unsigned int v1, unsigned int v2, unsigned int v3, unsigned int v4)
{
	// add a single quad, doubling the array if necessary
	if (m_NumQuads==m_NumQuadsAllocated) {
		if (!doubleQuadBuffers()) {
			return false;
		}
	}

	m_Quads[m_NumQuads*4+0] = v1;
	m_Quads[m_NumQuads*4+1] = v2;
	m_Quads[m_NumQuads*4+2] = v3;
	m_Quads[m_NumQuads*4+3] = v4;

	m_NumQuads++;
	return true;
}

int ContourGeometry::getNumVerts()
{
	return m_NumVertices;
}

int ContourGeometry::getNumTris()
{
	return m_NumTriangles;
}

void ContourGeometry::setDefaults()
{
	// set default values for all member variables
	m_Normal1 = 0;
	m_Normal2 = 0;
	m_Vertex1 = 0;
	m_Vertex2 = 0;
	m_Color1 = 0;
	m_Color2 = 0;
	m_NumVertices = 0;
	m_NumVerticesAllocated = 0;

	m_Triangles = 0;
	m_NumTriangles = 0;
	m_NumTrianglesAllocated = 0;

	m_Quads = 0;
	m_NumQuads = 0;
	m_NumQuadsAllocated = 0;

	m_Isovalue = 0.0f;

	//	m_HardwareAccelerated = false;
	m_InterpolationDone = false;
	m_UseColors = false;
	//	m_WireframeRender = false;
	//	m_SurfWithWire = false;
}

void ContourGeometry::setUseColors(bool useColors)
{
	m_UseColors = useColors;
}

void ContourGeometry::setIsovalue(float isovalue)
{
	m_Isovalue = isovalue;
}

void ContourGeometry::setSingleColor(float R, float G, float B)
{
	int i;

	// reassign all the vertex colors
	for (i=0; i < (int)m_NumVertices; i++) {
		m_Color1[i*3+0] = R;
		m_Color1[i*3+1] = G;
		m_Color1[i*3+2] = B;
	}
}

bool ContourGeometry::useColors()
{
	return m_UseColors;
}

bool ContourGeometry::allocateVertexBuffers(unsigned int initialSize)
{
	// if we already have a buffer big enough, just return true and use that buffer
	if (initialSize>m_NumVerticesAllocated) {
		destroyVertexBuffers();
		m_NumVertices = 0;
		m_InterpolationDone = false;
		return forceAllocateVertexBuffers(initialSize);
	}
	else {
		m_NumVertices = 0;
		m_InterpolationDone = false;
		return true;
	}
}

bool ContourGeometry::forceAllocateVertexBuffers(unsigned int size)
{
	// allocate the buffers without checking to see if we already have a big enough buffer
	m_Vertex1 = new float[size*4];
	m_Vertex2 = new float[size*4];
	m_Normal1 = new float[size*3];
	m_Normal2 = new float[size*3];
	m_Color1 = new float[size*3];
	m_Color2 = new float[size*3];

	if (m_Vertex1 && m_Vertex2 && m_Normal1 && m_Normal2 && m_Color1 && m_Color2) {
		m_NumVerticesAllocated = size;
		return true;
	}
	else {
		destroyVertexBuffers();
		return false;
	}
}

bool ContourGeometry::doubleVertexBuffers()
{
	// double the size of the vertex buffer and copy the vertices to the new buffer
	if (
		!doubleFloatArray(m_Vertex1, m_NumVerticesAllocated*4) ||
		!doubleFloatArray(m_Vertex2, m_NumVerticesAllocated*4) ||
		!doubleFloatArray(m_Normal1, m_NumVerticesAllocated*3) ||
		!doubleFloatArray(m_Normal2, m_NumVerticesAllocated*3) ||
		!doubleFloatArray(m_Color1, m_NumVerticesAllocated*3) ||
		!doubleFloatArray(m_Color2, m_NumVerticesAllocated*3)) {
		destroyVertexBuffers();
		return false;
	}
	else { // success
		m_NumVerticesAllocated*=2;
		return true;
	}
}

void ContourGeometry::destroyVertexBuffers()
{
	// free vertex buffer memory
	delete [] m_Vertex1;
	m_Vertex1 = 0;

	delete [] m_Vertex2;
	m_Vertex2 = 0;

	delete [] m_Normal1;
	m_Normal1 = 0;

	delete [] m_Normal2;
	m_Normal2 = 0;

	delete [] m_Color1;
	m_Color1 = 0;

	delete [] m_Color2;
	m_Color2 = 0;

	m_NumVerticesAllocated = 0;
	m_NumVertices = 0;
}

bool ContourGeometry::allocateTriangleBuffers(unsigned int initialSize)
{
	// if we already have a buffer big enough, just return true and use that buffer
	if (initialSize>m_NumTrianglesAllocated) {
		destroyTriangleBuffers();
		m_NumTriangles = 0;
		m_InterpolationDone = false;
		return forceAllocateTriangleBuffers(initialSize);
	}
	else {
		m_NumTriangles = 0;
		m_InterpolationDone = false;
		return true;
	}
}

bool ContourGeometry::forceAllocateTriangleBuffers(unsigned int size)
{
	// allocate the buffers without checking to see if we already have a big enough buffer
	m_Triangles = new unsigned int[size*3];

	if (m_Triangles) {
		m_NumTrianglesAllocated = size;
		return true;
	}
	else {
		destroyTriangleBuffers();
		return false;
	}
}

bool ContourGeometry::doubleTriangleBuffers()
{
	// double the size of the triangle buffer and copy the triangles to the new buffer
	if (
		!doubleIntArray(m_Triangles, m_NumTrianglesAllocated*3)) {
		destroyTriangleBuffers();
		return false;
	}
	else { // success
		m_NumTrianglesAllocated*=2;
		return true;
	}
}

void ContourGeometry::destroyTriangleBuffers()
{
	// free triangle buffer memory
	delete [] m_Triangles;
	m_Triangles = 0;

	m_NumTrianglesAllocated = 0;
	m_NumTriangles = 0;
}

bool ContourGeometry::allocateQuadBuffers(unsigned int initialSize)
{
	// if we already have a buffer big enough, just return true and use that buffer
	if (initialSize>m_NumQuadsAllocated) {
		destroyQuadBuffers();
		m_NumQuads = 0;
		//m_InterpolationDone = false;
		return forceAllocateQuadBuffers(initialSize);
	}
	else {
		m_NumQuads= 0;
		//m_InterpolationDone = false;
		return true;
	}
}

bool ContourGeometry::forceAllocateQuadBuffers(unsigned int size)
{
	// allocate the buffers without checking to see if we already have a big enough buffer
	m_Quads = new unsigned int[size*4];

	if (m_Quads) {
		m_NumQuadsAllocated = size;
		return true;
	}
	else {
		destroyQuadBuffers();
		return false;
	}
}

bool ContourGeometry::doubleQuadBuffers()
{
	// double the size of the quad buffer and copy the quads to the new buffer
	if (!doubleIntArray(m_Quads, m_NumQuadsAllocated*4)) {
		destroyQuadBuffers();
		return false;
	}
	else { // success
		m_NumQuadsAllocated*=2;
		return true;
	}
}

void ContourGeometry::destroyQuadBuffers()
{
	// free quad buffer memory
	delete [] m_Quads;
	m_Quads = 0;

	m_NumQuadsAllocated = 0;
	m_NumQuads = 0;
}

bool ContourGeometry::doubleFloatArray(float*& array, unsigned int size)
{
	// double a single float array, copy over the values, and delete the old array
	float* newArray = new float[size*2];
	if (!newArray) {
		return false;
	}
	unsigned int c;
	for (c=0; c<size; c++) {
		newArray[c] = array[c];
	}
	delete [] array;
	array = newArray;
	return true;
}

bool ContourGeometry::doubleIntArray(unsigned int*& array, unsigned int size)
{
	// double a single uint array, copy over the values, and delete the old array
	unsigned int* newArray = new unsigned int[size*2];
	if (!newArray) {
		return false;
	}
	unsigned int c;
	for (c=0; c<size; c++) {
		newArray[c] = array[c];
	}
	delete [] array;
	array = newArray;
	return true;
}


bool ContourGeometry::addTriangle(unsigned int v1, unsigned int v2, unsigned int v3)
{
	// add a single triangle, doubling the array if necessary
	if (m_NumTriangles==m_NumTrianglesAllocated) {
		if (!doubleTriangleBuffers()) {
			return false;
		}
	}

	m_Triangles[m_NumTriangles*3+0] = v1;
	m_Triangles[m_NumTriangles*3+1] = v2;
	m_Triangles[m_NumTriangles*3+2] = v3;

	m_NumTriangles++;
	return true;
}

void ContourGeometry::doInterpolation(float isovalue)
{
	if (!m_InterpolationDone) {
		unsigned int c;
		if (m_UseColors) {
			for (c=0; c<m_NumVertices; c++) {
				interpArray(m_Vertex1+c*4, m_Vertex2+c*4, m_Normal1+c*3, m_Normal2+c*3, m_Vertex1[c*4+3], m_Vertex2[c*4+3], isovalue);
			}
		}
		else {
			for (c=0; c<m_NumVertices; c++) {
				interpArray(m_Vertex1+c*4, m_Vertex2+c*4, m_Normal1+c*3, m_Normal2+c*3, m_Vertex1[c*4+3], m_Vertex2[c*4+3], isovalue);
			}
		}
		m_InterpolationDone = true;
	}
}

static void cross(float* dest, const float* v1, const float* v2)
{

	dest[0] = v1[1]*v2[2] - v1[2]*v2[1];
	dest[1] = v1[2]*v2[0] - v1[0]*v2[2];
	dest[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

static void normalize(float* v)
{
	float len = (float)sqrt(
		v[0] * v[0] +
		v[1] * v[1] +
		v[2] * v[2]);
	if (len!=0.0) {
		v[0]/=len; //*biggestDim;
		v[1]/=len; //*biggestDim;
		v[2]/=len; //*biggestDim;
	}
	else {
		v[0] = 1.0;
	}
}

static void add(float* dest, const float* v)
{
	dest[0] += v[0];
	dest[1] += v[1];
	dest[2] += v[2];
}

static void set(float* dest, const float* v)
{
	dest[0] = v[0];
	dest[1] = v[1];
	dest[2] = v[2];
}

void ContourGeometry::CalculateQuadSmoothNormals()
{
	unsigned int c, v0, v1, v2, v3;
	float normal[3];

	float zero[3] = {0.0f, 0.0f, 0.0f};

	for (c=0; c<m_NumVertices; c++) {
		set(m_Normal1+c*3, zero);
	}

		
	// for each Quadangle
	for (c=0; c<m_NumQuads; c++) {
		v0 = m_Quads[c*4+0];
		v1 = m_Quads[c*4+1];
		v2 = m_Quads[c*4+2];
		v3 = m_Quads[c*4+3];
		CalculateQuadNormal(normal, v0, v1, v3);
		add(m_Normal1+v0*3, normal);
		add(m_Normal1+v1*3, normal);
		add(m_Normal1+v3*3, normal);
		CalculateQuadNormal(normal, v2, v3, v1);
		add(m_Normal1+v2*3, normal);
		add(m_Normal1+v3*3, normal);
		add(m_Normal1+v1*3, normal);
	}
	
	// normalize the vectors
	for (c=0; c<m_NumVertices; c++) {
		normalize(m_Normal1+c*3);
	}
}

void ContourGeometry::CalculateQuadNormal(float* norm, unsigned int v0, unsigned int v1, unsigned int v2)
{
	float vec1[3], vec2[3];
	vec1[0] = vec2[0] = -m_Vertex1[v0*3+0];
	vec1[1] = vec2[1] = -m_Vertex1[v0*3+1];
	vec1[2] = vec2[2] = -m_Vertex1[v0*3+2];
	vec1[0] += m_Vertex1[v1*3+0];
	vec1[1] += m_Vertex1[v1*3+1];
	vec1[2] += m_Vertex1[v1*3+2];
	
	vec2[0] += m_Vertex1[v2*3+0];
	vec2[1] += m_Vertex1[v2*3+1];
	vec2[2] += m_Vertex1[v2*3+2];


	cross(norm, vec1, vec2);
}

}
