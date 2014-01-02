/*
  Copyright 2002-2003 The University of Texas at Austin
  
    Authors: Anthony Thane <thanea@ices.utexas.edu>
             Vinay Siddavanahalli <skvinay@cs.utexas.edu>
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

// Matrix.cpp: implementation of the Matrix class.
//
//////////////////////////////////////////////////////////////////////

#include "Matrix.h"
#include "Vector.h"
#include "Ray.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

namespace FastContouring
{

Matrix::Matrix()
{
	set(1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix::Matrix(
	float m00, float m01, float m02, float m03,
	float m10, float m11, float m12, float m13,
	float m20, float m21, float m22, float m23,
	float m30, float m31, float m32, float m33
	)
{
	m[0]=m00; m[1]=m10; m[2]=m20; m[3]=m30;
	m[4]=m01; m[5]=m11; m[6]=m21; m[7]=m31;
	m[8]=m02; m[9]=m12; m[10]=m22; m[11]=m32;
	m[12]=m03; m[13]=m13; m[14]=m23; m[15]=m33;
}

Matrix::~Matrix()
{

}

Matrix::Matrix(const Matrix& copy)
{
	set(copy);
}

Matrix& Matrix::operator=(const Matrix& copy)
{
	return set(copy);
}


Matrix& Matrix::set (
	float m00, float m01, float m02, float m03,
	float m10, float m11, float m12, float m13,
	float m20, float m21, float m22, float m23,
	float m30, float m31, float m32, float m33
	)
{
	m[0]=m00; m[1]=m10; m[2]=m20; m[3]=m30;
	m[4]=m01; m[5]=m11; m[6]=m21; m[7]=m31;
	m[8]=m02; m[9]=m12; m[10]=m22; m[11]=m32;
	m[12]=m03; m[13]=m13; m[14]=m23; m[15]=m33;
	return *this;
}

Matrix& Matrix::set(const Matrix& copy)
{
	if (this!=&copy) {
		set(
			copy.m[0], copy.m[4], copy.m[8], copy.m[12],
			copy.m[1], copy.m[5], copy.m[9], copy.m[13],
			copy.m[2], copy.m[6], copy.m[10], copy.m[14],
			copy.m[3], copy.m[7], copy.m[11], copy.m[15]
			);
	}
	return *this;
}

Matrix& Matrix::reset()
{
	set(
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0);
	return *this;
}

const float* Matrix::getMatrix() const
{
	return m;
}

Vector Matrix::operator*(const Vector& vec) const
{
	return Vector(
		m[0]*vec[0]+m[4]*vec[1]+m[8]*vec[2]+m[12]*vec[3],
		m[1]*vec[0]+m[5]*vec[1]+m[9]*vec[2]+m[13]*vec[3],
		m[2]*vec[0]+m[6]*vec[1]+m[10]*vec[2]+m[14]*vec[3],
		m[3]*vec[0]+m[7]*vec[1]+m[11]*vec[2]+m[15]*vec[3]
		);
}

Ray Matrix::operator*(const Ray& ray) const
{
	return Ray((*this)*ray.m_Origin, (*this)*ray.m_Dir);
}

Matrix Matrix::operator*(const Matrix& mat) const
{
	return Matrix(
		m[0]*mat.m[0]+m[4]*mat.m[1]+m[8]*mat.m[2]+m[12]*mat.m[3],
		m[0]*mat.m[4]+m[4]*mat.m[5]+m[8]*mat.m[6]+m[12]*mat.m[7],
		m[0]*mat.m[8]+m[4]*mat.m[9]+m[8]*mat.m[10]+m[12]*mat.m[11],
		m[0]*mat.m[12]+m[4]*mat.m[13]+m[8]*mat.m[14]+m[12]*mat.m[15],

		m[1]*mat.m[0]+m[5]*mat.m[1]+m[9]*mat.m[2]+m[13]*mat.m[3],
		m[1]*mat.m[4]+m[5]*mat.m[5]+m[9]*mat.m[6]+m[13]*mat.m[7],
		m[1]*mat.m[8]+m[5]*mat.m[9]+m[9]*mat.m[10]+m[13]*mat.m[11],
		m[1]*mat.m[12]+m[5]*mat.m[13]+m[9]*mat.m[14]+m[13]*mat.m[15],

		m[2]*mat.m[0]+m[6]*mat.m[1]+m[10]*mat.m[2]+m[14]*mat.m[3],
		m[2]*mat.m[4]+m[6]*mat.m[5]+m[10]*mat.m[6]+m[14]*mat.m[7],
		m[2]*mat.m[8]+m[6]*mat.m[9]+m[10]*mat.m[10]+m[14]*mat.m[11],
		m[2]*mat.m[12]+m[6]*mat.m[13]+m[10]*mat.m[14]+m[14]*mat.m[15],

		m[3]*mat.m[0]+m[7]*mat.m[1]+m[11]*mat.m[2]+m[15]*mat.m[3],
		m[3]*mat.m[4]+m[7]*mat.m[5]+m[11]*mat.m[6]+m[15]*mat.m[7],
		m[3]*mat.m[8]+m[7]*mat.m[9]+m[11]*mat.m[10]+m[15]*mat.m[11],
		m[3]*mat.m[12]+m[7]*mat.m[13]+m[11]*mat.m[14]+m[15]*mat.m[15]
		);
}

Matrix& Matrix::preMultiplication(const Matrix& mat)
{
	return set(mat*(*this));
}

Matrix& Matrix::postMultiplication(const Matrix& mat)
{
	return set((*this)*mat);
}

Matrix Matrix::inverse() const
{
	Matrix ret;
	float det = determinant();
	if (det!=0.0) {
		ret.set(
			(get(1, 2)*get(2, 3)*get(3, 1) - get(1, 3)*get(2, 2)*get(3, 1) + get(1, 3)*get(2, 1)*get(3, 2) - get(1, 1)*get(2, 3)*get(3, 2) - get(1, 2)*get(2, 1)*get(3, 3) + get(1, 1)*get(2, 2)*get(3, 3))/det,
			(get(0, 3)*get(2, 2)*get(3, 1) - get(0, 2)*get(2, 3)*get(3, 1) - get(0, 3)*get(2, 1)*get(3, 2) + get(0, 1)*get(2, 3)*get(3, 2) + get(0, 2)*get(2, 1)*get(3, 3) - get(0, 1)*get(2, 2)*get(3, 3))/det,
			(get(0, 2)*get(1, 3)*get(3, 1) - get(0, 3)*get(1, 2)*get(3, 1) + get(0, 3)*get(1, 1)*get(3, 2) - get(0, 1)*get(1, 3)*get(3, 2) - get(0, 2)*get(1, 1)*get(3, 3) + get(0, 1)*get(1, 2)*get(3, 3))/det,
			(get(0, 3)*get(1, 2)*get(2, 1) - get(0, 2)*get(1, 3)*get(2, 1) - get(0, 3)*get(1, 1)*get(2, 2) + get(0, 1)*get(1, 3)*get(2, 2) + get(0, 2)*get(1, 1)*get(2, 3) - get(0, 1)*get(1, 2)*get(2, 3))/det,

			(get(1, 3)*get(2, 2)*get(3, 0) - get(1, 2)*get(2, 3)*get(3, 0) - get(1, 3)*get(2, 0)*get(3, 2) + get(1, 0)*get(2, 3)*get(3, 2) + get(1, 2)*get(2, 0)*get(3, 3) - get(1, 0)*get(2, 2)*get(3, 3))/det,
			(get(0, 2)*get(2, 3)*get(3, 0) - get(0, 3)*get(2, 2)*get(3, 0) + get(0, 3)*get(2, 0)*get(3, 2) - get(0, 0)*get(2, 3)*get(3, 2) - get(0, 2)*get(2, 0)*get(3, 3) + get(0, 0)*get(2, 2)*get(3, 3))/det,
			(get(0, 3)*get(1, 2)*get(3, 0) - get(0, 2)*get(1, 3)*get(3, 0) - get(0, 3)*get(1, 0)*get(3, 2) + get(0, 0)*get(1, 3)*get(3, 2) + get(0, 2)*get(1, 0)*get(3, 3) - get(0, 0)*get(1, 2)*get(3, 3))/det,
			(get(0, 2)*get(1, 3)*get(2, 0) - get(0, 3)*get(1, 2)*get(2, 0) + get(0, 3)*get(1, 0)*get(2, 2) - get(0, 0)*get(1, 3)*get(2, 2) - get(0, 2)*get(1, 0)*get(2, 3) + get(0, 0)*get(1, 2)*get(2, 3))/det,

			(get(1, 1)*get(2, 3)*get(3, 0) - get(1, 3)*get(2, 1)*get(3, 0) + get(1, 3)*get(2, 0)*get(3, 1) - get(1, 0)*get(2, 3)*get(3, 1) - get(1, 1)*get(2, 0)*get(3, 3) + get(1, 0)*get(2, 1)*get(3, 3))/det,
			(get(0, 3)*get(2, 1)*get(3, 0) - get(0, 1)*get(2, 3)*get(3, 0) - get(0, 3)*get(2, 0)*get(3, 1) + get(0, 0)*get(2, 3)*get(3, 1) + get(0, 1)*get(2, 0)*get(3, 3) - get(0, 0)*get(2, 1)*get(3, 3))/det,
			(get(0, 1)*get(1, 3)*get(3, 0) - get(0, 3)*get(1, 1)*get(3, 0) + get(0, 3)*get(1, 0)*get(3, 1) - get(0, 0)*get(1, 3)*get(3, 1) - get(0, 1)*get(1, 0)*get(3, 3) + get(0, 0)*get(1, 1)*get(3, 3))/det,
			(get(0, 3)*get(1, 1)*get(2, 0) - get(0, 1)*get(1, 3)*get(2, 0) - get(0, 3)*get(1, 0)*get(2, 1) + get(0, 0)*get(1, 3)*get(2, 1) + get(0, 1)*get(1, 0)*get(2, 3) - get(0, 0)*get(1, 1)*get(2, 3))/det,

			(get(1, 2)*get(2, 1)*get(3, 0) - get(1, 1)*get(2, 2)*get(3, 0) - get(1, 2)*get(2, 0)*get(3, 1) + get(1, 0)*get(2, 2)*get(3, 1) + get(1, 1)*get(2, 0)*get(3, 2) - get(1, 0)*get(2, 1)*get(3, 2))/det,
			(get(0, 1)*get(2, 2)*get(3, 0) - get(0, 2)*get(2, 1)*get(3, 0) + get(0, 2)*get(2, 0)*get(3, 1) - get(0, 0)*get(2, 2)*get(3, 1) - get(0, 1)*get(2, 0)*get(3, 2) + get(0, 0)*get(2, 1)*get(3, 2))/det,
			(get(0, 2)*get(1, 1)*get(3, 0) - get(0, 1)*get(1, 2)*get(3, 0) - get(0, 2)*get(1, 0)*get(3, 1) + get(0, 0)*get(1, 2)*get(3, 1) + get(0, 1)*get(1, 0)*get(3, 2) - get(0, 0)*get(1, 1)*get(3, 2))/det,
			(get(0, 1)*get(1, 2)*get(2, 0) - get(0, 2)*get(1, 1)*get(2, 0) + get(0, 2)*get(1, 0)*get(2, 1) - get(0, 0)*get(1, 2)*get(2, 1) - get(0, 1)*get(1, 0)*get(2, 2) + get(0, 0)*get(1, 1)*get(2, 2))/det
		);
	}
	// if det==0.0, ret will be identity
	return ret;
}

Matrix Matrix::inverseTranspose() const
{
	Matrix ret;
	float det = determinant();
	if (det!=0.0) {
		ret.set(
			(get(1, 2)*get(2, 3)*get(3, 1) - get(1, 3)*get(2, 2)*get(3, 1) + get(1, 3)*get(2, 1)*get(3, 2) - get(1, 1)*get(2, 3)*get(3, 2) - get(1, 2)*get(2, 1)*get(3, 3) + get(1, 1)*get(2, 2)*get(3, 3))/det,
			(get(1, 3)*get(2, 2)*get(3, 0) - get(1, 2)*get(2, 3)*get(3, 0) - get(1, 3)*get(2, 0)*get(3, 2) + get(1, 0)*get(2, 3)*get(3, 2) + get(1, 2)*get(2, 0)*get(3, 3) - get(1, 0)*get(2, 2)*get(3, 3))/det,
			(get(1, 1)*get(2, 3)*get(3, 0) - get(1, 3)*get(2, 1)*get(3, 0) + get(1, 3)*get(2, 0)*get(3, 1) - get(1, 0)*get(2, 3)*get(3, 1) - get(1, 1)*get(2, 0)*get(3, 3) + get(1, 0)*get(2, 1)*get(3, 3))/det,
			(get(1, 2)*get(2, 1)*get(3, 0) - get(1, 1)*get(2, 2)*get(3, 0) - get(1, 2)*get(2, 0)*get(3, 1) + get(1, 0)*get(2, 2)*get(3, 1) + get(1, 1)*get(2, 0)*get(3, 2) - get(1, 0)*get(2, 1)*get(3, 2))/det,

			(get(0, 3)*get(2, 2)*get(3, 1) - get(0, 2)*get(2, 3)*get(3, 1) - get(0, 3)*get(2, 1)*get(3, 2) + get(0, 1)*get(2, 3)*get(3, 2) + get(0, 2)*get(2, 1)*get(3, 3) - get(0, 1)*get(2, 2)*get(3, 3))/det,
			(get(0, 2)*get(2, 3)*get(3, 0) - get(0, 3)*get(2, 2)*get(3, 0) + get(0, 3)*get(2, 0)*get(3, 2) - get(0, 0)*get(2, 3)*get(3, 2) - get(0, 2)*get(2, 0)*get(3, 3) + get(0, 0)*get(2, 2)*get(3, 3))/det,
			(get(0, 3)*get(2, 1)*get(3, 0) - get(0, 1)*get(2, 3)*get(3, 0) - get(0, 3)*get(2, 0)*get(3, 1) + get(0, 0)*get(2, 3)*get(3, 1) + get(0, 1)*get(2, 0)*get(3, 3) - get(0, 0)*get(2, 1)*get(3, 3))/det,
			(get(0, 1)*get(2, 2)*get(3, 0) - get(0, 2)*get(2, 1)*get(3, 0) + get(0, 2)*get(2, 0)*get(3, 1) - get(0, 0)*get(2, 2)*get(3, 1) - get(0, 1)*get(2, 0)*get(3, 2) + get(0, 0)*get(2, 1)*get(3, 2))/det,

			(get(0, 2)*get(1, 3)*get(3, 1) - get(0, 3)*get(1, 2)*get(3, 1) + get(0, 3)*get(1, 1)*get(3, 2) - get(0, 1)*get(1, 3)*get(3, 2) - get(0, 2)*get(1, 1)*get(3, 3) + get(0, 1)*get(1, 2)*get(3, 3))/det,
			(get(0, 3)*get(1, 2)*get(3, 0) - get(0, 2)*get(1, 3)*get(3, 0) - get(0, 3)*get(1, 0)*get(3, 2) + get(0, 0)*get(1, 3)*get(3, 2) + get(0, 2)*get(1, 0)*get(3, 3) - get(0, 0)*get(1, 2)*get(3, 3))/det,
			(get(0, 1)*get(1, 3)*get(3, 0) - get(0, 3)*get(1, 1)*get(3, 0) + get(0, 3)*get(1, 0)*get(3, 1) - get(0, 0)*get(1, 3)*get(3, 1) - get(0, 1)*get(1, 0)*get(3, 3) + get(0, 0)*get(1, 1)*get(3, 3))/det,
			(get(0, 2)*get(1, 1)*get(3, 0) - get(0, 1)*get(1, 2)*get(3, 0) - get(0, 2)*get(1, 0)*get(3, 1) + get(0, 0)*get(1, 2)*get(3, 1) + get(0, 1)*get(1, 0)*get(3, 2) - get(0, 0)*get(1, 1)*get(3, 2))/det,
			
			(get(0, 3)*get(1, 2)*get(2, 1) - get(0, 2)*get(1, 3)*get(2, 1) - get(0, 3)*get(1, 1)*get(2, 2) + get(0, 1)*get(1, 3)*get(2, 2) + get(0, 2)*get(1, 1)*get(2, 3) - get(0, 1)*get(1, 2)*get(2, 3))/det,
			(get(0, 2)*get(1, 3)*get(2, 0) - get(0, 3)*get(1, 2)*get(2, 0) + get(0, 3)*get(1, 0)*get(2, 2) - get(0, 0)*get(1, 3)*get(2, 2) - get(0, 2)*get(1, 0)*get(2, 3) + get(0, 0)*get(1, 2)*get(2, 3))/det,
			(get(0, 3)*get(1, 1)*get(2, 0) - get(0, 1)*get(1, 3)*get(2, 0) - get(0, 3)*get(1, 0)*get(2, 1) + get(0, 0)*get(1, 3)*get(2, 1) + get(0, 1)*get(1, 0)*get(2, 3) - get(0, 0)*get(1, 1)*get(2, 3))/det,
			(get(0, 1)*get(1, 2)*get(2, 0) - get(0, 2)*get(1, 1)*get(2, 0) + get(0, 2)*get(1, 0)*get(2, 1) - get(0, 0)*get(1, 2)*get(2, 1) - get(0, 1)*get(1, 0)*get(2, 2) + get(0, 0)*get(1, 1)*get(2, 2))/det
		);
	}
	// if det==0.0, ret will be identity
	return ret;
}

float Matrix::determinant() const
{
	double ret;
	ret = 
		get(0, 3) * get(1, 2) * get(2, 1) * get(3, 0)-get(0, 2) * get(1, 3) * get(2, 1) * get(3, 0)-get(0, 3) * get(1, 1) * get(2, 2) * get(3, 0)+get(0, 1) * get(1, 3) * get(2, 2) * get(3, 0)+
		get(0, 2) * get(1, 1) * get(2, 3) * get(3, 0)-get(0, 1) * get(1, 2) * get(2, 3) * get(3, 0)-get(0, 3) * get(1, 2) * get(2, 0) * get(3, 1)+get(0, 2) * get(1, 3) * get(2, 0) * get(3, 1)+
		get(0, 3) * get(1, 0) * get(2, 2) * get(3, 1)-get(0, 0) * get(1, 3) * get(2, 2) * get(3, 1)-get(0, 2) * get(1, 0) * get(2, 3) * get(3, 1)+get(0, 0) * get(1, 2) * get(2, 3) * get(3, 1)+
		get(0, 3) * get(1, 1) * get(2, 0) * get(3, 2)-get(0, 1) * get(1, 3) * get(2, 0) * get(3, 2)-get(0, 3) * get(1, 0) * get(2, 1) * get(3, 2)+get(0, 0) * get(1, 3) * get(2, 1) * get(3, 2)+
		get(0, 1) * get(1, 0) * get(2, 3) * get(3, 2)-get(0, 0) * get(1, 1) * get(2, 3) * get(3, 2)-get(0, 2) * get(1, 1) * get(2, 0) * get(3, 3)+get(0, 1) * get(1, 2) * get(2, 0) * get(3, 3)+
		get(0, 2) * get(1, 0) * get(2, 1) * get(3, 3)-get(0, 0) * get(1, 2) * get(2, 1) * get(3, 3)-get(0, 1) * get(1, 0) * get(2, 2) * get(3, 3)+get(0, 0) * get(1, 1) * get(2, 2) * get(3, 3);
	return (float)ret;
}

Matrix Matrix::rotationX(float angle)
{
	float ca = (float)cos(angle);
	float sa = (float)sin(angle);
	return Matrix(1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, ca, sa, 0.0f,
		0.0f, -sa, ca, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix Matrix::rotationY(float angle)
{
	float ca = (float)cos(angle);
	float sa = (float)sin(angle);
	return Matrix(ca, 0.0f, -sa, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		sa, 0.0f, ca, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix Matrix::rotationZ(float angle)
{
	float ca = (float)cos(angle);
	float sa = (float)sin(angle);
	return Matrix(ca, sa, 0.0f, 0.0f,
		-sa, ca, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix Matrix::translation(float x, float y, float z)
{
	return Matrix(1.0f, 0.0f, 0.0f, x,
		0.0f, 1.0f, 0.0f, y,
		0.0f, 0.0f, 1.0f, z,
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix Matrix::translation(const Vector& vec)
{
	return Matrix(1.0f, 0.0f, 0.0f, vec[0],
		0.0f, 1.0f, 0.0f, vec[1],
		0.0f, 0.0f, 1.0f, vec[2],
		0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix Matrix::scale(float x, float y, float z)
{
	return Matrix(x, 0.0f, 0.0f, 0.0f,
		0.0f, y, 0.0f, 0.0f,
		0.0f, 0.0f, z, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
}

}
