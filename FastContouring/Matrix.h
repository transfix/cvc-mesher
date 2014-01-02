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

// Matrix.h: interface for the Matrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MATRIX_H__AB7171AD_869D_4D29_A455_5919551F0B67__INCLUDED_)
#define AFX_MATRIX_H__AB7171AD_869D_4D29_A455_5919551F0B67__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

namespace FastContouring
{

class Quaternion;
class Vector;
class Ray;

class Matrix  
{
public:
	Matrix();
	Matrix(
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33
		);
	Matrix(const Quaternion& quat);
	virtual ~Matrix();
	Matrix(const Matrix& copy);
	Matrix& operator=(const Matrix& copy);


	Matrix& set (
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33
		);
	Matrix& set(const Matrix& copy);
	Matrix& reset();
	inline float get(int row, int column) const;
	const float* getMatrix() const;

	Vector operator*(const Vector& vec) const;
	Ray operator*(const Ray& ray) const;
	Matrix operator*(const Matrix& mat) const;
	Matrix& preMultiplication(const Matrix& mat);
	Matrix& postMultiplication(const Matrix& mat);

	Matrix inverse() const;
	Matrix inverseTranspose() const;

	float determinant() const;

	static Matrix rotationX(float angle);
	static Matrix rotationY(float angle);
	static Matrix rotationZ(float angle);
	static Matrix translation(float x, float y, float z);
	static Matrix translation(const Vector& vec);
	static Matrix scale(float x, float y, float z);

protected:
	float m[16];

};

float Matrix::get(int row, int column) const
{
	return m[row + column*4];
}

}

#endif // !defined(AFX_MATRIX_H__AB7171AD_869D_4D29_A455_5919551F0B67__INCLUDED_)
