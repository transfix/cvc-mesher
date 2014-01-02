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

// Quaternion.h: interface for the Quaternion class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_QUATERNION_H__4A5485F3_5ADE_437D_A2C9_6D864A63C23B__INCLUDED_)
#define AFX_QUATERNION_H__4A5485F3_5ADE_437D_A2C9_6D864A63C23B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Tuple.h"

namespace FastContouring
{

class Vector;
class Matrix;
class Ray;

class Quaternion : public Tuple  
{
public:
	Quaternion();
	virtual ~Quaternion();
	Quaternion(const Quaternion& copy);
	Quaternion& operator=(const Quaternion& copy);
	Quaternion(float w, float x, float y, float z);

	Quaternion& set(float w, float x, float y, float z);
	Quaternion& set(float* array);
	Quaternion& set(const Quaternion& copy);

	Quaternion operator+(const Quaternion& quat) const;
	Quaternion operator-(const Quaternion& quat) const;
	Quaternion operator-() const;
	Quaternion operator*(const Quaternion& quat) const;
	Quaternion operator*(float scalar) const;
	Quaternion& operator*=(float scalar);
	Quaternion operator/(float scalar) const;
	Quaternion& operator/=(float scalar);

	Quaternion& preMultiply(const Quaternion& quat);
	Quaternion& postMultiply(const Quaternion& quat);
	Quaternion& rotate(float angle, float x, float y, float z);
	Quaternion& normalize();

	Quaternion conjugate() const;
	Quaternion inverse() const;
	float norm() const;

	Vector applyRotation(const Vector& vec) const;
	Ray applyRotation(const Ray& ray) const;
	Matrix buildMatrix() const;
	Quaternion power(double scalar) const;
	
	static Quaternion log(const Quaternion& a);
	static Quaternion exponent(const Quaternion& a);

	static Quaternion rotation(float angle, float x, float y, float z);
	static Quaternion rotation(float angle, const Vector& axis);
	static Quaternion interpolate(const Quaternion& a, const Quaternion& b, float alpha);
	static Quaternion startCubicInterpolate(const Quaternion& a, const Quaternion& b, const Quaternion& c, float alpha);
	static Quaternion endCubicInterpolate(const Quaternion& a, const Quaternion& b, const Quaternion& c, float alpha);
	static Quaternion cubicInterpolate(const Quaternion& a, const Quaternion& b, const Quaternion& c, const Quaternion& d, float alpha);

protected:
	static Quaternion determineQi(const Quaternion& aiminus1, const Quaternion& ai, const Quaternion& aiplus1);

	explicit Quaternion(const Vector& vec);

};

}

#endif // !defined(AFX_QUATERNION_H__4A5485F3_5ADE_437D_A2C9_6D864A63C23B__INCLUDED_)
