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

// Vector.h: interface for the Vector class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VECTOR_H__C73C6BDB_2D75_4770_B86A_E3329C07A817__INCLUDED_)
#define AFX_VECTOR_H__C73C6BDB_2D75_4770_B86A_E3329C07A817__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Tuple.h"

namespace FastContouring
{

class Vector : public Tuple  
{
public:
	Vector(float x, float y, float z, float w);
	Vector(float* array);
	Vector();
	virtual ~Vector();
	Vector(const Vector& copy);
	Vector& operator=(const Vector& copy);

	Vector& set(float x, float y, float z, float w);
	Vector& set(float* array);
	Vector& set(const Vector& copy);

	Vector cross(const Vector& vec) const;
	Vector& crossEquals(const Vector& vec);
	float dot(const Vector& vec) const;

	Vector operator+(const Vector vec) const;
	Vector& operator+=(const Vector vec);
	Vector operator-(const Vector vec) const;
	Vector& operator-=(const Vector vec);

	Vector operator*(float scalar) const;
	Vector& operator*=(float scalar);

	Vector operator-() const;

	Vector& normalize();
	float norm();

	bool isBad();

	static Vector badVector();
	virtual Vector* clone() const;

	static Vector interpolate(const Vector& a, const Vector& b, float alpha);
	static Vector cubicInterpolate(const Vector& a, const Vector& b, const Vector& c, const Vector& d, float alpha);

};

}

#endif // !defined(AFX_VECTOR_H__C73C6BDB_2D75_4770_B86A_E3329C07A817__INCLUDED_)
