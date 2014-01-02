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

// Vector.cpp: implementation of the Vector class.
//
//////////////////////////////////////////////////////////////////////

#include "Vector.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

namespace FastContouring
{

const float EPS = 0.00001f;

Vector::Vector(float x, float y, float z, float w) : Tuple(x,y,z,w)
{
}

Vector::Vector() : Tuple()
{
}

Vector::Vector(float* array)
{
	set(array);
}

Vector::~Vector()
{

}

Vector::Vector(const Vector& copy): Tuple(copy)
{
}

Vector& Vector::operator=(const Vector& copy)
{
	if (this!=&copy) {
		set(copy);
	}
	return *this;
}


Vector& Vector::set(float x, float y, float z, float w)
{
	Tuple::set(x,y,z,w);
	return *this;
}

Vector& Vector::set(float* array)
{
	Tuple::set(array);
	return *this;
}

Vector& Vector::set(const Vector& copy)
{
	Tuple::set(copy);
	return *this;
}

Vector Vector::cross(const Vector& vec) const
{
	return Vector(
		p[1]*vec[2] - p[2]*vec[1],
		p[2]*vec[0] - p[0]*vec[2],
		p[0]*vec[1] - p[1]*vec[0],		
		0.0f		
		);
}

Vector& Vector::crossEquals(const Vector& vec)
{
	return set(
		p[1]*vec[2] - p[2]*vec[1],
		p[2]*vec[0] - p[0]*vec[2],
		p[0]*vec[1] - p[1]*vec[0],		
		0.0f		
		);
}

float Vector::dot(const Vector& vec) const
{
	return p[0]*vec[0] + p[1]*vec[1] + p[2]*vec[2] + p[3]*vec[3]; 
}


Vector Vector::operator+(const Vector vec) const
{
	return Vector(
		p[0]+vec[0],
		p[1]+vec[1],
		p[2]+vec[2],
		p[3]+vec[3]);
}

Vector& Vector::operator+=(const Vector vec)
{
	return set(
		p[0]+vec[0],
		p[1]+vec[1],
		p[2]+vec[2],
		p[3]+vec[3]);
}

Vector Vector::operator-(const Vector vec) const
{
	return Vector(
		p[0]-vec[0],
		p[1]-vec[1],
		p[2]-vec[2],
		p[3]-vec[3]);
}

Vector& Vector::operator-=(const Vector vec)
{
	return set(
		p[0]-vec[0],
		p[1]-vec[1],
		p[2]-vec[2],
		p[3]-vec[3]);
}

Vector Vector::operator*(float scalar) const
{
	return Vector(p[0]*scalar, p[1]*scalar, p[2]*scalar, p[3]);
}

Vector& Vector::operator*=(float scalar)
{
	return set(p[0]*scalar, p[1]*scalar, p[2]*scalar, p[3]);
}

Vector Vector::operator-() const
{
	return Vector(-p[0], -p[1], -p[2], p[3]);
}

Vector& Vector::normalize()
{
	if ((float)fabs(p[3])<=EPS) {
		float length = (float)sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
		return set(p[0]/length,p[1]/length,p[2]/length,0.0f);
	}
	else {
		return set(p[0]/p[3], p[1]/p[3], p[2]/p[3], 1.0f);
	}
}

float Vector::norm()
{
	return (float)sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}

bool Vector::isBad()
{
	return (p[0]==0.0f && p[1]==0.0f && p[2]==0.0f && p[3]==0.0f);
}

Vector Vector::badVector()
{
	return Vector(0.0f, 0.0f, 0.0f, 0.0f);
}

Vector* Vector::clone() const
{
	return new Vector(*this);
}

Vector Vector::interpolate(const Vector& a, const Vector& b, float alpha)
{
	return Vector(
		a[0]*(1.0f-alpha)+b[0]*alpha,
		a[1]*(1.0f-alpha)+b[1]*alpha,
		a[2]*(1.0f-alpha)+b[2]*alpha,
		a[3]*(1.0f-alpha)+b[3]*alpha
		);
}

Vector Vector::cubicInterpolate(const Vector& a, const Vector& b, const Vector& c, const Vector& d, float alpha)
{
	// for now, just do the linear interpolation
	return interpolate(b,c,alpha);
}

}
