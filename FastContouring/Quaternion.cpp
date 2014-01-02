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

// Quaternion.cpp: implementation of the Quaternion class.
//
//////////////////////////////////////////////////////////////////////

#include "Quaternion.h"
#include "Vector.h"
#include "Ray.h"
#include "Matrix.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

namespace FastContouring
{

Quaternion::Quaternion(float w, float x, float y, float z) : Tuple(w,x,y,z)
{
}

Quaternion::Quaternion() : Tuple(1.0, 0.0, 0.0, 0.0)
{

}

Quaternion::Quaternion(const Vector& vec)
{
	p[0] = 0.0f;
	p[1] = vec[0];
	p[2] = vec[1];
	p[3] = vec[2];
}

Quaternion::~Quaternion()
{

}

Quaternion::Quaternion(const Quaternion& copy): Tuple(copy)
{

}

Quaternion& Quaternion::operator=(const Quaternion& copy)
{
	if (this!=&copy) {
		set(copy);
	}
	return *this;
}

Quaternion& Quaternion::set(float w, float x, float y, float z)
{
	Tuple::set(w,x,y,z);
	return *this;
}

Quaternion& Quaternion::set(float* array)
{
	Tuple::set(array);
	return *this;
}

Quaternion& Quaternion::set(const Quaternion& copy)
{
	Tuple::set(copy);
	return *this;
}

Quaternion Quaternion::operator+(const Quaternion& quat) const
{
	return Quaternion(p[0]+quat.p[0],
		p[1]+quat.p[1],
		p[2]+quat.p[2],
		p[3]+quat.p[3]);
}

Quaternion Quaternion::operator-(const Quaternion& quat) const
{
	return Quaternion(p[0]-quat.p[0],
		p[1]-quat.p[1],
		p[2]-quat.p[2],
		p[3]-quat.p[3]);
}

Quaternion Quaternion::operator-() const
{
	return Quaternion(-p[0],
		-p[1],
		-p[2],
		-p[3]);
}

Quaternion Quaternion::operator*(const Quaternion& quat) const
{
	return Quaternion(
		p[0]*quat[0] - p[1]*quat[1] - p[2]*quat[2] - p[3]*quat[3],
		p[0]*quat[1] + p[1]*quat[0] + p[2]*quat[3] - p[3]*quat[2],
		p[0]*quat[2] - p[1]*quat[3] + p[2]*quat[0] + p[3]*quat[1],
		p[0]*quat[3] + p[1]*quat[2] - p[2]*quat[1] + p[3]*quat[0]
		);
}

Quaternion& Quaternion::preMultiply(const Quaternion& quat)
{
	return set(quat*(*this));
}

Quaternion& Quaternion::postMultiply(const Quaternion& quat)
{
	return set((*this)*quat);
}

Quaternion& Quaternion::rotate(float angle, float x, float y, float z)
{
	return preMultiply(rotation(angle, x, y, z));
}

Quaternion Quaternion::operator*(float scalar) const
{
	return Quaternion(p[0]*scalar, p[1]*scalar, p[2]*scalar, p[3]*scalar);
}

Quaternion& Quaternion::operator*=(float scalar)
{
	return set(p[0]*scalar, p[1]*scalar, p[2]*scalar, p[3]*scalar);
}

Quaternion Quaternion::operator/(float scalar) const
{
	return Quaternion(p[0]/scalar, p[1]/scalar, p[2]/scalar, p[3]/scalar);
}

Quaternion& Quaternion::operator/=(float scalar)
{
	return set(p[0]/scalar, p[1]/scalar, p[2]/scalar, p[3]/scalar);
}

Quaternion& Quaternion::normalize()
{
	float n = norm();
	return (*this)/=n;
}

Quaternion Quaternion::conjugate() const
{
	return Quaternion(p[0], -p[1], -p[2], -p[3]);
}

Quaternion Quaternion::inverse() const
{
	return conjugate()/norm();
}

float Quaternion::norm() const
{
	return (float)sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
}

Vector Quaternion::applyRotation(const Vector& vec) const
{
	Quaternion result = (*this) * Quaternion(vec) * (conjugate());
	return Vector(result[1], result[2], result[3], vec[3]);
}

Ray Quaternion::applyRotation(const Ray& ray) const
{
	Quaternion origin = (*this) * Quaternion(ray.m_Origin) * (conjugate());
	Quaternion dir = (*this) * Quaternion(ray.m_Dir) * (conjugate());
	return Ray(Vector(origin[1], origin[2], origin[3], ray.m_Origin[3]),
		Vector(dir[1], dir[2], dir[3], ray.m_Dir[3]));
}

Matrix Quaternion::buildMatrix() const
{
	float w = p[0];
	float x = p[1];
	float y = p[2];
	float z = p[3];
	return Matrix(
		1.0f-2.0f*y*y-2.0f*z*z,	2.0f*x*y-2.0f*w*z,		2.0f*x*z + 2.0f*w*y, 0.0f,
		2.0f*x*y + 2.0f*w*z,	1.0f - 2.0f*x*x - 2.0f*z*z,	2.0f*y*z - 2.0f*w*x, 0.0f,
        2.0f*x*z - 2.0f*w*y,	2.0f*y*z + 2.0f*w*x,		1.0f - 2.0f*x*x - 2.0f*y*y, 0.0f,
		0.0f,0.0f,0.0f,1.0f
		);
}

Quaternion Quaternion::power(double scalar) const
{
	float Dest[4];
	// determine theta
	double theta;
	if (p[0]>=0.9999f) {
		theta = 0;
	}
	else if (p[0]<=-0.9999f) {
		theta = 2.0*3.1415926535897932384626433832795;
	}
	else {
		theta = acos(p[0]);
	}

	// determine u
	double u[3];
	double scale = p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
	scale = sqrt(scale);
	if (p[1]==0.0f && p[2]==0.0f && p[3]==0.0f) {
		u[0] = 0.0;
		u[1] = 0.0;
		u[2] = 0.0;
	}
	else {
		u[0] = p[1]/scale;
		u[1] = p[2]/scale;
		u[2] = p[3]/scale;
	}

	Dest[0] = (float)cos(scalar*theta);
	Dest[1] = (float)(u[0] * sin(scalar*theta));
	Dest[2] = (float)(u[1] * sin(scalar*theta));
	Dest[3] = (float)(u[2] * sin(scalar*theta));
	return Quaternion( Dest[0], Dest[1], Dest[2], Dest[3] );
}

Quaternion Quaternion::log(const Quaternion& a)
{
	Quaternion b(a);
	b.normalize();

	float Dest[4];
	// determine theta
	double theta;
	if (b.p[0]>=0.9999f) {
		theta = 0;
	}
	else if (b.p[0]<=-0.9999f) {
		theta = 2.0*3.1415926535897932384626433832795;
	}
	else {
		theta = acos(b.p[0]);
	}

	// determine u
	double u[3];
	double scale = b.p[1]*b.p[1]+b.p[2]*b.p[2]+b.p[3]*b.p[3];
	scale = sqrt(scale);
	if (b.p[1]==0.0f && b.p[2]==0.0f && b.p[3]==0.0f) {
		u[0] = 0.0;
		u[1] = 0.0;
		u[2] = 0.0;
	}
	else {
		u[0] = b.p[1]/scale;
		u[1] = b.p[2]/scale;
		u[2] = b.p[3]/scale;
	}

	Dest[0] = (float)0.0;
	Dest[1] = (float)(u[0] * theta);
	Dest[2] = (float)(u[1] * theta);
	Dest[3] = (float)(u[2] * theta);
	return Quaternion( Dest[0], Dest[1], Dest[2], Dest[3] );
}

Quaternion Quaternion::exponent(const Quaternion& a)
{
	float Dest[4];
	double scale, u[3];
	scale = sqrt(a.p[1]*a.p[1]+a.p[2]*a.p[2]+a.p[3]*a.p[3]);
	if (a.p[1]==0.0f && a.p[2]==0.0f && a.p[3]==0.0f) {
		u[0] = 0.0;
		u[1] = 0.0;
		u[2] = 0.0;
	}
	else {
		u[0] = a.p[1]/scale;
		u[1] = a.p[2]/scale;
		u[2] = a.p[3]/scale;
	}
	
	Dest[0] = (float)((double)exp((double)a.p[0]) * cos(scale));
	Dest[1] = (float)((double)exp((double)a.p[0]) * u[0] * sin(scale));
	Dest[2] = (float)((double)exp((double)a.p[0]) * u[1] * sin(scale));
	Dest[3] = (float)((double)exp((double)a.p[0]) * u[2] * sin(scale));
	return Quaternion( Dest[0], Dest[1], Dest[2], Dest[3] );
}

Quaternion Quaternion::rotation(float angle, float x, float y, float z)
{
	float len = (float)sqrt(x*x+y*y+z*z);
	if (len!=0.0) {
		len = (float)(sin(angle/2.0f)/len);
		return Quaternion((float) cos(angle/2.0f), x*len, y*len, z*len);
	}
	else {
		return Quaternion();
	}
}

Quaternion Quaternion::rotation(float angle, const Vector& axis)
{
	float len = (float)sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
	if (len!=0.0) {
		len = (float)(sin(angle/2.0f)/len);
		return Quaternion((float) cos(angle/2.0f), axis[0]*len, axis[1]*len, axis[2]*len);
	}
	else {
		return Quaternion();
	}
}

Quaternion Quaternion::interpolate(const Quaternion& a, const Quaternion& b, float alpha)
{
	Quaternion q = (a.inverse()*b);
	return a*q.power(alpha);
}

Quaternion Quaternion::startCubicInterpolate(const Quaternion& a, const Quaternion& b, const Quaternion& c, float alpha)
{
	return cubicInterpolate(a, a, b, c, alpha);
}

Quaternion Quaternion::endCubicInterpolate(const Quaternion& a, const Quaternion& b, const Quaternion& c, float alpha)
{
	return cubicInterpolate(a, b, c, c, alpha);
}

Quaternion Quaternion::cubicInterpolate(const Quaternion& a, const Quaternion& b, const Quaternion& c, const Quaternion& d, float alpha)
{
	// compute q1 and q2
	Quaternion q1 = determineQi(a, b, c);
	Quaternion q2 = determineQi(b, c, d);

	// perform interpolation
	return interpolate( interpolate(b,c,alpha),
		interpolate(q1, q2, alpha),
		2.0f*alpha*(1.0f-alpha));
}

Quaternion Quaternion::determineQi(const Quaternion& aiminus1, const Quaternion& ai, const Quaternion& aiplus1)
{
	Quaternion temp = -(log(ai.inverse()*aiplus1) + log(ai.inverse()*aiminus1))/4.0;
	Quaternion qi = ai * exponent(temp);
	return qi;
}

}

