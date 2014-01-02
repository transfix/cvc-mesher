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

// Ray.cpp: implementation of the Ray class.
//
//////////////////////////////////////////////////////////////////////

#include "Ray.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

namespace FastContouring
{

Ray::Ray() : m_Origin(0.0f, 0.0f, 0.0f, 1.0f), m_Dir(0.0f, 0.0f, 1.0f, 0.0f)
{
}

Ray::Ray(const Vector& origin, const Vector& dir) : m_Origin(origin), m_Dir(dir)
{
	
}

Ray::~Ray()
{

}

Vector Ray::getPointOnRay(float t) const
{
	return m_Origin+m_Dir*t;
}

float Ray::nearestTOnXAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	//	float distance = ray.distanceToXAxis(Origin);
	float t = -(ray.m_Origin[1]*ray.m_Dir[1] + ray.m_Origin[2]*ray.m_Dir[2])/
		((ray.m_Dir[1]*ray.m_Dir[1]+ray.m_Dir[2]*ray.m_Dir[2]) );
	return t;
}

float Ray::nearestTOnYAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	//	float distance = ray.distanceToYAxis(Origin);
	float t = -(ray.m_Origin[0]*ray.m_Dir[0] + ray.m_Origin[2]*ray.m_Dir[2])/
		((ray.m_Dir[0]*ray.m_Dir[0]+ray.m_Dir[2]*ray.m_Dir[2]));
	return t;
}

float Ray::nearestTOnZAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	//	float distance = ray.distanceToZAxis(Origin);
	float t = -(ray.m_Origin[0]*ray.m_Dir[0] + ray.m_Origin[1]*ray.m_Dir[1])/
		((ray.m_Dir[1]*ray.m_Dir[1]+ray.m_Dir[0]*ray.m_Dir[0]));
	return t;
}

Vector Ray::nearestPointOnXAxis(Vector Origin) const
{
	Origin[3] = 0;
	float t = nearestTOnXAxis(Origin);
	Vector result = getPointOnRay(t);
	// project to axis
	result[1] = Origin[1];
	result[2] = Origin[2];
	//result+=Origin;
	return result;
}

Vector Ray::nearestPointOnYAxis(Vector Origin) const
{
	Origin[3] = 0;
	float t = nearestTOnYAxis(Origin);
	Vector result = getPointOnRay(t);
	// project to axis
	result[0] = Origin[0];
	result[2] = Origin[2];
	//result+=Origin;
	return result;
}

Vector Ray::nearestPointOnZAxis(Vector Origin) const
{
	Origin[3] = 0;
	float t = nearestTOnZAxis(Origin);
	Vector result = getPointOnRay(t);
	// project to axis
	result[0] = Origin[0];
	result[1] = Origin[1];
	return result;
}

float Ray::distanceToXAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	return (float)fabs( 
		( ray.m_Origin[2]*ray.m_Dir[1]-ray.m_Origin[1]*m_Dir[2] ) /
		(float)sqrt( ray.m_Dir[2]*ray.m_Dir[2] + ray.m_Dir[1]*ray.m_Dir[1] )
		);
}

float Ray::distanceToYAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	return (float)fabs( 
		( ray.m_Origin[2]*ray.m_Dir[0]-ray.m_Origin[0]*m_Dir[2] ) /
		(float)sqrt( ray.m_Dir[2]*ray.m_Dir[2] + ray.m_Dir[0]*ray.m_Dir[0] )
		);
}

float Ray::distanceToZAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	return (float)fabs( 
		( ray.m_Origin[0]*ray.m_Dir[1]-ray.m_Origin[1]*m_Dir[0] ) /
		(float)sqrt( ray.m_Dir[0]*ray.m_Dir[0] + ray.m_Dir[1]*ray.m_Dir[1] )
		);
}

}
