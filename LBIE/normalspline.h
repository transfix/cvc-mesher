#include <math.h>
  
#include <stdio.h>
#include <stdlib.h>
  
#include <string.h>

namespace LBIE
{

float   InitialAntiCausalCoefficient(float *, int, float);	
float	InitialCausalCoefficient(float *, int, float, float);
void	ConvertToInterpolationCoefficients(float *, int, float *, int ,float);

void    TransImg2Spline(float *, float *, int, int, int);

double  BS_Fun(double);

double  BS_GraFun(double);
					
void GradientAtPoint(float *,float , float , float , int, int , int, float *);

}
