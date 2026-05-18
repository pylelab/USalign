/* File parsing and basic geometry operations */
#ifndef TMalign_GeometryTools_h
#define TMalign_GeometryTools_h 1

#include "basic_fun.h" // For reading gzip and bz2 compressed files

using namespace std;

const double Extra=1.0e-4; // pseudocount
const float PI=3.141592653589793;

// subtract, norm, vectorsum, CoordinateRotation
void subtract(double *c1, double *c2, double *cc)
{
    cc[0] = c1[0] - c2[0];
    cc[1] = c1[1] - c2[1];
    cc[2] = c1[2] - c2[2];
}

bool norm(double *c)
{
    double magnitude = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
    if (magnitude<Extra) return false;
    c[0] /=magnitude;
    c[1] /=magnitude;
    c[2] /=magnitude;
    return true;
}

bool norm(double *c,double *cc)
{
    double magnitude = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
    if (magnitude<Extra) return false;
    cc[0]=c[0] /magnitude;
    cc[1]=c[1] /magnitude;
    cc[2]=c[2] /magnitude;
    return true;
}

void vectorsum(double *c1, double *c2, double *cc)
{
    cc[0] = c1[0] + c2[0];
    cc[1] = c1[1] + c2[1];
    cc[2] = c1[2] + c2[2];
    return;
}


inline double deg2rad(double deg)
{
   return deg*PI/180;
}

/***************************CoordinateRotation********************
* For a PointA, Given the Coordinates of 2 points(axisA, axisB) 
* in line of Rotation Axis and the Rotation Angle (degree)
* Generate the coordinate of PointB after operation
* Author: C.Y.
* Date 2006.12.2.  
* Author: C.Z.
* Date 2026.5.13  
*****************************END**********************************/
bool CoordinateRotation(double *pointA, double *axisA, double *axisB,
    double angle, double *pointB)
{
    double axis[3];
    axis[0]=axisB[0]-axisA[0]; 
    axis[1]=axisB[1]-axisA[1];
    axis[2]=axisB[2]-axisA[2]; 
   
    double ouc[3];
    if(!norm(axis, ouc))
    {
        pointB[0]=pointA[0];
        pointB[1]=pointA[1];
        pointB[2]=pointA[2];
        return false;
    }
    double c=cos(deg2rad(angle));
    double s=sin(deg2rad(angle));
    double t=1-c;
    
    double rotmtx[3][3];
    rotmtx[0][0]=t*ouc[0]*ouc[0]+c;
    rotmtx[0][1]=t*ouc[0]*ouc[1]+s*ouc[2];
    rotmtx[0][2]=t*ouc[0]*ouc[2]-s*ouc[1];
   
    rotmtx[1][0]=t*ouc[1]*ouc[0]-s*ouc[2];
    rotmtx[1][1]=t*ouc[1]*ouc[1]+c;
    rotmtx[1][2]=t*ouc[1]*ouc[2]+s*ouc[0];
   
    rotmtx[2][0]=t*ouc[2]*ouc[0]+s*ouc[1];
    rotmtx[2][1]=t*ouc[2]*ouc[1]-s*ouc[0];
    rotmtx[2][2]=t*ouc[2]*ouc[2]+c;
   
    double point_A[3];
    point_A[0]=pointA[0]-axisA[0];
    point_A[1]=pointA[1]-axisA[1];
    point_A[2]=pointA[2]-axisA[2];
   
    pointB[0]=dot(&rotmtx[0][0],point_A)+axisA[0];
    pointB[1]=dot(&rotmtx[1][0],point_A)+axisA[1];
    pointB[2]=dot(&rotmtx[2][0],point_A)+axisA[2];
    return true;
}

#endif
