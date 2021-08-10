/**************************************************************************
 * FILE: mexAH.cpp														  *
 * DESCRIPTION: Make advection matrix AV   	  							  *
 *                                                             			  *
 * AUTHOR: Geir-Arne Fuglstad <geirarne.fuglstad@gmail.com>               *
 **************************************************************************/

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "mex.h"

// Decleration of structs
// Decleration of structs
struct Point2d{
	double x,y;
	Point2d(double x, double y):x(x),y(y){}
	Point2d(){}
	Point2d operator+(const Point2d& rhs){
		return Point2d(x+rhs.x, y+rhs.y);
	}
	Point2d operator-(const Point2d& rhs){
		return Point2d(x-rhs.x, y-rhs.y);
	}
};

struct Point3d{
	double x,y,z;
	Point3d(double x, double y, double z):x(x),y(y),z(z){}
	Point3d(){}
	Point3d operator+(const Point3d& rhs) const{
		return Point3d(x+rhs.x, y+rhs.y, z+rhs.z);
	}
	Point3d operator-(const Point3d& rhs) const{
		return Point3d(x-rhs.x, y-rhs.y, z-rhs.z);
	}
	double norm() const{
		return std::sqrt(std::pow(x, 2.0)+std::pow(y, 2.0)+std::pow(z, 2.0));
	}
	Point3d cross(const Point3d& rhs) const{
		return Point3d(y*rhs.z-z*rhs.y, z*rhs.x-x*rhs.z, x*rhs.y-y*rhs.x);
	}
	Point3d operator*(double val) const{
		return Point3d(x*val, y*val, z*val);
	}
	double dot(const Point3d& rhs) const{
		return (x*rhs.x+y*rhs.y+z*rhs.z);
	}
};

struct Point2i{
	int x,y;
	Point2i(int x,int y):x(x), y(y){}
	Point2i(){}
	bool operator<(const Point2i& rhs) const{
		return ((x < rhs.x) || ((x == rhs.x) && (y < rhs.y)));
	}
	bool operator==(const Point2i& rhs) const{
		return (x == rhs.x && y == rhs.y);
	}
	double operator[](int i) const{
		if(i == 0)
			return x;
		else
			return y;
	}
};

struct Point3{
	int x,y,z;
	Point3(int x, int y, int z):x(x), y(y), z(z){}
	Point3(){}
	double operator[](int i) const{
		if(i == 0)
			return x;
		else if(i == 1)
			return y;
		else
			return z;
	}
};

// Forward decleration of functions
void makeAV(std::vector<double>& iVec, std::vector<double>& jVec, std::vector<double>& vVec, const std::vector<Point3d>& loc, const std::vector<Point3d>& cLoc, const std::vector<Point3>& tt, const std::vector<Point3>& tv, double* vxa, double* vya, double* vza);
std::vector<Point3d> calcCentroids(const std::vector<Point3>& tv, const std::vector<Point3d>& loc);

/********************************************************************************
 * function: mexFuction															*
 * Description: Entry point for call from MATLAB								*
 * INPUT:																		*
 *       int      nlhs:   Number of output arguments defined in the pointers.	*
 *       mxArray* plhs[]: Vector to pointers to output arguments.				*
 *       int      nrhs:   Number of input arguments defined in the pointers.	*
 * const mxArray* prhs[]: Vector to const pointers to input argumens.			*
 ********************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*
         * Limit input checking to just verify that the correct number
         * of input and output parameters are specified.
         */
	// Check for proper number of input arguments
	if(nrhs != 8)
	{
		mexErrMsgTxt("EIGHT (8) input arguments required.");
	}
	// Check for proper number of output arguments
	if(nlhs != 3)
	{
		mexErrMsgTxt("THREE (3) output arguments required.");
	}

	// Get x and y sizes of grid
	double* nT = mxGetPr(prhs[0]);
	double* nV = mxGetPr(prhs[1]);

	// Get physical x and y size of cell
	double* locP = mxGetPr(prhs[2]);
	double* tvP = mxGetPr(prhs[3]);
	double* ttP = mxGetPr(prhs[4]);

	// Get pointers to v vector
	double* vx = mxGetPr(prhs[5]);
	double* vy = mxGetPr(prhs[6]);
	double* vz = mxGetPr(prhs[7]);

	// Waste some time to convert to format in c++ code
	int numT = (int)(nT[0]+0.5);
	int numV = (int)(nV[0]+0.5);

	std::vector<Point3d> vLoc;
	for(int i = 0; i < 3*numV; i+=3){
		vLoc.push_back(Point3d(locP[i], locP[i+1], locP[i+2]));
    }

    std::vector<Point3> tt;
    for(int i = 0; i < 3*numT; i+=3){
		tt.push_back(Point3((int)(ttP[i]-0.5), (int)(ttP[i+1]-0.5), (int)(ttP[i+2]-0.5)));
    }

    std::vector<Point3> tv;
    for(int i = 0; i < 3*numT; i+=3){
		tv.push_back(Point3((int)(tvP[i]-0.5), (int)(tvP[i+1]-0.5), (int)(tvP[i+2]-0.5)));
    }

    // Calc centroids
    std::vector<Point3d> cLoc = calcCentroids(tv, vLoc);

	// Make AV matrix
	std::vector<double> iVec, jVec, vVec;
	makeAV(iVec, jVec, vVec, vLoc, cLoc, tt, tv, vx, vy, vz);

	// Initialize storage for result
	plhs[0] = mxCreateDoubleMatrix(iVec.size(), 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(jVec.size(), 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(vVec.size(), 1, mxREAL);

	// Populate output vectors and convert to 1-indexing
	double* iPointer = mxGetPr(plhs[0]);
	double* jPointer = mxGetPr(plhs[1]);
	double* vPointer = mxGetPr(plhs[2]);
	for(unsigned int i = 0; i < iVec.size(); ++i){
		iPointer[i] = iVec[i]+1;
		jPointer[i] = jVec[i]+1;
		vPointer[i] = vVec[i];
    }
}


/*****************************************************
 * Make AV matrix									 *
 *    Represent sparse matrix by 3 vectors           *
 *****************************************************/
void makeAV(std::vector<double>& iVec, std::vector<double>& jVec, std::vector<double>& vVec, const std::vector<Point3d>& loc, const std::vector<Point3d>& cLoc, const std::vector<Point3>& tt, const std::vector<Point3>& tv, double* vxa, double* vya, double* vza){
	// Make sure nothing bad happens...
	iVec.clear();
	jVec.clear();
	vVec.clear();

	// Iterate through triangles
	for(unsigned int tr = 0; tr < tt.size(); ++tr){
		// Iterate through edges in triangle
		for(unsigned int edg = 0; edg < 3; ++edg){
			// Ignore boundary edges
			if(tt[tr][edg] == -1){
				continue;
			}

			// Form a quadrilateral on edge
			Point3d xS = loc[tv[tr][(edg+1)%3]];
			Point3d xE = cLoc[tt[tr][edg]];
			Point3d xN = loc[tv[tr][(edg+2)%3]];
			Point3d xW = cLoc[tr];

			// Find normal vectors of each part
			Point3d n1 = (xS-xW).cross(xN-xS);
			Point3d n2 = (xS-xN).cross(xE-xS);

			// Find bisecting vector to use
			n1 = n1*(1.0/n1.norm());
			n2 = n2*(1.0/n2.norm());
			Point3d nS = (n1+n2)*0.5;

			// Local axes
			Point3d e2 = xN-xS;
			Point3d e1 = e2.cross(nS);
			e1 = e1*(1.0/e1.norm());
			e2 = e2*(1.0/e2.norm());

			// Project vector field onto this local system
			Point3d vecField = Point3d(vxa[3*tr+edg], vya[3*tr+edg], vza[3*tr+edg]);
			Point2d vecF = Point2d(e1.dot(vecField), e2.dot(vecField));

			// Calculate normal vector to side
			Point2d n = Point2d(e2.dot(xN-xS), 0.0);

			// Calculate weight
			double w = vecF.x*n.x+vecF.y*n.y;

			// Add to matrix
			iVec.push_back(tr);
			vVec.push_back(w);
			if(w >= 0){
				jVec.push_back(tr);
			}
			else{
				jVec.push_back(tt[tr][edg]);
			}
		}
	}
}

/*****************************************************
 * Make vector of centroids							 *
 *****************************************************/
std::vector<Point3d> calcCentroids(const std::vector<Point3>& tv, const std::vector<Point3d>& loc){
	// Initialize vector to hold centroids
	std::vector<Point3d> cLoc;

	// Do it
	for(unsigned int i = 0; i < tv.size(); ++i){
		cLoc.push_back((loc[tv[i].x] + loc[tv[i].y] + loc[tv[i].z])*(1.0/3.0));
	}

	return cLoc;
}
