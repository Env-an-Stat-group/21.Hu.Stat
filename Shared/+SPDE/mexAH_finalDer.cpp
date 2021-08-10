/**************************************************************************
 * FILE: mexAH.cpp														  *
 * DESCRIPTION: Make precision matrix according to FVM discretization     *
 *              of SPDE for sphere.                                       *
 * AUTHOR: Geir-Arne Fuglstad <geirarne.fuglstad@gmail.com>               *
 **************************************************************************/

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>

#include "mex.h"

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
void makeAH(std::vector<double>& iVec, std::vector<double>& jVec, std::vector<double>& vVec, const std::vector<Point3d>& cLoc, const std::vector<Point3d>& loc, const std::vector<Point3>& tt, const std::vector<Point3>& tv, const std::vector<std::vector<int> >& vt, double* vxa, double* vya, double* vza, double needDiff, double* vxDa, double* vyDa, double* vzDa, const std::vector<int> &sIdx, const std::vector<int> &tIdx, double* trW, double reScaleH);
double triArea(const Point3d& pA,const Point3d& pB, const Point3d& pC);
std::vector<Point3d> calcCentroids(const std::vector<Point3>& tv, const std::vector<Point3d>& loc);
std::vector<std::vector<int> > createVT(const std::vector<Point3>& tv);

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
	if(nrhs != 17)
	{
		mexErrMsgTxt("Seventeen (17) input arguments required.");
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

	// Get pointers to H matrices
	double* vx = mxGetPr(prhs[5]);
	double* vy = mxGetPr(prhs[6]);
	double* vz = mxGetPr(prhs[7]);

	// Get number saying if we need to differentiate
	double needDiff = *mxGetPr(prhs[8]);

	// Get vector field to differentiate with respect to
	double* vxD = mxGetPr(prhs[9]);
	double* vyD = mxGetPr(prhs[10]);
	double* vzD = mxGetPr(prhs[11]);
    
    // Get vertex reconstruction parts
    double* stIdx = mxGetPr(prhs[12]);
    double* trIdx = mxGetPr(prhs[13]);
    double* trW   = mxGetPr(prhs[14]);
    double* tmpLen= mxGetPr(prhs[15]);
    
    // Choose parametrization
    double reScaleH = *mxGetPr(prhs[16]);

	// Waste some time to convert to format in c++ code
	int numT = (int)(nT[0]+0.5);
	int numV = (int)(nV[0]+0.5);
    int numW = (int)(tmpLen[0]+0.5);
    
    // Make index vector
    std::vector<int> tIdx;
    for(int i = 0; i < numW; ++i){
        tIdx.push_back((int)(trIdx[i]+0.5));
    }
    std::vector<int> sIdx;
    for(int i = 0; i < numV+1; ++i){    
        sIdx.push_back((int)(stIdx[i]+0.5));
    }

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

    // Make vt map
	std::vector<std::vector<int> > vt = createVT(tv);

	// Calculate centroids
	std::vector<Point3d> cLoc = calcCentroids(tv, vLoc);

	// Make AH matrix
	std::vector<double> iVec, jVec, vVec;
	makeAH(iVec, jVec, vVec, cLoc, vLoc, tt, tv, vt, vx, vy, vz, needDiff, vxD, vyD, vzD, sIdx, tIdx, trW, reScaleH);

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
 * Make AH matrix									 *
 *    Represent sparse matrix by 3 vectors           *
 *****************************************************/
void makeAH(std::vector<double>& iVec, std::vector<double>& jVec, std::vector<double>& vVec, const std::vector<Point3d>& cLoc, const std::vector<Point3d>& loc, const std::vector<Point3>& tt, const std::vector<Point3>& tv, const std::vector<std::vector<int> >& vt, double* vxa, double* vya, double* vza, double needDiff, double* vxDa, double* vyDa, double* vzDa, const std::vector<int> &sIdx, const std::vector<int> &tIdx, double* trW, double reScaleH){
	// Make sure nothing bad happens...
	iVec.clear();
	jVec.clear();
	vVec.clear();

	// Iterate through triangles
	for(unsigned int tr = 0; tr < tt.size(); ++tr){
		// Iterate through edges in triangle
		for(unsigned int edg = 0; edg < 3; ++edg){
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

			// Local coordinates for easy point
			Point2d yS = Point2d(0.0, 0.0);
			Point2d yN = Point2d(0.0, e2.dot(xN-xS));
            
            // This is a bit tricky. Basically shift and fold the triangles
            // up so that they lie in the tangent space where the vector
            // is evaluated.

			// Project points and correct length according to angle
			double corr = (n1.norm()*nS.norm())/n1.dot(nS);
			Point2d yW = Point2d(e1.dot(xW-xS)*corr, e2.dot(xW-xS));
			Point2d yE = Point2d(e1.dot(xE-xS)*corr, e2.dot(xE-xS));

			// Project anisotropy vector field into plane as well;
            // With correction for length
            Point3d edgeC = (xN+xS)*0.5;
            double corr2 = (edgeC.norm()*nS.norm())/edgeC.dot(nS);
			Point3d vecField = Point3d(vxa[3*tr+edg], vya[3*tr+edg], vza[3*tr+edg]);
			Point2d vecF = Point2d(e1.dot(vecField)*corr2, e2.dot(vecField));

			// Project derivative vector field into plane as well
			Point3d vecDField = Point3d(vxDa[3*tr+edg], vyDa[3*tr+edg], vzDa[3*tr+edg]);
			Point2d vecDF = Point2d(e1.dot(vecDField)*corr2, e2.dot(vecDField));

			// Make a H matrix corresponding to this
			double fac = std::sqrt(1+std::pow(vecF.x, 2.0)+std::pow(vecF.y, 2.0));
            if(reScaleH < 0.5){
                fac = 1.0;
            }
			double H11, H12, H22;
			if(needDiff < 0.5){			
				H11 = (1.0+std::pow(vecF.x, 2.0))/fac;
				H12 = vecF.x*vecF.y/fac;
				H22 = (1.0+std::pow(vecF.y, 2.0))/fac;
			}
			else{
                if(reScaleH > 0.5){
                    H11 = -(1+std::pow(vecF.x, 2.0))*(vecF.x*vecDF.x+vecF.y*vecDF.y)/std::pow(fac, 3.0) + 2.0*vecF.x*vecDF.x/fac;
                    H12 = -(vecF.x*vecF.y)*(vecF.x*vecDF.x+vecF.y*vecDF.y)/std::pow(fac, 3.0) + (vecDF.x*vecF.y+vecF.x*vecDF.y)/fac;
                    H22 = -(1+std::pow(vecF.y, 2.0))*(vecF.x*vecDF.x+vecF.y*vecDF.y)/std::pow(fac, 3.0) + 2.0*vecF.y*vecDF.y/fac;
                } else{
                    H11 = 2.0*vecF.x*vecDF.x;
                    H12 = vecDF.x*vecF.y + vecF.x*vecDF.y;
                    H22 = 2.0*vecF.y*vecDF.y;
                }
			}
            
            // Help variables
            double lamA = 1.0-yW.y;
            double lamB = yE.y;
            double E = yN.y-yS.y;
            double H = yE.x-yW.x;
            
            // Help flux
            double kN = H11;
            double kT = (lamA+lamB-1)*E/H*H11 + H12;
            
            // Coefficients
            double wW = -kN/H*E;
            double wE = kN/H*E;
            double wN = kT/E*E;
            double wS = -kT/E*E;
            
			// Coefficient w.r.t. this curr triangle
			iVec.push_back(tr);
			jVec.push_back(tr);
			vVec.push_back(wW);

			// Coefficient w.r.t. neighbour
			iVec.push_back(tr);
			jVec.push_back(tt[tr][edg]);
			vVec.push_back(wE);
            
            // Vertex reconstruction using provided interpolation algorithm
            int northV = tv[tr][(edg+2)%3];
            for(int i = sIdx[northV]; i < sIdx[northV+1]; ++i){
                // North vertex
                iVec.push_back(tr);
                jVec.push_back(tIdx[i]);
                vVec.push_back(wN*trW[i]);
            }
            
            int southV = tv[tr][(edg+1)%3];
            for(int i = sIdx[southV]; i < sIdx[southV+1]; ++i){
                // North vertex
                iVec.push_back(tr);
                jVec.push_back(tIdx[i]);
                vVec.push_back(wS*trW[i]);
            }
		}
	}
}

/*****************************************************
 * Area of triangle                                  *
 *    Copy of Marco Zuliani - zuliani@ece.ucsb.edu   *
 *    file for MATLAB                                *
 *****************************************************/
double triArea(const Point3d& pA,const Point3d& pB, const Point3d& pC){
	// Calculate length of sides
	double s[3];
	s[0] = (pA-pB).norm();
	s[1] = (pC-pA).norm();
	s[2] = (pC-pB).norm();

	// Sort
	std::sort(s, s+3);

	// Do stabilized calculations
	double tmp = s[1] + s[0];
	double v1 = s[2] + tmp;
	tmp = s[2] - s[1];
	double v2 = s[0] - tmp;
	double v3 = s[0] + tmp;
	tmp = s[1] - s[0];
	double v4 = s[2] + tmp;
	return 0.25*std::sqrt(v1*v2*v3*v4);
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

/******************************************************
 * Create vertex to triangle map                      *
 ******************************************************/
std::vector<std::vector<int> > createVT(const std::vector<Point3>& tv){
	// Make array containing all reverse mappings
	std::vector<Point2i> tmpVT;
	for(unsigned int i = 0; i < tv.size(); ++i){
		tmpVT.push_back(Point2i(tv[i].x, i));
		tmpVT.push_back(Point2i(tv[i].y, i));
		tmpVT.push_back(Point2i(tv[i].z, i));
	}

	// Sort according to first coordinate and then second cordinate
	std::sort(tmpVT.begin(), tmpVT.end());

	// Remove duplicates
	std::unique(tmpVT.begin(), tmpVT.end());

	// Make neighbour list
	std::vector<std::vector<int> > tmpList;
	tmpList.resize(tmpVT[tmpVT.size()-1].x+1);
	for(unsigned int i = 0; i < tmpVT.size(); ++i){
		tmpList[tmpVT[i].x].push_back(tmpVT[i].y);
	}

	return tmpList;
}
