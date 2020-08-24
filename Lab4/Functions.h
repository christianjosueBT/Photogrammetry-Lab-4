// **********************************************************
// **********************************************************
// **********************************************************
// functions.h
// **********************************************************
// **********************************************************
// **********************************************************

#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <math.h>
#include <string>
#include <iomanip>

#define EIGEN_NO_DEBUG
#define PI 3.14159265358979323846264
#include <Eigen\Dense>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Class containing the exterior orientation parameters
class parameters
{
public:
	double x, y, z, omega, phi, kappa;
};
// Class that contains the 3D information of a point
class point
{
public:
	string id;
	double x, y, z;
};

// Functions
void read(const string& filename, vector<point>& coords);
void EOP(const vector<point>& coords, const vector<point>& coorected, const int& num, parameters& xhat);
void designMatrixA(const vector<point>& coords, const parameters& x0, MatrixXd& R, const double& C, MatrixXd& A);
void measurement(const vector<point>& coords, VectorXd& l);
void misclosure(const parameters& x0, const VectorXd& l, const vector<point>& corrected, const double& C, VectorXd& w);
void build_x0(const vector<point>& coords, const vector<point>& ground, const double& C, parameters& xhat);
void rotation(const parameters& xhat, MatrixXd& R);
void similarity(const vector<point>& image, const vector<point>& ground, VectorXd& xhat);
void correlation(const MatrixXd& Cxhat, MatrixXd& rho);
void outVec(const string& filename, const VectorXd& vec);
void outMat(const string& filename, const MatrixXd& mat);

// EOP partials
void Nx(const point& coords, const parameters& xhat, const MatrixXd& R, double& out);
void Ny(const point& coords, const parameters& xhat, const MatrixXd& R, double& out);
void D(const point& coords, const parameters& xhat, const MatrixXd& R, double& out);
void dxdx(const point& coords, const parameters& xhat, const MatrixXd& R, const double& U, const double& W, const double& C, double& out);
void dxdy(const point& coords, const parameters& xhat, const MatrixXd& R, const double& U, const double& W, const double& C, double& out);
void dxdz(const point& coords, const parameters& xhat, const MatrixXd& R, const double& U, const double& W, const double& C, double& out);
void dxdomega(const point& coords, const parameters& xhat, const MatrixXd& R, const double& U, const double& W, const double& C, double& out);
void dxdphi(const point& coords, const parameters& xhat, const double& U, const double& V, const double& W, const double& C, double& out);
void dxdkappa(const double& V, const double& W, const double& C, double& out);
void dydx(const point& coords, const parameters& xhat, const MatrixXd& R, const double& V, const double& W, const double& C, double& out);
void dydy(const point& coords, const parameters& xhat, const MatrixXd& R, const double& V, const double& W, const double& C, double& out);
void dydz(const point& coords, const parameters& xhat, const MatrixXd& R, const double& V, const double& W, const double& C, double& out);
void dydomega(const point& coords, const parameters& xhat, const MatrixXd& R, const double& V, const double& W, const double& C, double& out);
void dydphi(const point& coords, const parameters& xhat, const double& U, const double& V, const double& W, const double& C, double& out);
void dydkappa(const double& U, const double& W, const double& C, double& out);