// **********************************************************
// **********************************************************
// **********************************************************
// functions.cpp
// **********************************************************
// **********************************************************
// **********************************************************

#include "Functions.h"

// Function to read data from a file to a vector of class point
void read(const string& filename, vector<point>& coords) {

	point p;
	ifstream in;
	in.open(filename.c_str());
	if (in.fail()) {
		cout << "Error reading file!" << endl;
		return;
	}
	while (in >> p.id >> p.x >> p.y >> p.z) {
		coords.push_back(p);
	}
	return;
};
// Function to compute the EOP
void EOP(const vector<point>& coords, const vector<point>& corrected, const int& num, parameters& xhat) {

	MatrixXd P;
	P.setIdentity(coords.size() * 2, coords.size() * 2);

	parameters x0;
	double lambda;
	double focal = 153.167;
	double  U, W, V;
	build_x0(corrected, coords, focal, x0);

	VectorXd l, w;
	measurement(corrected, l);
	MatrixXd A, delta, N, u, R;
	delta.setOnes(6, 1);
	xhat = x0;
	
	rotation(x0, R);
	designMatrixA(coords, x0, R, focal, A);

	misclosure(x0, l, coords, focal, w);

	N = A.transpose() * P * A;
	u = A.transpose() * P * w;
	delta = -(N.inverse()) * u;

	double threshold = abs(delta(0));
	while (threshold >= 0.0000001) {

		rotation(x0, R);
		designMatrixA(coords, x0, R, focal, A);
		misclosure(x0, l, coords, focal, w);
		N = A.transpose() * P * A;
		u = A.transpose() * P * w;
		delta = -N.inverse() * u;

		xhat.x = x0.x + delta(0, 0);
		xhat.y = x0.y + delta(1, 0);
		xhat.z = x0.z + delta(2, 0);
		xhat.omega = x0.omega + delta(3, 0);
		xhat.phi = x0.phi + delta(4, 0);
		xhat.kappa = x0.kappa + delta(5, 0);

		x0 = xhat;

		threshold = abs(delta(0));
		for (int i = 0; i < delta.rows(); i++) {
			if (threshold < abs(delta(i)))
				threshold = abs(delta(i));
		}
	}

	// residuals
	VectorXd vhat = A * delta + w;
	string out_string;
	stringstream out;
	out << "Residuals_" << num << ".txt";
	out_string = out.str();
	outVec(out_string, vhat);

	// aposteriori variance factor
	VectorXd aposteriori = (vhat.transpose() * P * vhat) / (A.rows() - A.cols());
	double value = (vhat.transpose() * P * vhat)(0, 0);
	double apost = value / (A.rows() - A.cols());
	stringstream out1;
	out1 << "Aposteriori_" << num << ".txt";
	out_string = out1.str();
	outVec(out_string, aposteriori);

	// covariance matrix
	MatrixXd Cxhat = MatrixXd::Zero(delta.size(), delta.size());
	Cxhat = apost * N.inverse();
	stringstream out2;
	out2 << "CxHat_" << num << ".txt";
	out_string = out2.str();
	outMat(out_string, Cxhat);

	// correlation coefficient matrix
	MatrixXd rho;
	correlation(apost * Cxhat, rho);
	stringstream out3;
	out3 << "Correlation_" << num << ".txt";
	out_string = out3.str();
	outMat(out_string, rho);

	// standard deviation
	VectorXd stdSquared = Cxhat.diagonal();
	stringstream out4;
	out4 << "StDevSquared_" << num << ".txt";
	out_string = out4.str();
	outVec(out_string, stdSquared);

	return;
};
// Function to build the design matrix
void designMatrixA(const vector<point>& coords, const parameters& x0, MatrixXd& R, const double& C, MatrixXd& A) {

	A.resize(2 * coords.size(), 6);
	A.setZero(A.rows(), A.cols());
	rotation(x0, R);

	for (int i = 0; i < coords.size(); i++) {

		double nx;
		Nx(coords[i], x0, R, nx);
		double ny;
		Ny(coords[i], x0, R, ny);
		double d;
		D(coords[i], x0, R, d);
		double out;
		dxdx(coords[i], x0, R, nx, d, C, out);
		A(2 * i, 0) = out;
		dxdy(coords[i], x0, R, nx, d, C, out);
		A(2 * i, 1) = out;
		dxdz(coords[i], x0, R, nx, d, C, out);
		A(2 * i, 2) = out;
		dxdomega(coords[i], x0, R, nx, d, C, out);
		A(2 * i, 3) = out;
		dxdphi(coords[i], x0, nx, ny, d, C, out);
		A(2 * i, 4) = out;
		dxdkappa(ny, d, C, out);
		A(2 * i, 5) = out;
		dydx(coords[i], x0, R, ny, d, C, out);
		A(2 * i + 1, 0) = out;
		dydy(coords[i], x0, R, ny, d, C, out);
		A(2 * i + 1, 1) = out;
		dydz(coords[i], x0, R, ny, d, C, out);
		A(2 * i + 1, 2) = out;
		dydomega(coords[i], x0, R, ny, d, C, out);
		A(2 * i + 1, 3) = out;
		dydphi(coords[i], x0, nx, ny, d, C, out);
		A(2 * i + 1, 4) = out;
		dydkappa(nx, d, C, out);
		A(2 * i + 1, 5) = out;
	}
	return;
};
// Function to compute the measurement vector
void measurement(const vector<point>& coords, VectorXd& l) {

	l.resize(coords.size() * 2);
	for (int i = 0; i < coords.size(); i++)
	{
		l[2 * i] = coords[i].x;
		l[2 * i + 1] = coords[i].y;
	}
	return;
};
// Function to compute the misclosure vector
void misclosure(const parameters& x0, const VectorXd& l, const vector<point>& corrected, const double& C, VectorXd& w) {

	w.resize(l.size());
	MatrixXd R;
	rotation(x0, R);

	for (int i = 0; i < l.size() / 2; i++) {

		double Uval;
		Nx(corrected[i], x0, R, Uval);
		double Vval;
		Ny(corrected[i], x0, R, Vval);
		double Wval;
		D(corrected[i], x0, R, Wval);

		w(2 * i) = (-1 * C * Uval / Wval) - l[2 * i];
		w(2 * i + 1) = (-1 * C * (Vval / Wval)) - l[2 * i + 1];
	}
	return;
};
// Function to compute the initial approximate values
void build_x0(const vector<point>& coords, const vector<point>& ground, const double& C, parameters& xhat) {

	double H = 2000; 
	xhat.omega = 0;
	xhat.phi = 0;

	VectorXd similarity_par;
	similarity(coords, ground, similarity_par);

	double a = similarity_par(0);
	double b = similarity_par(1);

	double kappa = 2 * PI - atan2(b, a);
	xhat.kappa = kappa - (2 * PI) * floor(kappa / (2 * PI));
	double lambda = sqrt(pow(a, 2) + pow(b, 2));
	xhat.x = similarity_par[2];
	xhat.y = similarity_par[3];
	xhat.z = H;

	return;
}
// Function to compute the rotation matrix
void rotation(const parameters& xhat, MatrixXd& R)
{
	R.resize(3, 3);
	R << cos(xhat.phi) * cos(xhat.kappa), -cos(xhat.phi) * sin(xhat.kappa), sin(xhat.phi),
		cos(xhat.omega)* sin(xhat.kappa) + sin(xhat.omega) * sin(xhat.phi) * sin(xhat.kappa),
		cos(xhat.omega)* cos(xhat.kappa) - sin(xhat.omega) * sin(xhat.phi) * sin(xhat.kappa),
		-sin(xhat.omega) * cos(xhat.phi),
		sin(xhat.omega)* sin(xhat.kappa) - cos(xhat.omega) * sin(xhat.phi) * cos(xhat.kappa),
		sin(xhat.omega)* cos(xhat.kappa) + cos(xhat.omega) * sin(xhat.phi) * sin(xhat.kappa),
		cos(xhat.omega)* cos(xhat.phi);
	return;
};
// Function to compute the similarity transformation vector
void similarity(const vector<point>& image, const vector<point>& ground, VectorXd& xhat) {

	MatrixXd A(2 * image.size(), 4);
	A.setZero(A.rows(), A.cols());

	for (int i = 0; i < image.size(); i++)
	{
		A(2 * i, 0) = image[i].x;
		A(2 * i, 1) = image[i].y;
		A(2 * i, 2) = 1;
		A(2 * i, 3) = 0;

		A(2 * i + 1, 0) = image[i].y;
		A(2 * i + 1, 1) = -1 * image[i].x;
		A(2 * i + 1, 2) = 0;
		A(2 * i + 1, 3) = 1;
	}

	VectorXd l(A.rows());
	for (int i = 0; i < ground.size(); i++)
	{
		l[2 * i] = ground[i].x;
		l[2 * i + 1] = ground[i].y;
	}
	MatrixXd N = (A.transpose() * A).inverse();
	VectorXd u = A.transpose() * l;

	xhat = N * u;
	return;
};
// Function to compute the correlation matrix
void correlation(const MatrixXd& Cxhat, MatrixXd& rho) {

	rho.resize(Cxhat.rows(), Cxhat.cols());

	for (int i = 0; i < rho.rows(); i++) {
		for (int j = 0; j < rho.cols(); j++) {
			rho(i, j) = Cxhat(i, j) / (sqrt(Cxhat(i, i)) * sqrt(Cxhat(j, j)));
		}
	}
	return;
};
// Function to write the contents of a vector to an output file
void outVec(const string& filename, const VectorXd& vec) {

	ofstream out;
	out.open(filename.c_str());
	if (out.fail()) {
		cout << "Error writing file!" << endl;
		return;
	}
	for (int i = 0; i < vec.size(); i++) {
		out << fixed << setprecision(8) << vec(i) << "\n";
	}
	out.close();
	return;
};
// Function to write the contents of a matrix to an output file
void outMat(const string& filename, const MatrixXd& mat) {

	ofstream out;
	out.open(filename.c_str());
	if (out.fail()) {
		cout << "Error writing file!" << endl;
		return;
	}
	for (int i = 0; i < mat.rows(); i++) {
		for (int j = 0; j < mat.cols(); j++) {
			out << fixed << setprecision(8) << mat(i, j) << "\t";
		}
		out << "\n";
	}
	out.close();
	return;
};



// EOP partials from the appendix given in the handout
void Nx(const point& coords, const parameters& xhat, const MatrixXd& R, double& out) {

	out = R(0, 0) * (coords.x - xhat.x) + R(1, 0) * (coords.y - xhat.y) + R(2, 0) * (coords.z - xhat.z);
	return;
}
void Ny(const point& coords, const parameters& xhat, const MatrixXd& R, double& out) {

	out = R(0, 1) * (coords.x - xhat.x) + R(1, 1) * (coords.y - xhat.y) + R(2, 1) * (coords.z - xhat.z);
	return;
}
void D(const point& coords, const parameters& xhat, const MatrixXd& R, double& out) {

	out = R(0, 2) * (coords.x - xhat.x) + R(1, 2) * (coords.y - xhat.y) + R(2, 2) * (coords.z - xhat.z);
	return;
};
void dxdx(const point& coords, const parameters& xhat, const MatrixXd& R, const double& U, const double& W, const double& C, double& out) {

	out = C * (R(0, 0) * W - R(0, 2) * U) / pow(W, 2);
	return;
}
void dxdy(const point& coords, const parameters& xhat, const MatrixXd& R, const double& U, const double& W, const double& C, double& out) {

	out = C * (R(1, 0) * W - R(1, 2) * U) / pow(W, 2);
	return;
}
void dxdz(const point& coords, const parameters& xhat, const MatrixXd& R, const double& U, const double& W, const double& C, double& out) {

	out = C * (R(2, 0) * W - R(2, 2) * U) / pow(W, 2);
	return;
}
void dxdomega(const point& coords, const parameters& xhat, const MatrixXd& R, const double& U, const double& W, const double& C, double& out) {

	out = (-C / pow(W, 2)) * (W * (-R(2, 0) * (coords.y - xhat.y) + R(1, 0) * (coords.z - xhat.z))
		+ U * (R(2, 2) * (coords.y - xhat.y) - R(1, 2) * (coords.z - xhat.z)));
	return;
}
void dxdphi(const point& coords, const parameters& xhat, const double& U, const double& V, const double& W, const double& C, double& out) {

	out = (-C / pow(W, 2)) * (-pow(W, 2) * cos(xhat.kappa) + U * (-U * cos(xhat.kappa) + V * sin(xhat.kappa)));
	return;
}
void dxdkappa(const double& V, const double& W, const double& C, double& out) {

	out = -C * (V / W);
	return;
}
void dydx(const point& coords, const parameters& xhat, const MatrixXd& R, const double& V, const double& W, const double& C, double& out) {

	out = C * (R(0, 1) * W - R(0, 2) * V) / pow(W, 2);
	return;
}
void dydy(const point& coords, const parameters& xhat, const MatrixXd& R, const double& V, const double& W, const double& C, double& out) {

	out = C * (R(1, 1) * W - R(1, 2) * V) / pow(W, 2);
	return;
}
void dydz(const point& coords, const parameters& xhat, const MatrixXd& R, const double& V, const double& W, const double& C, double& out) {

	out = C * (R(2, 1) * W - R(2, 2) * V) / pow(W, 2);
	return;
}
void dydomega(const point& coords, const parameters& xhat, const MatrixXd& R, const double& V, const double& W, const double& C, double& out) {

	out = (-C / pow(W, 2)) * (W * (-R(2, 1) * (coords.y - xhat.y) + R(1, 1) * (coords.z - xhat.z))
		+ V * (R(2, 2) * (coords.y - xhat.y) - R(1, 2) * (coords.z - xhat.z)));
	return;
}
void dydphi(const point& coords, const parameters& xhat, const double& U, const double& V, const double& W, const double& C, double& out) {

	out = (-C / pow(W, 2)) * (pow(W, 2) * sin(xhat.kappa) + V * (U * cos(xhat.kappa) + V * sin(xhat.kappa)));
	return;
}
void dydkappa(const double& U, const double& W, const double& C, double& out) {

	out = (C * U) / W;
	return;
}