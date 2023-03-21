#ifndef ENGINE_H
#define ENGINE_H

// Definitions
#include "myLib.h"
#include "MLayers.h"
#include "Utilities.h"
#include "QuadL.h"
#include "Bessel.h"
#include "Vector.h"

typedef struct Field Field;
struct Field{
	complex double x, y, z;
};

typedef struct PlaneWaveFields PlaneWaveFields;
struct PlaneWaveFields{
	Field E, H;
};

typedef struct FarField FarField;
struct FarField{
	complex double theta, phi;
};

typedef struct Dipole Dipole;
struct Dipole{
	complex double Ilx, Ily, Ilz;
	double x_, y_, z_;
	double theta0, phi0;
};

// Functions
complex double GTest_integrand(complex double t, void *args); // For testing only!
int selectLayer(double z, MLayers *myConfig); // For testing only!
double findMax(MLayers *myConfig); // For testing only!
double findDetourParameters(double rho, double z, double z_, MLayers *myConfig); // For testing only!

//
complex double Gamma_d(int n, char p, complex double k_rho, MLayers *myConfig);
complex double Gamma_u(int n, char p, complex double k_rho, MLayers *myConfig);
void TLGF(double z, double z_, complex double k_rho, char p, char s, int option, 
	MLayers *myConfig, complex double *V, complex double *I);

PlaneWaveFields PlaneWave(complex double ETheta, complex double EPhi,
	double theta_i, double phi_i, double x, double y, double z, MLayers *myConfig);
FarField FarFieldDipoleJ(Dipole J, double theta_s, double phi_s, double r, MLayers *myConfig);
FarField FarFieldDipoleM(Dipole M, double theta_s, double phi_s, double r, MLayers *myConfig);

complex double GEJxx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEJxy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEJxz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEJyx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEJyy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEJyz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEJzx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEJzy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEJzz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);

complex double GEMxx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEMxy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEMxz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEMyx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEMyy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEMyz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEMzx(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEMzy(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);
complex double GEMzz(double x, double x_, double y, double y_, double z, double z_, 
	double tol, int kmax, MLayers *myConfig);

#endif