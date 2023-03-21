//
#include "TestBench.h"
#include "Utilities.h"
#include "QuadL.h"
#include "MLayers.h"
#include "Engine.h"

void TestUtilities(){

	complex double j=csqrt(-1.0);
	complex double z=2.4+j*0.7;
	showComplex(z);

	double x=roundn(0.123456789, 4);
	printf("%21.14E\n", x);

	int N=1000;
	Timer T;
	setTimer(&T);
	for (int i=0; i<N; i++){
		usleep(1000.0);
		progressBar(i, N, "");
	}
	unsetTimer(&T);

	printf("%21.14E\n", deg2rad(180.0));
	printf("%21.14E\n", rad2deg(pi));

}

complex double func(complex double z, void *args){
	assert(args==NULL);
	double beta=10.0;
    return cexp(-beta*beta*z*z)+csin(z);
}

void TestQuadL(){

	double a=-2.0;
	double b=+4.0;
	double tol=1.0E-14;
	int kmax=10;

	int flag=0;
	complex double ans=0.0;
	ans = QuadLAdaptive1D_32(func, NULL, a, b, 
		tol, kmax, 0, 0, &flag);
	// 4.14742169407021E-01
	showComplex(ans);
	printf("flag = %d\n", flag);

}

void TestMLayer(){

	// Allocate The Memory
	int N=4;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigPaulus(&myConfig);
	logConfig(&myConfig);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestGammaDown(){

	// Allocate The Memory
	int N=3;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	// ConfigPaulus(&myConfig);
	// ConfigGoldKretschmann(633.0E-9, &myConfig);
	ConfigOttoGraphene(1.0E12, &myConfig);
	logConfig(&myConfig);

	// Test Gamma
	int Ns=1001;
	double theta_min=deg2rad(0.0);
	double theta_max=deg2rad(90.0);
	double dtheta = (theta_max-theta_min)/(Ns-1.0);
	double theta=0.0;
	complex double Gamma=0.0;
	complex double k_rho=0.;
	FILE *file=fopen("Data/Reflection/Gamma.dat", "w");
	for (int i=0; i<Ns; i++){
		theta = theta_min+i*dtheta;
		k_rho = myConfig.Layers[0].k*sin(theta);
		assert(fprintf(file, "%21.14E ", rad2deg(theta)));
		Gamma = Gamma_d(1, 'e', k_rho, &myConfig);
		assert(fprintf(file,"%21.14E ", cabs(Gamma)*cabs(Gamma)));
		Gamma = Gamma_d(1, 'h', k_rho, &myConfig);
		assert(fprintf(file,"%21.14E ", cabs(Gamma)*cabs(Gamma)));
		assert(fprintf(file,"\n"));
	}
	fclose(file);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestGammaUp(){

	// Allocate The Memory
	int N=3;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	// ConfigPaulus(&myConfig);
	// ConfigGoldKretschmann(633.0E-9, &myConfig);
	ConfigOttoGraphene(1.0E12, &myConfig);
	logConfig(&myConfig);

	// Test Gamma
	int Ns=1001;
	double theta_min=deg2rad(0.0);
	double theta_max=deg2rad(90.0);
	double dtheta = (theta_max-theta_min)/(Ns-1.0);
	double theta=0.0;
	complex double Gamma=0.0;
	complex double k_rho=0.;
	FILE *file=fopen("Data/Reflection/Gamma.dat", "w");
	for (int i=0; i<Ns; i++){
		theta = theta_min+i*dtheta;
		k_rho = myConfig.Layers[N-1].k*sin(theta);
		assert(fprintf(file, "%21.14E ", rad2deg(theta)));
		Gamma = Gamma_u(N, 'e', k_rho, &myConfig);
		assert(fprintf(file,"%21.14E ", cabs(Gamma)*cabs(Gamma)));
		Gamma = Gamma_u(N, 'h', k_rho, &myConfig);
		assert(fprintf(file,"%21.14E ", cabs(Gamma)*cabs(Gamma)));
		assert(fprintf(file,"\n"));
	}
	fclose(file);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestSelectLayer(){

	// Allocate The Memory
	double nm=1.0E-9;
	int N=4;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigPaulus(&myConfig);
	logConfig(&myConfig);

	printf("%d\n", selectLayer(+7500*nm, &myConfig));
	printf("%d\n", selectLayer(+750*nm, &myConfig));
	printf("%d\n", selectLayer(+200*nm, &myConfig));
	printf("%d\n", selectLayer(+0*nm, &myConfig));
	printf("%d\n", selectLayer(-250*nm, &myConfig));
	printf("%d\n", selectLayer(-750*nm, &myConfig));
	printf("%d\n", selectLayer(-7500*nm, &myConfig));
	
	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestTLGF(){

	// Allocate The Memory
	double nm=1.0E-9;
	int N=4;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigPaulus(&myConfig);
	logConfig(&myConfig);

	int option=0;
	int Ns=1000;
	double z_=+750*nm;
	complex double k_rho=0.7*myConfig.k0;
	double z_min=-1000.0*nm;
	double z_max=+1000.0*nm;
	double z, dz=(z_max-z_min)/(Ns-1.0);
	complex double V, I;
	
	FILE *file=fopen("Data/TLGF/Data.dat", "w");
	for (int i=0; i<Ns; i++){
		z = z_min+i*dz;
		assert(fprintf(file, "%21.14E ", z/nm));
		TLGF(z, z_, k_rho, 'e', 'v', option, &myConfig, &V, &I);
		assert(fprintf(file, "%21.14E %21.14E ", creal(V), cimag(V)));
		assert(fprintf(file, "%21.14E %21.14E ", creal(I)*eta0, cimag(I)*eta0));
		TLGF(z, z_, k_rho, 'e', 'i', option, &myConfig, &V, &I);
		assert(fprintf(file, "%21.14E %21.14E ", creal(V)/eta0, cimag(V)/eta0));
		assert(fprintf(file, "%21.14E %21.14E ", creal(I), cimag(I)));
		TLGF(z, z_, k_rho, 'h', 'v', option, &myConfig, &V, &I);
		assert(fprintf(file, "%21.14E %21.14E ", creal(V), cimag(V)));
		assert(fprintf(file, "%21.14E %21.14E ", creal(I)*eta0, cimag(I)*eta0));
		TLGF(z, z_, k_rho, 'h', 'i', option, &myConfig, &V, &I);
		assert(fprintf(file, "%21.14E %21.14E ", creal(V)/eta0, cimag(V)/eta0));
		assert(fprintf(file, "%21.14E %21.14E ", creal(I), cimag(I)));
		assert(fprintf(file,"\n"));
	}
	fclose(file);
	
	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestSI(){

	// Allocate The Memory
	double nm=1.0E-9;
	int N=4;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigPaulus(&myConfig);
	logConfig(&myConfig);

	int Ns=1000;
	double x_=0.0*nm;
	double y_=0.0*nm;
	double z_=+750*nm;
	double z_min=-1000.0*nm;
	double z_max=+1000.0*nm;
	double rho=myConfig.lambda0;
	double phi=deg2rad(45.0);
	double tol=1.0E-4;
	int kmax=10;

	double x=rho*cos(phi);
	double y=rho*sin(phi);
	double z, dz=(z_max-z_min)/(Ns-1.0);

	FILE *file=fopen("Data/SI/Data.dat", "w");
	complex double E;

	Timer T;
	setTimer(&T);
	for (int i=0; i<Ns; i++){
		z = z_min+i*dz;
		fprintf(file, "%21.14E ", z/nm);

		E = GEJxx(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEJxy(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEJxz(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEJyx(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEJyy(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEJyz(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEJzx(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEJzy(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEJzz(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));

		E = GEMxx(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEMxy(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEMxz(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEMyx(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEMyy(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEMyz(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEMzx(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEMzy(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		E = GEMzz(x, x_, y, y_, z, z_, tol, kmax, &myConfig);
		fprintf(file, "%21.14E ", isinf(log10(cabs(E))) ? 0.0 : log10(cabs(E)));
		fprintf(file, "\n");

		progressBar(i, Ns, "");
	}
	unsetTimer(&T);
	fclose(file);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestChew(){

	// Allocate The Memory
	int N=7;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigChew(&myConfig);
	logConfig(&myConfig);

	int Ns=1001;
	double x_=0.0;
	double y_=0.0;
	double z_=-1.4;
	double theta0=deg2rad(20.0);
	double phi0=deg2rad(30.0);
	double x_min=-3.0;
	double x_max=+3.0;
	double y=+1.0;
	double z=-0.3;
	double tol=1.0E-4;
	int kmax=15;

	double x, dx=(x_max-x_min)/(Ns-1.0);
	complex double Jx, Jy, Jz;
	Jx = sin(theta0)*cos(phi0);
	Jy = sin(theta0)*sin(phi0);
	Jz = cos(theta0);
	complex double Mx, My, Mz;
	Mx = sin(theta0)*cos(phi0);
	My = sin(theta0)*sin(phi0);
	Mz = cos(theta0);

	Timer T;
	complex double Ex, Ey, Ez;

	FILE *file1=fopen("Data/Chew/DataEJ.dat", "w");
	setTimer(&T);
	for (int i=0; i<Ns; i++){
		x = x_min+i*dx;
		fprintf(file1, "%21.14E ", x);
		Ex = GEJxx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jx+
			 GEJxy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jy+
			 GEJxz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jz;
		Ey = GEJyx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jx+
			 GEJyy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jy+
			 GEJyz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jz;
		Ez = GEJzx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jx+
			 GEJzy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jy+
			 GEJzz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jz;
		fprintf(file1, "%21.14E ", cabs(Ex));
		fprintf(file1, "%21.14E ", cabs(Ey));
		fprintf(file1, "%21.14E ", cabs(Ez));
		fprintf(file1, "\n");
		progressBar(i, Ns, "");
	}
	unsetTimer(&T);
	fclose(file1);

	FILE *file2=fopen("Data/Chew/DataEM.dat", "w");
	setTimer(&T);
	for (int i=0; i<Ns; i++){
		x = x_min+i*dx;
		fprintf(file2, "%21.14E ", x);
		Ex = GEMxx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Mx+
			 GEMxy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*My+
			 GEMxz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Mz;
		Ey = GEMyx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Mx+
			 GEMyy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*My+
			 GEMyz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Mz;
		Ez = GEMzx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Mx+
			 GEMzy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*My+
			 GEMzz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Mz;
		fprintf(file2, "%21.14E ", cabs(Ex));
		fprintf(file2, "%21.14E ", cabs(Ey));
		fprintf(file2, "%21.14E ", cabs(Ez));
		fprintf(file2, "\n");
		progressBar(i, Ns, "");
	}
	unsetTimer(&T);
	fclose(file2);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestNadson(){

	// Definitions
	double um=1.0E-6;
	double nm=1.0E-9;
	double lambda0=633.0*nm;

	// Allocate The Memory
	int N=3;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigGoldKretschmannNadson(lambda0, &myConfig);
	logConfig(&myConfig);

	int Ns=1001;
	double x_=0.0;
	double y_=0.0;
	double z_=+20.0*nm;
	double theta0=deg2rad(0.0);
	double phi0=deg2rad(0.0);
	double x_min=-2000.0*nm;
	double x_max=+2000.0*nm;
	double y=0.0*nm;
	double z_min=-2000.0*nm;
	double z_max=+2000.0*nm;
	double tol=1.0E-4;
	int kmax=15;

	double x, dx=(x_max-x_min)/(Ns-1.0);
	double z, dz=(z_max-z_min)/(Ns-1.0);
	complex double Jx, Jy, Jz;
	Jx = sin(theta0)*cos(phi0);
	Jy = sin(theta0)*sin(phi0);
	Jz = cos(theta0);

	Timer T;
	complex double Ex, Ey, Ez;
	double E;

	FILE *file=fopen("Data/Nadson/Data.dat", "w");
	setTimer(&T);
	for (int m=0; m<Ns; m++){
		x = x_min+m*dx;
		for (int n=0; n<Ns; n++){
			z = z_min+n*dz;
			Ex = GEJxx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jx+
				 GEJxy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jy+
				 GEJxz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jz;
			Ey = GEJyx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jx+
				 GEJyy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jy+
				 GEJyz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jz;
			Ez = GEJzx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jx+
				 GEJzy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jy+
				 GEJzz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jz;
			E = sqrt(creal(Ex)*creal(Ex)+creal(Ey)*creal(Ey)+creal(Ez)*creal(Ez));
			fprintf(file, "%21.14E %21.14E %21.14E\n", x/um, z/um, 20*log10(E));
		}
		fprintf(file, "\n");
		progressBar(m, Ns, "");
	}
	unsetTimer(&T);
	fclose(file);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestPaulusNearField(){

	// Definitions
	double um=1.0E-6;
	double nm=1.0E-9;

	// Allocate The Memory
	int N=4;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigPaulus(&myConfig);
	logConfig(&myConfig);

	int Ns=1001;
	double x_=0.0;
	double y_=0.0;
	double z_=+750.0*nm;
	double theta0=deg2rad(0.0);
	double phi0=deg2rad(0.0);
	double x_min=-2000.0*nm;
	double x_max=+2000.0*nm;
	double y=0.0*nm;
	double z_min=-2000.0*nm;
	double z_max=+2000.0*nm;
	double tol=1.0E-4;
	int kmax=15;

	double x, dx=(x_max-x_min)/(Ns-1.0);
	double z, dz=(z_max-z_min)/(Ns-1.0);
	complex double Jx, Jy, Jz;
	Jx = sin(theta0)*cos(phi0);
	Jy = sin(theta0)*sin(phi0);
	Jz = cos(theta0);

	Timer T;
	complex double Ex, Ey, Ez;
	double E;

	FILE *file=fopen("Data/Paulus/Data.dat", "w");
	setTimer(&T);
	for (int m=0; m<Ns; m++){
		x = x_min+m*dx;
		for (int n=0; n<Ns; n++){
			z = z_min+n*dz;
			Ex = GEJxx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jx+
				 GEJxy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jy+
				 GEJxz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jz;
			Ey = GEJyx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jx+
				 GEJyy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jy+
				 GEJyz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jz;
			Ez = GEJzx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jx+
				 GEJzy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jy+
				 GEJzz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*Jz;
			E = sqrt(creal(Ex)*creal(Ex)+creal(Ey)*creal(Ey)+creal(Ez)*creal(Ez));
			fprintf(file, "%21.14E %21.14E %21.14E\n", x/um, z/um, 20*log10(E));
		}
		fprintf(file, "\n");
		progressBar(m, Ns, "");
	}
	unsetTimer(&T);
	fclose(file);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestPlaneWave(){

	// Definitions
	double um=1.0E-6;
	double nm=1.0E-9;
	double lambda0=633.0*nm;

	// Allocate The Memory
	int N=4;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigPlasmonicWG(lambda0, &myConfig);
	logConfig(&myConfig);

	int Ns=1001;
	complex double ETheta=1.0;
	complex double EPhi=0.0;
	double theta_i=deg2rad(71.08);
	double phi_i=deg2rad(180.0);
	double x_min=-800.0*nm;
	double x_max=+800.0*nm;
	double y=0.0*nm;
	double z_min=-600.0*nm;
	double z_max=+1000.0*nm;

	double x, dx=(x_max-x_min)/(Ns-1.0);
	double z, dz=(z_max-z_min)/(Ns-1.0);

	Timer T;
	complex double Ex, Ey, Ez;
	double E;

	PlaneWaveFields results;
	FILE *file=fopen("Data/PlaneWave/Data.dat", "w");
	setTimer(&T);
	for (int m=0; m<Ns; m++){
		x = x_min+m*dx;
		for (int n=0; n<Ns; n++){
			z = z_min+n*dz;
			results = PlaneWave(ETheta, EPhi, theta_i, phi_i, x, y, z, &myConfig);
			Ex = results.E.x;
			Ey = results.E.y;
			Ez = results.E.z;
			E = sqrt(creal(Ex)*creal(Ex)+creal(Ey)*creal(Ey)+creal(Ez)*creal(Ez));
			fprintf(file, "%21.14E %21.14E %21.14E\n", x/um, z/um, E);
		}
		fprintf(file, "\n");
		progressBar(m, Ns, "");
	}
	unsetTimer(&T);
	fclose(file);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestFEKOPlaneWave(){

	// Allocate The Memory
	int N=3;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigFEKOPlaneWave(&myConfig);
	logConfig(&myConfig);

	int Ns=1001;
	complex double ETheta;
	complex double EPhi;
	double theta_i=deg2rad(30.0);
	double phi_i=deg2rad(60.0);
	double x=0.0;
	double y=0.0;
	double z_min=-2.0;
	double z_max=+2.0;

	double z, dz=(z_max-z_min)/(Ns-1.0);

	Timer T;
	PlaneWaveFields results;
	complex double Ex, Ey, Ez;
	complex double Hx, Hy, Hz;
	double E, H;

	// TM
	ETheta = 1.0;
	EPhi = 0.0;
	FILE *fileTM=fopen("Data/FEKOPlaneWave/DataTM.dat", "w");
	setTimer(&T);
	for (int i=0; i<Ns; i++){
		z = z_min+i*dz;
		results = PlaneWave(ETheta, EPhi, theta_i, phi_i, x, y, z, &myConfig);
		Ex = results.E.x;
		Ey = results.E.y;
		Ez = results.E.z;
		E = sqrt(cabs(Ex)*cabs(Ex)+cabs(Ey)*cabs(Ey)+cabs(Ez)*cabs(Ez));
		Hx = results.H.x;
		Hy = results.H.y;
		Hz = results.H.z;
		H = sqrt(cabs(Hx)*cabs(Hx)+cabs(Hy)*cabs(Hy)+cabs(Hz)*cabs(Hz));
		fprintf(fileTM, "%21.14E %21.14E %21.14E\n", z, E, H);
		progressBar(i, Ns, "");
	}
	unsetTimer(&T);
	fclose(fileTM);

	// TE
	ETheta = 0.0;
	EPhi = 1.0;
	FILE *fileTE=fopen("Data/FEKOPlaneWave/DataTE.dat", "w");
	setTimer(&T);
	for (int i=0; i<Ns; i++){
		z = z_min+i*dz;
		results = PlaneWave(ETheta, EPhi, theta_i, phi_i, x, y, z, &myConfig);
		Ex = results.E.x;
		Ey = results.E.y;
		Ez = results.E.z;
		E = sqrt(cabs(Ex)*cabs(Ex)+cabs(Ey)*cabs(Ey)+cabs(Ez)*cabs(Ez));
		Hx = results.H.x;
		Hy = results.H.y;
		Hz = results.H.z;
		H = sqrt(cabs(Hx)*cabs(Hx)+cabs(Hy)*cabs(Hy)+cabs(Hz)*cabs(Hz));
		fprintf(fileTM, "%21.14E %21.14E %21.14E\n", z, E, H);
		progressBar(i, Ns, "");
	}
	unsetTimer(&T);
	fclose(fileTE);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestFarField(){

	// Allocate The Memory
	int N=3;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigFEKOPlaneWave(&myConfig);
	logConfig(&myConfig);
	
	int Ns=1001;
	complex double theta0=deg2rad(90.0);
	complex double phi0=deg2rad(0.0);
	double x_=+0.0;
	double y_=+0.0;
	double z_=-0.5;
	complex double Il=1.0;
	Dipole J={Il*sin(theta0)*cos(phi0), 
			  Il*sin(theta0)*sin(phi0), 
			  Il*cos(theta0), 
			  x_, y_, z_, theta0, phi0};
	double r=100.0*myConfig.lambda0;
	double phi=deg2rad(90.0);
	double theta, dtheta=2.0*pi/(Ns-1.0);

	Timer T;
	FarField results;
	complex double E_theta, E_phi;

	FILE *file=fopen("Data/FarField/Data.dat", "w");
	setTimer(&T);
	for (int i=0; i<Ns; i++){
		theta = i*dtheta;
		if (theta<=pi){
			results = FarFieldDipoleJ(J, theta, phi, r, &myConfig);
		}else
		if ((theta>pi)&&(theta<=2.0*pi)){
			results = FarFieldDipoleJ(J, 2.0*pi-theta, phi+pi, r, &myConfig);
		}
		E_theta = results.theta;
		E_phi = results.phi;
		fprintf(file, "%21.14E %21.14E %21.14E %21.14E\n", theta, phi, cabs(E_theta), cabs(E_phi));
		progressBar(i, Ns, "");
	}
	unsetTimer(&T);
	fclose(file);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestFarFieldBragg(){

	// Definitions
	double nm=1.0E-9;
	// Allocate The Memory
	int N=19;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigBragg(&myConfig);
	logConfig(&myConfig);
	
	int Ns=10001;
	complex double theta0=deg2rad(135.0);
	complex double phi0=deg2rad(0.0);
	double x_=+0.0*nm;
	double y_=+0.0*nm;
	double z_=-15.0*nm;
	complex double Il=1.0;
	Dipole J={Il*sin(theta0)*cos(phi0), 
			  Il*sin(theta0)*sin(phi0), 
			  Il*cos(theta0), 
			  x_, y_, z_, theta0, phi0};
	double r=100.0*myConfig.lambda0;
	double phi=deg2rad(0.0);
	double theta, dtheta=2.0*pi/(Ns-1.0);

	Timer T;
	FarField results;
	complex double E_theta, E_phi;

	FILE *file=fopen("Data/Bragg/Data.dat", "w");
	setTimer(&T);
	for (int i=0; i<Ns; i++){
		theta = i*dtheta;
		if (theta<=pi){
			results = FarFieldDipoleJ(J, theta, phi, r, &myConfig);
		}else
		if ((theta>pi)&&(theta<=2.0*pi)){
			results = FarFieldDipoleJ(J, 2.0*pi-theta, phi+pi, r, &myConfig);
		}
		E_theta = results.theta;
		E_phi = results.phi;
		fprintf(file, "%21.14E %21.14E %21.14E %21.14E\n", theta, phi, cabs(E_theta), cabs(E_phi));
		progressBar(i, Ns, "");
	}
	unsetTimer(&T);
	fclose(file);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestFarFieldvsExact(){

	// Allocate The Memory
	int N=3;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	ConfigFEKOPlaneWave(&myConfig);
	// ConfigPaulus(&myConfig);
	logConfig(&myConfig);
	
	int Ns=1001;
	complex double theta0=deg2rad(0.0);
	complex double phi0=deg2rad(0.0);
	double x_=+0.0;
	double y_=+0.0;
	double z_=-0.05;
	complex double Il=1.0;
	Dipole J={Il*sin(theta0)*cos(phi0), 
			  Il*sin(theta0)*sin(phi0), 
			  Il*cos(theta0), 
			  x_, y_, z_, theta0, phi0};
	double r=10.0*myConfig.lambda0;
	double phi=deg2rad(0.0);
	double theta, dtheta=2.0*pi/(Ns-1.0);
	double x, y, z;

	int kmax=25;
	double tol=1.0E-6;

	Timer T;
	FarField results;
	complex double E_theta, E_phi;
	complex double Ex, Ey, Ez;
	double E_far, E_exact;

	FILE *file=fopen("Data/FarField/DataVsExact.dat", "w");
	setTimer(&T);
	for (int i=0; i<Ns; i++){
		results.theta = 0.0;
		results.phi = 0.0;

		theta = i*dtheta;
		if (theta<=pi){
			results = FarFieldDipoleJ(J, theta, phi, r, &myConfig);
			// results = FarFieldDipoleM(J, theta, phi, r, &myConfig);
		}else
		if ((theta>pi)&&(theta<=2.0*pi)){
			results = FarFieldDipoleJ(J, 2.0*pi-theta, phi+pi, r, &myConfig);
			// results = FarFieldDipoleM(J, 2.0*pi-theta, phi+pi, r, &myConfig);
		}

		E_theta = results.theta;
		E_phi = results.phi;
		E_far = sqrt(cabs(E_theta)*cabs(E_theta)+cabs(E_phi)*cabs(E_phi));
		fprintf(file, "%21.14E %21.14E %21.14E ", theta, phi, E_far);

		x = r*sin(theta)*cos(phi);
		y = r*sin(theta)*sin(phi);
		z = r*cos(theta);

		Ex = GEJxx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilx+
			 GEJxy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ily+
			 GEJxz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilz;
		Ey = GEJyx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilx+
			 GEJyy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ily+
			 GEJyz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilz;
		Ez = GEJzx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilx+
			 GEJzy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ily+
			 GEJzz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilz;
		
		// Ex = GEMxx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilx+
		// 	 GEMxy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ily+
		// 	 GEMxz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilz;
		// Ey = GEMyx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilx+
		// 	 GEMyy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ily+
		// 	 GEMyz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilz;
		// Ez = GEMzx(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilx+
		// 	 GEMzy(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ily+
		// 	 GEMzz(x, x_, y, y_, z, z_, tol, kmax, &myConfig)*J.Ilz;

		E_exact = sqrt(creal(Ex)*creal(Ex)+creal(Ey)*creal(Ey)+creal(Ez)*creal(Ez));
		fprintf(file, "%21.14E\n", E_exact);

		progressBar(i, Ns, "");
	}
	unsetTimer(&T);
	fclose(file);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}

void TestIntegrand(){

	// Allocate The Memory
	int N=4;
	MLayers myConfig=DefaultMLayers;
	myConfig.N = N;
	myConfig.Layers = (Layer*)malloc(myConfig.N*sizeof(Layer));

	// ConfigFEKOPlaneWave(&myConfig);
	ConfigPaulus(&myConfig);
	logConfig(&myConfig);
	
	int Ns=10001;
	double lambda0=myConfig.lambda0;
	double k0=myConfig.k0;
	double r=10.0*lambda0;
	double theta=deg2rad(60.0);
	double z_=+250E-9;
	double z=z_+r*cos(theta);
	double rho=r*sin(theta);
	double k_rho_min=0.0*k0;
	double k_rho_max=4.0*k0;
	double a=findMax(&myConfig); 
	double b=findDetourParameters(rho, z, z_, &myConfig);
	printf("%21.14E %21.14E\n", a/k0, b/k0);
	double k_rho, d_k_rho=(k_rho_max-k_rho_min)/(Ns-1.0);
	typedef struct Args Args;
	struct Args{
		double z, z_, rho;
		double a, b;
		MLayers *myConfig;
	};
	Args myArgs={z, z_, rho, a, b, &myConfig};
	complex double integrand;
	FILE *file=fopen("Data/SI_Integrand/Data.dat", "w");
	for (int i=0; i<Ns; i++){
		k_rho = k_rho_min+i*d_k_rho;
		integrand = GTest_integrand(k_rho, &myArgs);
		fprintf(file, "%21.14E %21.14E\n", k_rho/k0, cabs(integrand));
	}
	fclose(file);

	// Free The Memory
	free(myConfig.Layers);
	myConfig.Layers = NULL;
}


