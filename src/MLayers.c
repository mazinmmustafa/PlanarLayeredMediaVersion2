//
#include "MLayers.h"

void logConfig(MLayers *myConfig){	
	assert(myConfig->Layers!=NULL);
	assert(myConfig->N>0);
	FILE *file=fopen("Data/log.txt", "w");
	assert(file!=NULL);
	assert(fprintf(file, "N = %d\n", myConfig->N));
	assert(fprintf(file, "lambda0 = %21.14E\n", myConfig->lambda0));
	assert(fprintf(file, "freq = %21.14E\n", myConfig->freq));
	assert(fprintf(file, "k0 = %21.14E\n", myConfig->k0));
	assert(fprintf(file, "Gamma_u = (%21.14E, %21.14E)\n", creal(myConfig->Gamma_u), cimag(myConfig->Gamma_u)));
	assert(fprintf(file, "Gamma_d = (%21.14E, %21.14E)\n", creal(myConfig->Gamma_d), cimag(myConfig->Gamma_d)));
	assert(fprintf(file, "\n"));
	for (int n=0; n<myConfig->N; n++){
		assert(fprintf(file, "Layer %d\n", n+1));
		assert(fprintf(file, "zu = %21.14E\n", myConfig->Layers[n].zu));
		assert(fprintf(file, "zd = %21.14E\n", myConfig->Layers[n].zd));
		assert(fprintf(file, "d = %21.14E\n", myConfig->Layers[n].d));
		assert(fprintf(file, "eps = (%21.14E, %21.14E)\n", creal(myConfig->Layers[n].eps), cimag(myConfig->Layers[n].eps)));
		assert(fprintf(file, "mu = (%21.14E, %21.14E)\n", creal(myConfig->Layers[n].mu), cimag(myConfig->Layers[n].mu)));
		assert(fprintf(file, "sigma_s = (%21.14E, %21.14E)\n", creal(myConfig->Layers[n].sigma_s), cimag(myConfig->Layers[n].sigma_s)));
		assert(fprintf(file, "k = (%21.14E, %21.14E)\n", creal(myConfig->Layers[n].k), cimag(myConfig->Layers[n].k)));
		assert(fprintf(file, "eta = (%21.14E, %21.14E)\n", creal(myConfig->Layers[n].eta), cimag(myConfig->Layers[n].eta)));
		assert(fprintf(file, "\n"));
	}
	fclose(file);
}

void ConfigPaulus(MLayers *myConfig){
	assert(myConfig->N==4);
	// Definitions
	double nm=1.0E-9;
	complex double j=csqrt(-1.0);
	// General Parameters
	myConfig->lambda0 = 633.0*nm;
	myConfig->freq = c0/myConfig->lambda0;
	myConfig->k0 = 2.0*pi/myConfig->lambda0;
	myConfig->Gamma_u = 0.0;
	myConfig->Gamma_d = 0.0;
	int n=0;
	// Layer 1
	myConfig->Layers[n].zu = +2000.0*nm;
	myConfig->Layers[n].zd = +500.0*nm;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 2
	myConfig->Layers[n].zu = +500.0*nm;
	myConfig->Layers[n].zd = +0.0*nm;
	myConfig->Layers[n].eps = 2.0-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 3
	myConfig->Layers[n].zu = +0.0*nm;
	myConfig->Layers[n].zd = -500.0*nm;
	myConfig->Layers[n].eps = 10.0-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 4
	myConfig->Layers[n].zu = -500.0*nm;
	myConfig->Layers[n].zd = -2000.0*nm;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
}

void ConfigGoldKretschmann(double lambda0, MLayers *myConfig){
	assert(myConfig->N==3);
	// Definitions
	double nm=1.0E-9;
	complex double j=csqrt(-1.0);
	// General Parameters
	myConfig->lambda0 = lambda0;
	myConfig->freq = c0/myConfig->lambda0;
	myConfig->k0 = 2.0*pi/myConfig->lambda0;
	myConfig->Gamma_u = 0.0;
	myConfig->Gamma_d = 0.0;
	int n=0;
	// Layer 1
	myConfig->Layers[n].zu = +2000.0*nm;
	myConfig->Layers[n].zd = +0.0*nm;
	myConfig->Layers[n].eps = 2.3013-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 2
	myConfig->Layers[n].zu = +0.0*nm;
	myConfig->Layers[n].zd = -50.0*nm;
	myConfig->Layers[n].eps = -11.753-j*1.2596;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 3
	myConfig->Layers[n].zu = -50.0*nm;
	myConfig->Layers[n].zd = -2000.0*nm;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
}

void ConfigOttoGraphene(double freq, MLayers *myConfig){
	assert(myConfig->N==3);
	// Definitions
	double um=1.0E-6;
	complex double j=csqrt(-1.0);
	// General Parameters
	myConfig->lambda0 = c0/freq;
	myConfig->freq = c0/myConfig->lambda0;
	myConfig->k0 = 2.0*pi/myConfig->lambda0;
	myConfig->Gamma_u = 0.0;
	myConfig->Gamma_d = 0.0;
	int n=0;
	// Layer 1
	myConfig->Layers[n].zu = +1000.0*um;
	myConfig->Layers[n].zd = +0.0*um;
	myConfig->Layers[n].eps = (2.003-j*0.0)*(2.003-j*0.0);
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 2
	myConfig->Layers[n].zu = +0.0*um;
	myConfig->Layers[n].zd = -20.0*um;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = +0.000369059545723-j*0.015237384931248;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 3
	myConfig->Layers[n].zu = -20.0*um;
	myConfig->Layers[n].zd = -1000.0*um;
	myConfig->Layers[n].eps = (1.762-j*0.0)*(1.762-j*0.0);
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
}

void ConfigChew(MLayers *myConfig){
	assert(myConfig->N==7);
	// Definitions
	complex double j=csqrt(-1.0);
	// General Parameters
	myConfig->lambda0 = 1.0;
	myConfig->freq = c0/myConfig->lambda0;
	myConfig->k0 = 2.0*pi/myConfig->lambda0;
	myConfig->Gamma_u = 0.0;
	myConfig->Gamma_d = 0.0;
	int n=0;
	// Layer 1
	myConfig->Layers[n].zu = +10.0;
	myConfig->Layers[n].zd = +0.0;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 2
	myConfig->Layers[n].zu = +0.0;
	myConfig->Layers[n].zd = -0.2;
	myConfig->Layers[n].eps = 2.6-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 3
	myConfig->Layers[n].zu = -0.2;
	myConfig->Layers[n].zd = -0.5;
	myConfig->Layers[n].eps = 6.5-j*0.0;
	myConfig->Layers[n].mu = 3.2-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 4
	myConfig->Layers[n].zu = -0.5;
	myConfig->Layers[n].zd = -1.0;
	myConfig->Layers[n].eps = 4.2-j*0.0;
	myConfig->Layers[n].mu = 6.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 5
	myConfig->Layers[n].zu = -1.0;
	myConfig->Layers[n].zd = -1.3;
	myConfig->Layers[n].eps = 6.5-j*0.0;
	myConfig->Layers[n].mu = 3.2-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 6
	myConfig->Layers[n].zu = -1.3;
	myConfig->Layers[n].zd = -1.5;
	myConfig->Layers[n].eps = 2.6-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 7
	myConfig->Layers[n].zu = -1.5;
	myConfig->Layers[n].zd = -10.0;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
}

void ConfigGoldKretschmannNadson(double lambda0, MLayers *myConfig){
	assert(myConfig->N==3);
	// Definitions
	double nm=1.0E-9;
	complex double j=csqrt(-1.0);
	// General Parameters
	myConfig->lambda0 = lambda0;
	myConfig->freq = c0/myConfig->lambda0;
	myConfig->k0 = 2.0*pi/myConfig->lambda0;
	myConfig->Gamma_u = 0.0;
	myConfig->Gamma_d = 0.0;
	int n=0;
	// Layer 1
	myConfig->Layers[n].zu = +2000.0*nm;
	myConfig->Layers[n].zd = +0.0*nm;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 2
	myConfig->Layers[n].zu = +0.0*nm;
	myConfig->Layers[n].zd = -50.0*nm;
	myConfig->Layers[n].eps = -11.753-j*1.2596;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 3
	myConfig->Layers[n].zu = -50.0*nm;
	myConfig->Layers[n].zd = -2000.0*nm;
	myConfig->Layers[n].eps = 2.3013-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
}

void ConfigPlasmonicWG(double lambda0, MLayers *myConfig){
	assert(myConfig->N==4);
	// Definitions
	double nm=1.0E-9;
	complex double j=csqrt(-1.0);
	// General Parameters
	myConfig->lambda0 = lambda0;
	myConfig->freq = c0/myConfig->lambda0;
	myConfig->k0 = 2.0*pi/myConfig->lambda0;
	myConfig->Gamma_u = 0.0;
	myConfig->Gamma_d = 0.0;
	int n=0;
	// Layer 1
	myConfig->Layers[n].zu = +2000.0*nm;
	myConfig->Layers[n].zd = +50.0*nm;
	myConfig->Layers[n].eps = 2.3013-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 2
	myConfig->Layers[n].zu = +50.0*nm;
	myConfig->Layers[n].zd = +0.0*nm;
	myConfig->Layers[n].eps = -18.5698-j*0.397;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 3
	myConfig->Layers[n].zu = +0.0*nm;
	myConfig->Layers[n].zd = -30.0*nm;
	myConfig->Layers[n].eps = 2.07-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 4
	myConfig->Layers[n].zu = -30.0*nm;
	myConfig->Layers[n].zd = -2000.0*nm;
	myConfig->Layers[n].eps = 1.79-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
}

void ConfigFEKOPlaneWave(MLayers *myConfig){
	assert(myConfig->N==3);
	// Definitions
	complex double j=csqrt(-1.0);
	// General Parameters
	myConfig->lambda0 = 1.0;
	myConfig->freq = c0/myConfig->lambda0;
	myConfig->k0 = 2.0*pi/myConfig->lambda0;
	myConfig->Gamma_u = 0.0;
	myConfig->Gamma_d = 0.0;
	int n=0;
	// Layer 1
	myConfig->Layers[n].zu = +1.0;
	myConfig->Layers[n].zd = +0.0;
	myConfig->Layers[n].eps = 3.0-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 2
	myConfig->Layers[n].zu = +0.0;
	myConfig->Layers[n].zd = -0.1;
	myConfig->Layers[n].eps = 2.0-j*0.0;
	myConfig->Layers[n].mu = 2.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 3
	myConfig->Layers[n].zu = -0.1;
	myConfig->Layers[n].zd = -1.0;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].mu = 3.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
}

void ConfigBragg(MLayers *myConfig){
	assert(myConfig->N==19);
	// Definitions
	double nm=1.0E-9;
	complex double j=csqrt(-1.0);
	complex double epsL=2.1229-j*0.0;
	complex double epsH=4.5967-j*0.0;
	// General Parameters
	myConfig->lambda0 = 633.0*nm;
	myConfig->freq = c0/myConfig->lambda0;
	myConfig->k0 = 2.0*pi/myConfig->lambda0;
	myConfig->Gamma_u = 0.0;
	myConfig->Gamma_d = 0.0;
	int n=0;
	// Layer 1
	myConfig->Layers[n].zu = +2500.0*nm;
	myConfig->Layers[n].zd = +1523.0*nm;
	myConfig->Layers[n].eps = 2.3013-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 2
	myConfig->Layers[n].zu = +1523.0*nm;
	myConfig->Layers[n].zd = +1445.0*nm;
	myConfig->Layers[n].eps = epsH;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 3
	myConfig->Layers[n].zu = +1445.0*nm;
	myConfig->Layers[n].zd = +1329.0*nm;
	myConfig->Layers[n].eps = epsL;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 4
	myConfig->Layers[n].zu = +1329.0*nm;
	myConfig->Layers[n].zd = +1241.0*nm;
	myConfig->Layers[n].eps = epsH;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 5
	myConfig->Layers[n].zu = +1241.0*nm;
	myConfig->Layers[n].zd = +1115.0*nm;
	myConfig->Layers[n].eps = epsL;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 6
	myConfig->Layers[n].zu = +1115.0*nm;
	myConfig->Layers[n].zd = +1037.0*nm;
	myConfig->Layers[n].eps = epsH;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 7
	myConfig->Layers[n].zu = +1037.0*nm;
	myConfig->Layers[n].zd = +911.0*nm;
	myConfig->Layers[n].eps = epsL;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 8
	myConfig->Layers[n].zu = +911.0*nm;
	myConfig->Layers[n].zd = +833.0*nm;
	myConfig->Layers[n].eps = epsH;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 9
	myConfig->Layers[n].zu = +833.0*nm;
	myConfig->Layers[n].zd = +707.0*nm;
	myConfig->Layers[n].eps = epsL;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 10
	myConfig->Layers[n].zu = +707.0*nm;
	myConfig->Layers[n].zd = +629.0*nm;
	myConfig->Layers[n].eps = epsH;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 11
	myConfig->Layers[n].zu = +629.0*nm;
	myConfig->Layers[n].zd = +503.0*nm;
	myConfig->Layers[n].eps = epsL;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 12
	myConfig->Layers[n].zu = +503.0*nm;
	myConfig->Layers[n].zd = +425.0*nm;
	myConfig->Layers[n].eps = epsH;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 13
	myConfig->Layers[n].zu = +425.0*nm;
	myConfig->Layers[n].zd = +299.0*nm;
	myConfig->Layers[n].eps = epsL;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 14
	myConfig->Layers[n].zu = +299.0*nm;
	myConfig->Layers[n].zd = +221.0*nm;
	myConfig->Layers[n].eps = epsH;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 15
	myConfig->Layers[n].zu = +221.0*nm;
	myConfig->Layers[n].zd = +69.0*nm;
	myConfig->Layers[n].eps = epsL;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 16
	myConfig->Layers[n].zu = +69.0*nm;
	myConfig->Layers[n].zd = +42.0*nm;
	myConfig->Layers[n].eps = epsL;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 17
	myConfig->Layers[n].zu = +42.0*nm;
	myConfig->Layers[n].zd = +0.0*nm;
	myConfig->Layers[n].eps = -18.3511-j*0.4331;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 18
	myConfig->Layers[n].zu = +0.0*nm;
	myConfig->Layers[n].zd = -27.0*nm;
	myConfig->Layers[n].eps = epsL;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
	// Layer 19
	myConfig->Layers[n].zu = -27.0*nm;
	myConfig->Layers[n].zd = -500.0*nm;
	myConfig->Layers[n].eps = 1.0-j*0.0;
	myConfig->Layers[n].mu = 1.0-j*0.0;
	myConfig->Layers[n].sigma_s = 0.0-j*0.0;
	myConfig->Layers[n].k = myConfig->k0*csqrt(myConfig->Layers[n].mu*myConfig->Layers[n].eps);
	myConfig->Layers[n].eta = eta0*csqrt(myConfig->Layers[n].mu/myConfig->Layers[n].eps);
	myConfig->Layers[n].d = myConfig->Layers[n].zu-myConfig->Layers[n].zd;
	assert(myConfig->Layers[n].d>=0.0);
	n++;
}
