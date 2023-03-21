//
#include "QuadL.h"

extern void cgqf_f77_(int *rule, int *order, double *x, double *w);

void QuadL(int N, double *x, double *w){
	int rule=1;
	cgqf_f77_(&rule, &N, x, w);
}

void printQuadLRule(int Nmax){
	for (int N=2; N<=Nmax; N*=2){
		double x[N], w[N];
		QuadL(N, x, w);
		printf("double x%d[%d] = {\n", N, N);
		for (int i=0; i<N; i++){
			printf("%21.14E", x[i]);
			if (i<N){
				printf(",\n");
			}else{
				printf("\n");
			}
		}
		printf("};\n\n");
		printf("double w%d[%d] = {\n", N, N);
		for (int i=0; i<N; i++){
			printf("%21.14E", w[i]);
			if (i<N){
				printf(",\n");
			}else{
				printf("\n");
			}
		}
		printf("};\n\n");
	}
}

double x32[32] = {
-9.97263861849481E-01,
-9.85611511545268E-01,
-9.64762255587506E-01,
-9.34906075937740E-01,
-8.96321155766053E-01,
-8.49367613732570E-01,
-7.94483795967942E-01,
-7.32182118740290E-01,
-6.63044266930215E-01,
-5.87715757240762E-01,
-5.06899908932229E-01,
-4.21351276130636E-01,
-3.31868602282128E-01,
-2.39287362252137E-01,
-1.44471961582796E-01,
-4.83076656877384E-02,
 4.83076656877379E-02,
 1.44471961582797E-01,
 2.39287362252137E-01,
 3.31868602282128E-01,
 4.21351276130635E-01,
 5.06899908932229E-01,
 5.87715757240762E-01,
 6.63044266930215E-01,
 7.32182118740290E-01,
 7.94483795967943E-01,
 8.49367613732570E-01,
 8.96321155766052E-01,
 9.34906075937740E-01,
 9.64762255587507E-01,
 9.85611511545268E-01,
 9.97263861849482E-01,
};

double w32[32] = {
 7.01861000947075E-03,
 1.62743947309056E-02,
 2.53920653092616E-02,
 3.42738629130212E-02,
 4.28358980222266E-02,
 5.09980592623755E-02,
 5.86840934785365E-02,
 6.58222227763618E-02,
 7.23457941088485E-02,
 7.81938957870707E-02,
 8.33119242269467E-02,
 8.76520930044036E-02,
 9.11738786957633E-02,
 9.38443990808042E-02,
 9.56387200792752E-02,
 9.65400885147271E-02,
 9.65400885147275E-02,
 9.56387200792749E-02,
 9.38443990808046E-02,
 9.11738786957628E-02,
 8.76520930044037E-02,
 8.33119242269474E-02,
 7.81938957870695E-02,
 7.23457941088490E-02,
 6.58222227763622E-02,
 5.86840934785361E-02,
 5.09980592623757E-02,
 4.28358980222272E-02,
 3.42738629130216E-02,
 2.53920653092625E-02,
 1.62743947309049E-02,
 7.01861000947000E-03,
};

complex double QuadL1D_32(complex double func(complex double, void*),
        void *args, double a, double b){
    double hm=(b-a)/2.0;
    double hp=(b+a)/2.0;
    complex double I=0.0;
    for (int n=0; n<32; n++){
        I+=hm*w32[n]*func(hm*x32[n]+hp, args);
    }
    return I;
}

complex double QuadLAdaptive1D_32(complex double func(complex double, void*),
        void *args, double a, double b, double tol,
        int kmax, int k, complex double Ip, int *flag){
    assert(kmax>0);
    if (k==0){
        Ip = QuadL1D_32(func, args, a, b);
    }
    k++;
    double m=0.5*(a+b);
    complex double I1=QuadL1D_32(func, args, a, m);
    complex double I2=QuadL1D_32(func, args, m, b);
    complex double In=I1+I2;
    double err=cabs(In-Ip);
    if ((err>(tol*cabs(In)))&&(k<=kmax)){
        I1 = QuadLAdaptive1D_32(func, args, 
            a, m, tol, kmax, k, I1, flag);
        I2 = QuadLAdaptive1D_32(func, args, 
            m, b, tol, kmax, k, I2, flag);
        In = I1+I2;
    }
    if (k>kmax){
        (*flag) = 1;
        printf("Warning: No Convergence!\n");
    }else{
        (*flag) = 0;
    }
    return In;
}