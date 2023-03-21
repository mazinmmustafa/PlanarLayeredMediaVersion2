#ifndef QUADL_H
#define QUADL_H

// Definitions 
#include "myLib.h"

typedef struct QuadResults QuadResults;
struct QuadResults{
	complex double I;
	int flag;
};

// Functions
complex double QuadLAdaptive1D_32(complex double func(complex double, void*),
        void *args, double a, double b, double tol,
        int kmax, int k, complex double Ip, int *flag);
#endif