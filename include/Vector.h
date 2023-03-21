#ifndef VECTOR_H
#define VECTOR_H

// Includes
#include "myLib.h"

// Definitions
typedef struct Vector Vector;
struct Vector{
    double x, y, z;
};

typedef struct ComplexVector ComplexVector;
struct ComplexVector{
    complex double x, y, z;
};

// Functions
double magVector(Vector A);
double dotVector(Vector A, Vector B);
Vector crossVector(Vector A, Vector B);
Vector unitVector(Vector A);
Vector addVector(Vector A, Vector B);
Vector subVector(Vector A, Vector B);
Vector scaleVector(Vector A, double a);
int isEqualVector(Vector A, Vector B);
complex double dotComplexVector(ComplexVector A, ComplexVector B);
complex double dotComplexVector2Vector(ComplexVector A, Vector B);
ComplexVector addComplexVector(ComplexVector A, ComplexVector B);
ComplexVector subComplexVector(ComplexVector A, ComplexVector B);
ComplexVector scaleComplexVector(ComplexVector A, complex double a);
ComplexVector scaleVector2ComplexVector(Vector A, complex double a);

#endif