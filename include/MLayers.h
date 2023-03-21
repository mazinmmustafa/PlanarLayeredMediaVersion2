#ifndef MLAYERS_H
#define MLAYERS_H

// Definitions
#include "myLib.h"

typedef struct Layer Layer;
struct Layer{
    double zu, zd, d;
    complex double eps, mu, k, eta, sigma_s;
};

typedef struct MLayers MLayers;
struct MLayers{
    int N;
    double lambda0, freq, k0;
    complex double Gamma_u, Gamma_d;
    Layer *Layers;
};
#define DefaultMLayers {0, 0, 0, 0, 0, 0, NULL}

// Functions
void logConfig(MLayers *myConfig);
void ConfigPaulus(MLayers *myConfig);
void ConfigGoldKretschmann(double lambda0, MLayers *myConfig);
void ConfigOttoGraphene(double freq, MLayers *myConfig);
void ConfigChew(MLayers *myConfig);
void ConfigGoldKretschmannNadson(double lambda0, MLayers *myConfig);
void ConfigPlasmonicWG(double lambda0, MLayers *myConfig);
void ConfigFEKOPlaneWave(MLayers *myConfig);
void ConfigBragg(MLayers *myConfig);

#endif