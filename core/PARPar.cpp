#include "PARPar.h"

PARPar::PARPar(double AUC, std::vector<double> pars) : ODEPar(6, pars)
{
    _AUC = AUC;
    _pars.resize(numPars, 1.0);
    _pars[5] = 0.01;
}

PARPar::PARPar() : ODEPar(6) {}
