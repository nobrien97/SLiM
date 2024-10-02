#include "NARPar.h"

NARPar::NARPar(double AUC, std::vector<double> pars) : ODEPar(4, pars)
{
    _AUC = AUC;
    _pars.resize(numPars, 1.0);
}

NARPar::NARPar() : ODEPar(4) {}