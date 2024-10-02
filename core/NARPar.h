#include <vector>
#include <memory>
#include "odePar.h"
#pragma once
class NARPar : public ODEPar
{
private:
    double _AUC = 2.5353073611315153; // default value when all parameters are 1
public:
    NARPar(double AUC, std::vector<double> pars);
    NARPar();

    const double& AUC() const { return _AUC; }
    const double& aZ() const { return _pars[0]; }
    const double& bZ() const { return _pars[1]; }
    const double& KZ() const { return _pars[2]; }
    const double& KXZ() const { return _pars[3]; }

    double& AUC() { return _AUC; }
    double& aZ() { return _pars[0]; }
    double& bZ() { return _pars[1]; }
    double& KZ() { return _pars[2]; }
    double& KXZ() { return _pars[3]; }

};
