#include "odePar.h"

ODEPar::ODEPar(/* args */)
{
}

ODEPar::~ODEPar()
{
}

void ODEPar::setParValue(size_t i, double val)
{
    switch (i)
    {
    case 0:
        _aZ = val;
        break;
    case 1:
        _KZ = val;
        break;
    case 2:
        _bZ = val;
        break;
    case 3:
        _KXZ = val;
        break;
    default:
        break;
    }

    return;
}

std::vector<double> ODEPar::getPars()
{
    return std::vector<double>({ _AUC, _aZ, _bZ, _KZ, _KXZ });
} 

// Get an ODEPar from a vector of ODEPars
double ODEPar::getODEValFromVector(const ODEPar& target, const std::vector<std::unique_ptr<ODEPar>>& vec)
{
    for (const auto& ODE : vec)
    {
        if (target == *ODE.get())
            return ODE->AUC();
    }
    return 0;
} 
