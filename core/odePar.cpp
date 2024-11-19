#include "odePar.h"
#include "NARPar.h"
#include "PARPar.h"
#include "FFLC1Par.h"
#include "FFLI1Par.h"
#include "FFBHPar.h"


std::unique_ptr<ODEPar> ODEPar::MakeODEPtr(motif_enum motifType)
{
        {
        switch (motifType)
        {
            case NAR:
                return std::make_unique<NARPar>();
                break;
            case PAR:
                return std::make_unique<PARPar>();
                break;
            case FFLC1:
                return std::make_unique<FFLC1Par>();
                break;
            case FFLI1:
                return std::make_unique<FFLI1Par>();
                break;
            case FFBH:
                return std::make_unique<FFBHPar>();
                break;
            default:
                return nullptr;
                break;
        }
    }
}

std::unique_ptr<ODEPar> ODEPar::MakeODEPtr(motif_enum motifType, const ODEPar &initialODEPar)
{
    switch (motifType)
    {
        case NAR:
            return std::make_unique<NARPar>(initialODEPar._AUC, initialODEPar._pars);
            break;
        case PAR:
            return std::make_unique<PARPar>(initialODEPar._AUC, initialODEPar._pars);
            break;
        case FFLC1:
            return std::make_unique<FFLC1Par>(initialODEPar._AUC, initialODEPar._pars);
            break;
        case FFLI1:
            return std::make_unique<FFLI1Par>(initialODEPar._AUC, initialODEPar._pars);
            break;
        case FFBH:
            return std::make_unique<FFBHPar>(initialODEPar._AUC, initialODEPar._pars);
            break;
        default:
            return nullptr;
            break;
    }
}

// No-op default implementation: return an empty vector
std::vector<double> ODEPar::SolveODE()
{
    return std::vector<double>(0);
}

bool ODEPar::Compare(const ODEPar rhs)
{
    if (numPars != rhs.numPars)
    {
        return false;
    }
    int sum = 0;
    for (size_t i = 0; i < numPars; ++i)
    {
        sum += (_pars[i] == rhs._pars[i] ? 1 : 0);
    }
    return sum == numPars;
}

ODEPar::ODEPar(int pars) : numPars(pars), _pars(pars, 1.0) {}

ODEPar::ODEPar(int numPar, double AUC, std::vector<double> pars) : numPars(numPar)
{
    _pars.resize(pars.size());
    for (size_t i = 0; i < numPars; ++i)
    {
        _pars[i] = pars[i];
    }
    _AUC = AUC;

}

// Get an ODEPar from a vector of ODEPars
double ODEPar::getODEValFromVector(const ODEPar& target, const std::vector<std::unique_ptr<ODEPar>>& vec, bool incrementCount)
{
    for (auto& ODE : vec)
    {
        ODEPar source = *(ODE.get());
        if (target == source)
        {
            if (incrementCount) 
                ++*ODE;

            return source._AUC;
        }
    }
    return 0;
} 

std::vector<double> ODEPar::getPars(bool returnAUC)
{
    int n = numPars;
    int startIndex = 0;

    // If we're returning the AUC, make the vector bigger
    if (returnAUC)
    {
        startIndex = 1;
    }
    std::vector<double> result(n + startIndex);

    for (int i = 0; i < n; ++i)
    {
        result[i + startIndex] = _pars[i];
    }

    // First entry is AUC
    if (returnAUC)
    {
        result[0] = _AUC;
    }

    return result;
}

double ODEPar::getParValue(int i)
{
    return _pars[i];
}

void ODEPar::setParValue(int i, double val)
{
    _pars[i] = val;
    return;
}

void ODEPar::setParValue(std::vector<double> vals, bool firstElementIsAUC)
{
    int startIndex = 0;

    if (firstElementIsAUC)
    {
        startIndex++;
        setAUC(vals[0]);
    }    

    for (int i = startIndex; i < numPars; ++i)
    {
        setParValue(i, vals[i]);
    }
    return;
}

// Calculates area under the curve via the trapezoid method
#pragma omp declare simd
double ODEPar::AUC(const double &h, const double &a, const double &b)
{
	return ((a + b) * 0.5) * h;
}
