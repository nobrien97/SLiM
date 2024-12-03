#include "odePar.h"
#include "NARPar.h"
#include "PARPar.h"
#include "FFLC1Par.h"
#include "FFLI1Par.h"
#include "FFBHPar.h"


std::unique_ptr<ODEPar> ODEPar::MakeODEPtr(motif_enum motifType)
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

std::unique_ptr<ODEPar> ODEPar::MakeODEPtr(motif_enum motifType, const ODEPar &initialODEPar)
{
    switch (motifType)
    {
        case NAR:
            return std::make_unique<NARPar>(initialODEPar._solutionTraits, initialODEPar._pars);
            break;
        case PAR:
            return std::make_unique<PARPar>(initialODEPar._solutionTraits, initialODEPar._pars);
            break;
        case FFLC1:
            return std::make_unique<FFLC1Par>(initialODEPar._solutionTraits, initialODEPar._pars);
            break;
        case FFLI1:
            return std::make_unique<FFLI1Par>(initialODEPar._solutionTraits, initialODEPar._pars);
            break;
        case FFBH:
            return std::make_unique<FFBHPar>(initialODEPar._solutionTraits, initialODEPar._pars);
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

ODEPar::ODEPar(int numPar, int numTrait, std::vector<double> traits, std::vector<double> pars) : numPars(numPar), numTraits(numTrait)
{
    _pars.resize(pars.size());
    for (size_t i = 0; i < numPars; ++i)
    {
        _pars[i] = pars[i];
    }

    _solutionTraits.resize(traits.size());
    for (size_t j = 0; j < numTraits; ++j)
    {
        _solutionTraits[j] = traits[j]; 
    }

}

// Get an ODEPar from a vector of ODEPars
std::vector<double> ODEPar::getODEValFromVector(const ODEPar& target, const std::vector<std::unique_ptr<ODEPar>>& vec, bool incrementCount)
{
    for (auto& ODE : vec)
    {
        ODEPar source = *(ODE.get());
        if (target == source)
        {
            if (incrementCount) 
                ++*ODE;

            return source._solutionTraits;
        }
    }
    return {0};
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

// Calculates the response time and steady state concentration for a given ODE solution
std::vector<double> ODEPar::CalcSteadyState(const asc::Recorder &solution, const double& startTime, const int &solutionIndex, const bool &reverseOrder)
{
    // First get the steady state and halfway point
    std::vector<double> result{0.0, 0.0};
    double half = 0.0;
    static float epsilon = 0.001f;
    int steadyCount = 0;
    static int maxSteadyCount = 4;

    // Find start index in solution history
    int startIndex = (int)startTime * 10; // TODO: HACK: multiply by 10 because the sampling rate is 0.1

    if (startIndex >= solution.history.size()) {
        startIndex = solution.history.size() - 1;
    }

    // If we're looking forwards
    if (!reverseOrder)
    {
        for (int i = startIndex + 1; i < solution.history.size() - 1; ++i)
        {
            // concentrations for current and previous time points
            double c1 = solution.history[i-1][solutionIndex];
            double c2 = solution.history[i][solutionIndex];

            if (std::abs(c2 - c1) < epsilon) {
                steadyCount++;
                result[0] = c2; // steady state value

                if (steadyCount >= maxSteadyCount) {
                    break;
                }
            }
        }

        // If there's no steady state, return early
        if (result[0] <= epsilon) {
            return result;
        }

        // Find the response time
        half = result[0] * 0.5;

        // Figure out where the halfway point is
        for (int i = 1; i < solution.history.size() - 1; ++i)
        {
            double t1 = solution.history[i-1][0];
            double t2 = solution.history[i][0];
            double c1 = solution.history[i-1][solutionIndex];
            double c2 = solution.history[i][solutionIndex];

            if ((c1 < half && c2 >= half) || (c1 > half && c2 <= half)) {
                result[1] = Interpolate(t1, c1, t2, c2, half); // response time
                break;
            }
        }

        return result;
    }


    // if we're working backwards, the for loops are different
    for (int i = startIndex; i > 0; --i)
    {
        // concentrations for current and previous time points
        double c1 = solution.history[i][solutionIndex];
        double c2 = solution.history[i-1][solutionIndex];

        // we want to check when the steady state stops being stable, this is the 
        // point where the steady state starts (assuming steady state is reached at end of evaluation)
        if (std::abs(c2 - c1) > epsilon) {
            result[0] = c1;
            break;
        }
    }
    // If there's no steady state, return early
    if (result[0] <= epsilon) {
        return result;
    }

    // Find the response time
    half = result[0] * 0.5;
    // Figure out where the halfway point is
    for (int i = solution.history.size() - 1; i > 1; --i)
    {
        double t1 = solution.history[i][0];
        double t2 = solution.history[i-1][0];
        double c1 = solution.history[i][solutionIndex];
        double c2 = solution.history[i-1][solutionIndex];
        if ((c1 < half && c2 >= half) || (c1 > half && c2 <= half)) {
            result[1] = Interpolate(t1, c1, t2, c2, half);
            break;
        }
    }
    return result;
}

// Finds the delay from a certain point before n ODE starts changing
double ODEPar::CalcDelayTime(const asc::Recorder &solution, const double &startTime, const int &solutionIndex)
{
    double result = 0.0;
    static float epsilon = 0.001f;

    int startIndex = (int)startTime * 10; // TODO: HACK: multiply by 10 because the sampling rate is 0.1

    if (startIndex >= solution.history.size()) {
        startIndex = solution.history.size() - 1;
    }

    for (int i = startIndex + 1; i < solution.history.size(); ++i)
    {
        double t = solution.history[i][0];
        double c1 = solution.history[i-1][solutionIndex];
        double c2 = solution.history[i][solutionIndex];

        if (std::abs(c2 - c1) < epsilon) {
            result = t;
            break;
        }
    }

    return result;
}

// Maximum expression and the time at which it is reached
std::vector<double> ODEPar::CalcMaxExpression(const asc::Recorder &solution, const int &solutionIndex)
{
    std::vector<double> result {0.0, 0.0};
    double curMax = 0.0;
    double curTime = 0.0;
    for (int i = 0; i < solution.history.size(); ++i)
    {
        double curVal = solution.history[i][solutionIndex];
        if (curVal > curMax) {
            curMax = curVal;
            curTime = solution.history[i][0];
        }
    }

    result[0] = curMax;
    result[1] = curTime;
    return result;
}

// Amount of time spent above a threshold
double ODEPar::CalcTimeAboveThreshold(const asc::Recorder &solution, const double &threshold, const int &solutionIndex)
{
    static double DELTA = 0.1; // delta between two measurements
    double timeAboveThreshold = 0.0;

    for (int i = 0; i < solution.history.size(); ++i)
    {
        double curVal = solution.history[i][solutionIndex];

        if (curVal >= threshold) {
            timeAboveThreshold += DELTA;
        }
    }

    return timeAboveThreshold;
}