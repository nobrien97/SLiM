#include "FFLI1Par.h"

FFLI1Par::FFLI1Par(std::vector<double> traits, std::vector<double> pars) : ODEPar(numPars, numTraits, traits, pars)
{
    _pars.resize(numPars, 1.0);
	setParValue(pars, false);
	_solutionTraits.resize(numTraits);
	SetTraits(traits);
    //_pars[6] = 0.01; // constitutive promoter
}

FFLI1Par::FFLI1Par() : ODEPar(numPars) 
{
	_pars.resize(numPars, 1.0);
    //_pars[6] = 0.01; // constitutive promoter
	_solutionTraits.resize(numTraits);
}

std::vector<double> FFLI1Par::SolveODE()
{
    // Initialise output
    //std::vector<double> result(1, 0.0);

    // static components
	const double Xstart = 1.0; 
	const double Xstop = 6.0;

    // Starting X and Y
	int X = 0;

    // FFL I1 system
	auto FFLI1Derivative = [this, &Xstart, &Xstop, &X](const asc::state_t &curState, asc::state_t &nextState, double t)
	{
        // X <- XMult * (t > Xstart && t <= Xstop)
        // dY <- base * X + bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
        // dZ <- base * X + bZ * ((X * KY)^Hilln)/( (KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln) ) - aZ*Z
        
        double Xnew = X * XMult();
		double baseline = std::max(base() - 1.0, 0.0); // Adjust baseline so it is relative to the default value, 0 (instead of 1)

		nextState[0] = bY() * pow(Xnew, n()) / (pow(KY(), n()) + pow(Xnew, n())) - aY() * curState[0];
		nextState[1] = baseline + bZ() * pow(Xnew * KY(), n()) / ( (pow(KXZ(), n()) + pow(Xnew, n())) * (pow(KY(), n()) + pow(curState[0], n())) ) - aZ() * curState[1];    
    };

	// Set up the initial state
	asc::state_t state = { 0.0, 0.0 };
	double t = 0.0;
	double dt = 0.1;
	double t_end = 10.0;
	asc::RK4 integrator;
	asc::Recorder recorder;

	while (t < t_end)
	{
		// Add a small epsilon to get around t floating point inaccuracy
		X = ((t >= Xstart - 1e-5) && (t <= Xstop + 1e-5));
		recorder({t, (asc::value_t)X, state[0], state[1]});
		integrator(FFLI1Derivative, state, t, dt);
	}

	// Calculate traits
	std::vector<double> maxExp = ODEPar::CalcMaxExpression(recorder, Xstart, 3);
	SetMaxExpression(maxExp[0]);
	SetTimeToHalfMaxExpression(maxExp[1]);
	SetTimeAboveHalfMaxExpression(ODEPar::CalcTimeAboveThreshold(recorder, maxExp[0] * 0.5, 3)); // time above half max expression

	// Calculate AUC
	// double z = 0;
	// #pragma omp simd reduction(+:z)
	// for (uint i = 0; i < recorder.history.size()-1; ++i)
	// {
	// 	z += ODEPar::AUC(0.1, (double)recorder.history[i][3], (double)recorder.history[i + 1][3]);
	// }
	
	// Check that z is > 0, set AUC, return
	// z = (z >= 0) ? z : 0.0;
    // result[0] = z;
    // setAUC(z);
    return {TimeToHalfMaxExpression(), MaxExpression(), TimeAboveHalfMaxExpression()};
}
