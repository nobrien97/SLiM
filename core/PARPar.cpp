#include "PARPar.h"

PARPar::PARPar(std::vector<double> traits, std::vector<double> pars) : ODEPar(numPars, numTraits, traits, pars)
{
    _pars.resize(numPars, 1.0);
	setParValue(pars, false);
	_solutionTraits.resize(numTraits);
	SetTraits(traits);
    //_pars[4] = 0.01; // starting baseline expression
}

PARPar::PARPar() : ODEPar(numPars) 
{
	_pars.resize(numPars, 1.0);
    _pars[4] = 0.01; // set baseline to 0 to start
	_solutionTraits.resize(numTraits);
}

std::vector<double> PARPar::SolveODE()
{
    // Initialise output
    //std::vector<double> result(1, 0.0);

    // static components
	const double Xstart = 1.0; 
	const double Xstop = 6.0;
	int X = 0;

	// Lambdas for AUC and ODE system
	// Declare/define a lambda which defines the ODE system - this is going to be very ugly
	auto PARDerivative = [this, &Xstart, &Xstop, &X](const asc::state_t &val, asc::state_t &dxdt, double t)
	{
		double Xnew = X * XMult();
		double baseline = std::max(base() - 0.99, 0.0); // Adjust baseline so it is relative to the default value, 0.01 (instead of 1)
		// dZ <- base * X + bZ * (X^n/(KXZ^n + X^n)) * ((Z^n)/((KZ^n)+(Z^n))) - aZ*Z
		dxdt[0] = baseline + bZ() * pow(Xnew, n()) / (pow(KXZ(), n()) + pow(Xnew, n())) * (pow(val[0], n())/(pow(KZ(), n()) + pow(val[0], n()))) - aZ() * val[0];
	};

	// Set up the initial state
	asc::state_t state = { 0.0 };
	double t = 0.0;
	double dt = 0.1;
	double t_end = 10.0;
	asc::RK4 integrator;
	asc::Recorder recorder;

	while (t < t_end)
	{
		// Add a small epsilon to get around t floating point inaccuracy
		X = ((t >= Xstart - 1e-5) && (t <= Xstop + 1e-5));
		recorder({t, (asc::value_t)X, state[0]});
		integrator(PARDerivative, state, t, dt);
	}

	// Calculate traits
	std::vector<double> steadyState = ODEPar::CalcSteadyState(recorder, Xstart, Xstop, 2);
	SetSteadyState(steadyState[1]);
	SetResponseTime(steadyState[0]);	

	// NOTE: Cut response delay since it requires very specific molecular components
	// For the response delay, need to calculate the dt coming from the environmental input X
	//SetResponseDelay(ODEPar::CalcDelayTime(recorder, Xstart, Xstop, 2, aZ(), base()));

	// Calculate AUC
	// double z = 0;
	// #pragma omp simd reduction(+:z)
	// for (uint i = 0; i < recorder.history.size()-1; ++i)
	// {
	// 	z += ODEPar::AUC(0.1, (double)recorder.history[i][2], (double)recorder.history[i + 1][2]);
	// }
	
	// Check that z is > 0, set AUC, return
	// z = (z >= 0) ? z : 0.0;
    // result[0] = z;
    // setAUC(z);
    return {ResponseTime(), SteadyState()};
}
