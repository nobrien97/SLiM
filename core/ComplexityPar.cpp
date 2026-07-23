#include "ComplexityPar.h"

ComplexityPar::ComplexityPar(std::vector<double> traits, std::vector<double> pars) : ODEPar(numPars, numTraits, traits, pars)
{
    _pars.resize(numPars, 1.0);
	setParValue(pars, false);
	_solutionTraits.resize(numTraits);
	SetTraits(traits);
    //_pars[4] = 0.0; // set baseline to 0 to start
}

ComplexityPar::ComplexityPar() : ODEPar(numPars) 
{
	_pars.resize(numPars, 1.0);
	_solutionTraits.resize(numTraits);
    //_pars[4] = 0.0; // set baseline to 0 to start
}

std::vector<double> ComplexityPar::SolveODE()
{
    // Initialise output
    //std::vector<double> result(1, 0.0);

    // static components
	const double Xstart = 1.0; 
	const double Xstop = 6.0;
	int X = 0;

	// Compound Hill function
	auto Hx = [](const int i, double* Xh, int* p, double* Kh)
	{
		double result = 1.0;
		for (int j = 0; i < numTraits; ++j)
		{
			result *= p[i] * (Xh[j] / (Xh[j] + Kh[i])) + (1 - p[i]) * (Kh[i] / (Xh[j] + Kh[i]));
		}
		return result;
	};

	// Lambdas for AUC and ODE system
	// Declare/define a lambda which defines the ODE system
	auto CompDerivative = [this, &Hx](const asc::state_t &val, asc::state_t &dxdt, double t)
	{
		// Setup parameters to avoid pow as much as possible
		static int numParsPerModel = 5; 
		double h = _pars[0];
		double beta[numTraits], alpha[numTraits], Xh[numTraits];
		double K[numTraits * numTraits];
		int p[numTraits * numTraits];

		// Fill beta and alpha
		int offset = 1;
		for (int i = 0; i < numTraits; ++i)
		{
			// Offset by one because h is the first entry
			alpha[i] = _pars[i + numTraits + offset];
			beta[i] = _pars[i + offset];

			// Pre-solve X^h
			Xh[i] = pow(val[i], h);
		}

		// Fill K and p
		offset += numTraits * 2;
		for (int i = 0; i < numTraits * numTraits; ++i)
		{
			// Offset, presolve K^h
			K[i] = pow(_pars[i + offset], h);
			p[i] = _pars[i + numTraits * numTraits + offset];
		}

		// Solve each equation
		// offset by number of iterations times the number of parameters
		for (int i = 0; i < numTraits; ++i)
		{
			dxdt[i] = beta[i] * Hx(i, Xh, p, K) - alpha[i] * val[i];
		}
	};

	// Set up the initial state
	asc::state_t state(numTraits, 0.1);
	double t = 0.0;
	double dt = 0.1;
	double t_end = 10.0;
	asc::RK4 integrator;
	asc::Recorder recorder;

	while (t < t_end)
	{
		// Add a small epsilon to get around t floating point inaccuracy
		//X = ((t >= Xstart - 1e-5) && (t <= Xstop + 1e-5));
		recorder({t, (asc::value_t)X, state[0]});
		integrator(CompDerivative, state, t, dt);
	}

	// Calculate traits
	std::vector<double> expression = ODEPar::CalcTotalExpression(recorder);

	for (int i = 0; i < numTraits; ++i)
	{
		_solutionTraits[i] = expression[i];
	}

    return {expression};
}
