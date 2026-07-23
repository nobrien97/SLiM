#include <vector>
#include <memory>
#include "odePar.h"
#pragma once
class ComplexityPar : public ODEPar
{
private:
    double _AUC = 0.0;
public:
    ComplexityPar(std::vector<double> traits, std::vector<double> pars);
    ComplexityPar();

    std::vector<double> SolveODE() override;

    // Hardcoded number of parameters/traits
    const static int numTraits = 5;

    // Number of parameters is alpha and beta per trait, h, then a K and p per combination numTraits^2
    const static int numPars = numTraits + numTraits + 1 + 2 * (numTraits * numTraits);

    // Molecular components
    const double& AUC() const { return _AUC; }

    double& AUC() { return _AUC; }

    // Parameters
    /*
    In order of Gene 0 - nTraits
    h
    beta0, b1, b2, ..., bn
    alpha0, a1, a2, ..., an
    K00, K01, K02, ..., K0n, K10, ..., Knn
    p00, p01, p02, ..., p0n, p10, ..., pnn
    */


};
