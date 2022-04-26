#include <vector>
#include <memory>
#pragma once
class ODEPar
{
private:
    double _AUC = 0.627;
    double _aZ = 1.0;
	double _bZ = 1.0;
	double _KZ = 1.0;
	double _KXZ = 1.0;

public:
    ODEPar(double AUC, double aZ, double bZ, double KZ, double KXZ) : 
            _AUC(AUC), _aZ(aZ), _bZ(bZ), _KZ(KZ), _KXZ(KXZ) {};
    ODEPar();
    ~ODEPar();

    bool operator==(const ODEPar rhs) const 
    {
        return 
            _aZ == rhs._aZ &&
            _bZ == rhs._bZ &&
            _KZ == rhs._KZ &&
            _KXZ == rhs._KXZ
        ;
    }
    // Set a given value
    void setParValue(size_t i, double val);

    // Get an ODEPar from a vector of ODEPars
    static double getODEValFromVector(const ODEPar& target, const std::vector<std::unique_ptr<ODEPar>>& vec);

    // Get all the values from an ODEPar
    std::vector<double> getPars();

    const double& AUC() const { return _AUC; }
    const double& aZ() const { return _aZ; }
    const double& bZ() const { return _bZ; }
    const double& KZ() const { return _KZ; }
    const double& KXZ() const { return _KXZ; }

    double& AUC() { return _AUC; }
    double& aZ() { return _aZ; }
    double& bZ() { return _bZ; }
    double& KZ() { return _KZ; }
    double& KXZ() { return _KXZ; }

    void setAUC(double val) { _AUC = val; }

};
