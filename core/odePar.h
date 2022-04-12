class ODEPar
{
private:
    double finalZ;
    double aZ;
	double bZ;
	double KZ;
	double KXZ;

public:
    ODEPar(/* args */);
    ~ODEPar();
    bool operator==(const ODEPar& rhs) const 
    {
        return 
            aZ == rhs.aZ &&
            bZ == rhs.bZ &&
            KZ == rhs.KZ &&
            KXZ == rhs.KXZ
        ;
    }
    // Set a given value
    void setParValue(size_t i, double val);

    // Get an ODEPar from a vector of ODEPars
    static double getODEValFromVector(const ODEPar& target, const std::vector<ODEPar>& vec)
    {
        for (auto ODE : vec)
        {
            if (target == ODE)
                return ODE.finalZ;
        }
    } 
};

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
        aZ = val;
        break;
    case 1:
        KZ = val;
        break;
    case 2:
        bZ = val;
        break;
    case 3:
        KXZ = val;
        break;
    default:
        break;
    }

    return;
}


