#ifndef FAFITISCONCRETE_H
#define FAFITISCONCRETE_H
#include "material.h"
#include<cmath>

class FafitisConcrete : public Material
{
    //type definition which is pointer to the method of this object
    //that contains formula to calculate sigm for given function part
    typedef double (FafitisConcrete::*SigmaEpsilon)(double);

    //FUNCTION PARTS
    const int noOfFunctionParts = 2; //parabolic-linear

    //STRAIN THRESHOLD VALUES
    double epsilon_tu = 0.0;
    double epsilon_co = -0.002;
    double epsilon_cu = -0.0035;

    //Strength values
    //double fc = 34.47; //[MPa] = -5Ksi
    //double fc = 5; //[Ksi]
    double fc = 9; //[Ksi]
    double alfa = 0.85; //strength coefficient!

public:
    FafitisConcrete() : Material(FAFITIS_CONCRETE) { }

    int getModuleOfElasticity_E(void) {
        //wc - density of concrete
        //double wc = 2400; // [kg/m^3]
        //double Ec = std::pow(wc, 1.5)*0.043*std::sqrt(std::abs(fc)); //[MPa]
        //return Ec/1000; [GPa]
        double wc = 150.0; // [lb/ft^3]
        //in below formula we need to convert fc[ksi] to f'c [psi] = fc*1000
        double Ec = std::pow(wc, 1.5)*33*std::sqrt(std::abs(alfa*fc*1000)); //[psi]

        //return 2125;
        return 2900;
        //in below we need to conver back from [psi] -> [ksi] where 1ksi = 1000psi
        return (int) (Ec/1000); //[ksi]
    }

    double tensileStrainLimit_eps_tu(void)
    {
        return epsilon_tu;
    }

    double compressiveStrainLimit_eps_cu(void)
    {
        return epsilon_cu;
    }

    double pureCompressionStrainLimit_eps_co(void)
    {
        return epsilon_co;
    }

    double sigmaEpsilon(double epsilon) {
        if(epsilon <= 0 && epsilon >= epsilon_co)
        {
            return alfa*fc*(1000*epsilon + 250000*std::pow(epsilon, 2.0)); //[MPa]
        } else if(epsilon < epsilon_co && epsilon >= epsilon_cu)
        {
            return -alfa*fc; //[MPa]
        }
        //concrete doesn't bear tensile stresses
        return 0; //[MPa]
    }

    int numberOfFunctionParts(void)
    {
        return noOfFunctionParts;
    }

    int functionPartNumberFor(double epsilon)
    {
        if(epsilon < epsilon_cu)
        {
            throw StrainOutOfBoundsException("Przekroczone odkształcenia dla Fafitis_Concrete");

        } else if(epsilon < epsilon_co)
        {
            return 1; //linear function part f1
        } else if(epsilon <= 0)
        {
            return 0; //parabolic funcion part f0
        } else {
            //tension
            throw StrainOutOfBoundsException("Przekroczone odkształcenia dla Fafitis_Concrete");
        }
    }

    double functionPartStrainLimitA(int number)
    {
        switch(number) {
            case 0:
                //parabolic sigma-epsilon, function part f0
                return 0;
            case 1:
                //linear sigma-epsilon, function part f1
                return epsilon_co;
        }
    }

    double functionPartStrainLimitB(int number)
    {
        switch(number) {
            case 0:
                //parabolic sigma-epsilon, function part f0
                return epsilon_co;
            case 1:
                //linear sigma-epsilon, function part f1
                return epsilon_cu;
        }
    }

    double functionPartStrainLimitA(double epsilon)
    {
        if(epsilon < epsilon_cu || epsilon > 0)
        {
            throw StrainOutOfBoundsException("Przekroczone odkształcenia dla Fafitis_Concrete");
        }
        return functionPartStrainLimitA(functionPartNumberFor(epsilon));
    }

    double functionPartStrainLimitB(double epsilon)
    {
        if(epsilon < epsilon_cu || epsilon > 0)
        {
            throw StrainOutOfBoundsException("Przekroczone odkształcenia dla Fafitis_Concrete");
        }
        return functionPartStrainLimitB(
                    functionPartNumberFor(epsilon));
    }

    int functionPartDegree(int number)
    {
        switch(number) {
            case 0:
                return 2;
            case 1:
                return 0;
        }
    }

    int functionPartDegree(double epsilon)
    {
        if(epsilon < epsilon_cu || epsilon > 0)
            throw StrainOutOfBoundsException("Przekroczone odkształcenia dla Fafities_Concrete");

        return functionPartDegree(functionPartNumberFor(epsilon));
    }

    Qt::GlobalColor materialColor(void) {
        return Qt::lightGray;
    }
};

#endif // FAFITISCONCRETE_H
