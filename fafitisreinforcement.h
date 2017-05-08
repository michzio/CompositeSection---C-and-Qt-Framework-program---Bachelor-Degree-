#ifndef FAFITISREINFORCEMENT_H
#define FAFITISREINFORCEMENT_H
#include "material.h"

class FafitisReinforcement : public Material
{
    //Fafitis reinforcement is elasto-plastic
    //with yielding stress fy = 60ksi
    //and elastic modulus of elasticity Es = 30,000ksi

    //type definition which is pointer to the method of this object
    //that contains forumula to calculate sigma for given function part of sigma-epsilon relation
    typedef double (FafitisReinforcement::*SigmaEpsilon)(double);

    //number of function parts
    const int noOfFunctionParts = 3; //bilinear, elasto-plastic

    //modulus of elasticity Es
    static const int Es = 30000; //[ksi]

    //characteristic value f_yk
    static const int f_yk = 60; //[ksi]

    //ultimate strain e_uk
    double epsilon_uk = 0.025; //[-]
    double epsilon_yk; //[-] f_yk/Es

public:
    FafitisReinforcement() : Material(FAFITIS_REINFORCEMENT) {}

    double tensileStrainLimit_eps_tu(void)
    {
        return epsilon_uk; //[-]
    }

    double compressiveStrainLimit_eps_cu(void)
    {
        return -epsilon_uk; //[-]
    }

    double pureCompressionStrainLimit_eps_co()
    {
        return compressiveStrainLimit_eps_cu(); //[-]
    }

    int getModuleOfElasticity_E(void)
    {
        return Es; //[ksi]
    }

    double sigmaEpsilon(double epsilon) {

        double epsilon_yk = (double) f_yk/Es;

        if(epsilon > epsilon_uk && epsilon < -epsilon_uk)
        {
            throw StrainOutOfBoundsException("Przekroczenie odkształceń granicznych w Fafitis_Reinforcement");
        } else if(epsilon <= epsilon_uk && epsilon > epsilon_yk)
        {
            return f_yk; //[ksi]
        } else if(epsilon <= epsilon_yk && epsilon >= -epsilon_yk)
        {
            return epsilon*Es; //[ksi]
        } else if(epsilon < -epsilon_yk && epsilon >= -epsilon_uk) {
            return -f_yk; //[ksi]
        }
    }

    int numberOfFunctionParts(void) {
        return noOfFunctionParts;
    }

    int functionPartNumberFor(double epsilon)
    {
        double epsilon_yk = (double) f_yk/Es;

        if(epsilon > epsilon_uk && epsilon < -epsilon_uk)
        {
            throw StrainOutOfBoundsException("Przekroczenie odkształceń granicznych w Fafitis_Reinforcement");
        } else if(epsilon <= epsilon_uk && epsilon > epsilon_yk)
        {
           return 0;

        } else if(epsilon <= epsilon_yk && epsilon >= -epsilon_yk)
        {
           return 1;

        } else if(epsilon < -epsilon_yk && epsilon >= -epsilon_uk)
        {
           return 2;
        }
    }

    double functionPartStrainLimitA(int number)
    {
        double epsilon_yk = (double) f_yk/Es;

        switch(number) {
            case 0:
                return epsilon_uk;
            case 1:
                return epsilon_yk;
            case 2:
                return -epsilon_yk;
        }
    }

    double functionPartStrainLimitB(int number)
    {
        double epsilon_yk = (double) f_yk/Es;

        switch(number) {
            case 0:
                return epsilon_yk;
            case 1:
                return -epsilon_yk;
            case 2:
                return -epsilon_uk;
        }
    }

    double functionPartStrainLimitA(double epsilon)
    {
        if(epsilon < -epsilon_uk || epsilon > epsilon_uk)
            throw StrainOutOfBoundsException("Przekroczenie odkształceń w Fafitis_Reinforcement");

        return functionPartStrainLimitA(functionPartNumberFor(epsilon));
    }

    double functionPartStrainLimitB(double epsilon)
    {
        if(epsilon < -epsilon_uk || epsilon > epsilon_uk)
            throw StrainOutOfBoundsException("Przekroczenie odkształceń w Fafitis_Reinforcement");

        return functionPartStrainLimitB(functionPartNumberFor(epsilon));
    }

    int functionPartDegree(int number)
    {
        switch(number)
        {
            case 0:
                return 1;
            case 1:
                return 1;
            case 2:
                return 1;
        }
    }

    int functionPartDegree(double epsilon)
    {
        if(epsilon < -epsilon_uk || epsilon > epsilon_uk)
            throw StrainOutOfBoundsException("Przekroczenie odkształceń w Fafitis_Reinforcement");

        return functionPartDegree(functionPartNumberFor(epsilon));
    }

    Qt::GlobalColor materialColor(void) {
        return Qt::darkGray;
    }

};

#endif // FAFITISREINFORCEMENT_H
