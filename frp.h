#ifndef FRP_H
#define FRP_H

#include "material.h"

class FRP : public Material
{
    constexpr static const double EPSILON_ERROR_OFFSET = 0.00000001;
    constexpr static const double MODULUS_OF_ELASTICITY = 240;

    constexpr static const double EPSILON_TU = 0.006; //[-]
    constexpr static const double EPSILON_CU = -0.006; //[-]

    double e_tu;
    double e_cu;
    double E;

public:
    FRP() : Material(FIBER_REINFORCED_PLASTIC) {  e_tu = EPSILON_TU; e_cu = EPSILON_CU;  E = MODULUS_OF_ELASTICITY;}

    int getModuleOfElasticity_E(void)
    {
        return E;
    }

    double tensileStrainLimit_eps_tu(void)
    {
        return e_tu;
    }

    double compressiveStrainLimit_eps_cu(void)
    {
        return e_cu;
    }

    double pureCompressionStrainLimit_eps_co(void)
    {
        return e_cu;
    }

    Qt::GlobalColor materialColor(void)
    {
        return Qt::red;
    }

    double sigmaEpsilon(double epsilon) {
        if(epsilon > e_tu + EPSILON_ERROR_OFFSET || epsilon < e_cu - EPSILON_ERROR_OFFSET) {
            fprintf(stderr, "Odkształcenia epsilon = %g, e_tu = %g, e_cu = %g\n", epsilon, e_tu, e_cu);
            throw StrainOutOfBoundsException("Przekroczenie odkształceń w taśmie FRP!");
            return 0;
        }
        if(epsilon < 0) return 0;

        return epsilon*E*1000;

    }

    int numberOfFunctionParts(void) {
        return 1;
    }

    int functionPartNumberFor(double epsilon)
    {
        if(epsilon > e_tu + EPSILON_ERROR_OFFSET || epsilon < e_cu - EPSILON_ERROR_OFFSET) {
            throw StrainOutOfBoundsException("Przekroczenie odkształceń w taśmie FRP!");
        }
        return 0;
    }

    double functionPartStrainLimitA(int number) {

        switch(number) {
            case 0:
                return e_tu;
        }
    }

    double functionPartStrainLimitB(int number)
    {
        switch(number) {
            case 0:
                return e_cu;
        }
    }

    double functionPartStrainLimitA(double epsilon)
    {
        if(epsilon > e_tu + EPSILON_ERROR_OFFSET || epsilon < e_cu - EPSILON_ERROR_OFFSET) {
            throw StrainOutOfBoundsException("Przekroczenie odkształceń w taśmie FRP!");
        }

         return functionPartStrainLimitA(functionPartNumberFor(epsilon));
    }

    double functionPartStrainLimitB(double epsilon)
    {
        if(epsilon > e_tu + EPSILON_ERROR_OFFSET || epsilon < e_cu - EPSILON_ERROR_OFFSET) {
            throw StrainOutOfBoundsException("Przekroczenie odkształceń w taśmie FRP!");
        }

         return functionPartStrainLimitB(functionPartNumberFor(epsilon));
    }

    int functionPartDegree(int number)
    {
        switch(number)
        {
            case 0:
                return 1;
        }
    }

    int functionPartDegree(double epsilon)
    {
        if(epsilon > e_tu + EPSILON_ERROR_OFFSET || epsilon < e_cu - EPSILON_ERROR_OFFSET) {
            throw StrainOutOfBoundsException("Przekroczenie odkształceń w taśmie FRP!");
        }

         return functionPartDegree(functionPartNumberFor(epsilon));
    }

};

#endif // FRP_H
