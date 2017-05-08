#ifndef HOLLOW_H
#define HOLLOW_H

#include "material.h"

//Exception thrown when user pass strain epsilon value out of bounds
struct HollowCalculationAttemptException : public std::exception
{
   std::string s;
   HollowCalculationAttemptException(std::string ss) : s(ss) {}
   ~HollowCalculationAttemptException() throw () {} // Updated
   const char* what() const throw() { return s.c_str(); }
};

class Hollow : public Material
{
public:
    Hollow() : Material(HOLLOW) { }

    int getModuleOfElasticity_E(void)
    {
        //hollow material i.e. opening hasn't modulus of elasticity
         return 0;
    }

    double tensileStrainLimit_eps_tu()
    {
        throw HollowCalculationAttemptException("Hollow::eps_tu()");
    }

    double compressiveStrainLimit_eps_cu()
    {
        //hollow doesn't limit strain profile
        throw HollowCalculationAttemptException("Hollow::eps_cu()");
    }

    double pureCompressionStrainLimit_eps_co()
    {
        //hollow doesn't limit strain profile
        throw HollowCalculationAttemptException("Hollow::eps_co()");
    }

    double sigmaEpsilon(double epsilon)
    {
        //hollow doesn't bear any stresses
        throw HollowCalculationAttemptException("Hollow::sigmaEpsilon()");
    }

    int numberOfFunctionParts(void)
    {
        throw HollowCalculationAttemptException("Hollow::numberOfFunctionParts()");
    }

    int functionPartNumberFor(double epsilon) {
       throw HollowCalculationAttemptException("Hollow::functionPartNumberFor()");
    }

    double functionPartStrainLimitA(int number) {
        throw HollowCalculationAttemptException("Hollow::functionPartStrainLimitA()");
    }

    double functionPartStrainLimitB(int number) {
        throw HollowCalculationAttemptException("Hollow::functionPartStrainLimitB()");
    }

    double functionPartStrainLimitA(double epsilon) {
        throw HollowCalculationAttemptException("Hollow::functionPartStrainLimitA()");
    }

    double functionPartStrainLimitB(double epsilon) {
        throw HollowCalculationAttemptException("Hollow::functionPartStrainLimitB()");
    }

    int functionPartDegree(int number) {
        throw HollowCalculationAttemptException("Hollow::functionPartDegree()");
    }

    int functionPartDegree(double epsilon) {
        throw HollowCalculationAttemptException("Hollow::functionPartDegree()");
    }

    Qt::GlobalColor materialColor(void) {
        return Qt::white;
    }
};


#endif // HOLLOW_H
