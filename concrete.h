#ifndef CONCRETE_H
#define CONCRETE_H

#include "material.h"

//Strain values (parabolic-linear fig. 3.3 in PN-EN 1992-1-1) from Table 3.1 in PN-EN 1992-1-1
const double EPSILON_CO_LESS_THAN_50 = -0.002;
const double EPSILON_CU_LESS_THAN_50 = -0.0035;
                                        //-0.0038;
//EPSILON_TU IGNORED (ignoring concrete tensile strength)

//safety coefficients
const double ALFA_CC = 1.00;
const double GAMMA_C = //1.5;
                       1.0;

class Concrete : public Material
{
public:
    //type definition which is pointer to the method of this object
    //that contains forumla to calculate sigma for given function part of sigma-epsilon
    typedef double (Concrete::*SigmaEpsilon)(double);
private:

    //FUNCTION PARTS
    const int noOfFunctionParts = 2; //parabolic-linear

    //STRAIN THRESHOLD VALUES
    double epsilon_co; //strain corresponding to the compressive strength under pure compression
    double epsilon_cu; //strain corresponding to the ultimate compressive strength

    int f_ck; //characteristic strength of concrete [MPa]
    double n; //exponent in formula for parabolic sigma-epsilon relation
    int Ecm; //module of elasticity [GPa]
    int f_cm; //mean strength of concrete (f_ck + 8MPa) [MPa]
public:
    Concrete(int fck);

    void setCharacteristicStrength_fck(int);
    int getCharacteristicStrength_fck(void);
    double getDesignStrength_fcd(void);
    double getUltimateStrain_ecu(void);
    double getPureCompressionStrain_eco(void);
    double getExponentN(void);
    int getModuleOfElasticity_E(void);
    int getMeanStrength_fcm(void);
    //function returning tensile and compressive ultimate strain limits of current material
    double tensileStrainLimit_eps_tu(void);
    double compressiveStrainLimit_eps_cu(void);
    double pureCompressionStrainLimit_eps_co(void);

    //function takes epsilon (strain) and return sigma (stress)
    //depending on epsilon value uses parabolicSigmaEpsilon() or linearSigmaEpsilon()
    double sigmaEpsilon(double epsilon);

    //function parts
    int numberOfFunctionParts(void);
    int functionPartNumberFor(double epsilon);
    Concrete::SigmaEpsilon functionPart(int number);
    Concrete::SigmaEpsilon functionPart(double epsilon);
    //returns eps_lim_A, eps_lim_B for given function part number
    //limit A is upper strain limit, limit B is lower strain limit
    //we are assuming that compressive strains are negative
    //and tensile strains are positive
    double functionPartStrainLimitA(int number);
    double functionPartStrainLimitB(int number);
    double functionPartStrainLimitA(double epsilon);
    double functionPartStrainLimitB(double epsilon);
    int functionPartDegree(int number);
    int functionPartDegree(double epsilon);

private:
    //epsilon_co, epislon_cu are based on strength of concrete fck
    //function is invoked by setCharacteristicStrength_fck() and constructor Concrete()
    void setStrainThresholdValues(void);

    void setExponentN(void);

    void setModuleOfElasticity_E(void);

    void setMeanStrength_fcm(void);

    //sigma - epsilon (parabolic-linear) function parts
public:
    //1) parabolic relation sig-eps, 3.17, PN-EN 1992-1-1
    //   function takes epsilon (strain value) and returns sigma (stress value)
    //   for 0 < epsilon < epsilon_co
    double parabolicSigmaEpsilon(double epsilon);

    //2) linear relation sig-eps, 3.18, PN-EN 1992-1-1
    // function just returns sigma = f_cd (design value of ultimate compressive strength)
    // for epsilon_co < epsilon < epsilon_cu
    double linearSigmaEpsilon(double epsilon);

    Qt::GlobalColor materialColor(void) {
        return Qt::lightGray;
    }


};

//Exception thrown when user has choosen invalid fck value of concrete
struct InvalidF_ckException : public std::exception
{
   std::string s;
   InvalidF_ckException(std::string ss) : s(ss) {}
   ~InvalidF_ckException() throw () {} // Updated
   const char* what() const throw() { return s.c_str(); }
};


#endif // CONCRETE_H
