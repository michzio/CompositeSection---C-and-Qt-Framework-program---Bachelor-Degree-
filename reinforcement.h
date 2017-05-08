#ifndef REINFORCEMENT_H
#define REINFORCEMENT_H

#include <cmath>
#include "material.h"

//type of reinforcement steel
enum GradeOfSteel { St3SYb500, RB500W, Bst500S, B500SP };
enum ClassOfSteel { STEEL_A, STEEL_B, STEEL_C };

class Reinforcement : public Material
{
 static constexpr double ERROR_OFFSET = //0.0003;
         0.0001;

public:
    //type definition which is pointer to the method of this object
    //that contains formula to calculate sigma for given function part of sigma-epsilon relation
    typedef double (Reinforcement::*SigmaEpsilon)(double);
    //EX. USAGE
    //Reinforcement::SigmaEpsilon pSigmaEpsilon = &Reinforcement::elasticSigmaEpsilon;
    //albo:
    //Reinforcement::SigmaEpsilon pSigmaEpsilon = pReinforcement->functionPart(num);
    //(pReinforcement->*pSigmaEpsilon)();

 private:
    //FUNCTION PARTS
    const int noOfFunctionParts = 3; //bilinear

    //Yield Strendth limits allowed by PN-EN 1992
    static const int MIN_YIELD_STRENGTH = 400; //[MPa]
    static const int MAX_YIELD_STRENGTH = 600; //[MPa]

    //Modul of elasticity
    static const int Es = 200; //[GPa]

    //Characteristic Strains from PN-EN 1992, appendix C, table C.1
    constexpr static const double STRAIN_e_uk_A = 0.025;
    constexpr static const double STRAIN_e_uk_B = 0.05;
    constexpr static const double STRAIN_e_uk_C = 0.075;

    //Ductility k from PN-EN 1992, appednix C, Table C.1
    constexpr static const double DUCTILITY_A = 1.05;
    constexpr static const double DUCTILITY_B = 1.08;
    constexpr static const double DUCTILITY_C = 1.35;

    constexpr static const double EPSILON_UK = 0.02;

    //Safety coefficients
    constexpr static const double GAMMA_S = //1.15;
            1.00;

    int f_yk; //characteristic value f_yk [MPa]
    double epsilon_uk;
    double ductility; // k = (f_t/f_y)_k
    GradeOfSteel grade;
    ClassOfSteel classOfSteel;

public:
    //constructor takes f_y, grade
    Reinforcement(GradeOfSteel);
    int getCharacteristicYieldStrength(void);
    void setGradeOfSteel(GradeOfSteel);
    GradeOfSteel getGradeOfSteel(void);
    std::string getGradeAsString(void);
    double getDuctility(void);
    double getCharacteristicUltimateStrainInTension(void);
    double getDesignYieldStrength(void);
    double getCharacteristicTensileStrength(void);
    double getDesignTensileStrength(void);
    //function returning tensile and compressive ultimate strain limits of current material
    double tensileStrainLimit_eps_tu(void);
    double compressiveStrainLimit_eps_cu(void);
    double pureCompressionStrainLimit_eps_co(void);

    static ClassOfSteel classFromGradeOfSteel(GradeOfSteel grade);

    int getModuleOfElasticity_E(void);

    //bilinear sigma - epsilon relation
    double sigmaEpsilon(double epsilon);

    //function parts
    int numberOfFunctionParts(void);
    int functionPartNumberFor(double epsilon);
    Reinforcement::SigmaEpsilon functionPart(int number);
    Reinforcement::SigmaEpsilon functionPart(double epsilon);
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
    //based on grade of steel and its class
    void setCharacteristicYieldStrength(int);
    void setDuctilityOfSteel(void);
    void setCharacteristicUltimateStrainInTension(void);
public:
    double elasticSigmaEpsilon(double epsilon);
    double plasticCompressiveSigmaEpsilon(double epsilon);
    double plasticTensileSigmaEpsilon(double epsilon);
    Qt::GlobalColor materialColor(void) {
        return Qt::black;
    }
private:
    double plasticSigmaEpsilon(double epsilon);

};

//Exception thrown when user has choosen invalid fck value of concrete
struct InvalidF_ykException : public std::exception
{
   std::string s;
   InvalidF_ykException(std::string ss) : s(ss) {}
   ~InvalidF_ykException() throw () {} // Updated
   const char* what() const throw() { return s.c_str(); }
};
#endif // REINFORCEMENT_H
