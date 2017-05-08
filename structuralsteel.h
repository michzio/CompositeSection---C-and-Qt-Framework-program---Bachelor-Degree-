#ifndef STRUCTURALSTEEL_H
#define STRUCTURALSTEEL_H

#include "material.h"

//grade of steel based on Table 3.1 PN-EN 1993-1-1
//and definition of constants for f_y, f_u
enum GradeOfStructuralSteel { S235, S275, S355, S450 };

class StructuralSteel : public Material
{
 public:
    //type definition which is pointer to the method of this object
    //that contains formula to calculate sigma for given function part
    //of sigma-epsilon relation
    typedef double (StructuralSteel::*SigmaEpsilon)(double);
private:
    //FUNCTION PARTS
    const int noOfFunctionParts = 3; //bilinear

    //Yield and Ultimate strength limits allowed by PN-EN 1993
    static const int YIELD_STRENGTH_S235 = 235; //[MPa]
    static const int YIELD_STRENGTH_S275 = 275; //[MPa]
    static const int YIELD_STRENGTH_S355 = 355; //[MPa]
    static const int YIELD_STRENGTH_S450 = 440; //[MPa]
    static const int ULTIMATE_STRENGTH_S235 = 360; //[MPa]
    static const int ULTIMATE_STRENGTH_S275 = 430; //[MPa]
    static const int ULTIMATE_STRENGTH_S355 = 510; //[MPa]
    static const int ULTIMATE_STRENGTH_S450 = 550; //[MPa]

    constexpr static const double EPSILON_UK = 0.02;

    //Safety coefficients
    const double GAMMA_M = //1.10; //[-]
                            1.00;
    //Module of elasticity
    static const int E = 210; //[GPa]
                        //200; //[GPa] chen et al
    int f_yk; //characteristic value of yield strength [MPa]
    int f_uk; //characteristic value of ultimate strength [MPa]
    GradeOfStructuralSteel grade;
public:
    StructuralSteel(GradeOfStructuralSteel);
    void setGradeOfStructuralSteel(GradeOfStructuralSteel);
    GradeOfStructuralSteel getGradeOfStructuralSteel(void);
    std::string getGradeAsString(void);
    int getCharacteristicYieldStrength(void);
    int getCharacteristicUltimateStrength(void);
    double getDesignYieldStrength(void);
    double getDesignUltimateStrength(void);
    double getCharacteristicUltimateStrainInTension(void);
    double getCharacteristicYieldStrainInTension(void);
    double getDesignYieldStrainInTension(void);
    double getDesignUltimateStrainInTension(void);
    //function returning tensile and compressive ultimate strain limits of current material
    double tensileStrainLimit_eps_tu(void);
    double compressiveStrainLimit_eps_cu(void);
    double pureCompressionStrainLimit_eps_co(void);

    //overriding abstract method form Material class
    int getModuleOfElasticity_E();

    //bilinear sigma - epsilon fig. 5.8, PN-EN 1993
    double sigmaEpsilon(double);

    //function parts
    int numberOfFunctionParts(void);
    int functionPartNumberFor(double epsilon);
    StructuralSteel::SigmaEpsilon functionPart(int number);
    StructuralSteel::SigmaEpsilon functionPart(double epsilon);
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
    void setYieldStrength();
    void setUltimateStrength();
    double plasticSigmaEpsilon(double epsilon);
public:
    double elasticSigmaEpsilon(double epsilon);
    double plasticTensileSigmaEpsilon(double epsilon);
    double plasticCompressiveSigmaEpsilon(double epsilon);

    Qt::GlobalColor materialColor(void) {
        return Qt::darkGray;
    }

};

#endif // STRUCTURALSTEEL_H
