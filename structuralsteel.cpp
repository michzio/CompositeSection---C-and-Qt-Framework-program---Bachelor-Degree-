#include "structuralsteel.h"
#include<cmath>

StructuralSteel::StructuralSteel(GradeOfStructuralSteel g) : Material(STRUCTURAL_STEEL)
{
    setGradeOfStructuralSteel(g);
}

void StructuralSteel::setGradeOfStructuralSteel(GradeOfStructuralSteel g) {
    grade = g;
    setYieldStrength();
    setUltimateStrength();
}


GradeOfStructuralSteel StructuralSteel::getGradeOfStructuralSteel(void) {
    return grade;
}

std::string StructuralSteel::getGradeAsString(void)
{
    switch(this->grade)
    {
        case S235:
            return "S235";
        case S275:
            return "S275";
        case S355:
            return "S355";
        case S450:
            return "S450";
        default:
            return "S235"; //default value as grade with minimal strength (safety side)
    }
}

void StructuralSteel::setYieldStrength() {
    switch (grade) {
    case S235:
        f_yk = YIELD_STRENGTH_S235; //[MPa]
        ///TEST ONLY Chen T
        //f_yk = 239; //[MPa]
        ///
        break;
    case S275:
        f_yk = YIELD_STRENGTH_S275; //[MPa]
        ///TEST ONLY - Chen W12x120
        //f_yk = 300; //[MPa]
        /// TEST ONLY
        break;
    case S355:
        f_yk = YIELD_STRENGTH_S355; //[MPa]
        break;
    case S450:
        f_yk = YIELD_STRENGTH_S450; //[MPa]
        break;
    default:
        f_yk = YIELD_STRENGTH_S235; //[MPa]
        break;
    }
}

void StructuralSteel::setUltimateStrength(void) {

    switch (grade) {
    case S235:
        f_uk = ULTIMATE_STRENGTH_S235; //[MPa]
        break;
    case S275:
        f_uk = ULTIMATE_STRENGTH_S275; //[MPa]
        break;
    case S355:
        f_uk = ULTIMATE_STRENGTH_S355; //[MPa]
        break;
    case S450:
        f_uk = ULTIMATE_STRENGTH_S450; //[MPa]
        break;
    default:
        f_uk = ULTIMATE_STRENGTH_S235; //[MPa]
        break;
    }
}

int StructuralSteel::getCharacteristicYieldStrength(void)
{
    return f_yk; //[MPa]
}

int StructuralSteel::getCharacteristicUltimateStrength(void)
{
    return f_uk; //[MPa]
}

double StructuralSteel::getDesignUltimateStrength(void)
{
    return f_uk/GAMMA_M; //[MPa]
}

double StructuralSteel::getDesignYieldStrength(void)
{
    return f_yk/GAMMA_M; //[MPa]
}

double StructuralSteel::getCharacteristicYieldStrainInTension(void)
{
    return (double) f_yk/(E*1000); // f_yk [MPa], E [GPa], 1 GPa = 1000MPa, return [-]
}

double StructuralSteel::getDesignYieldStrainInTension(void)
{
    return (double) getDesignYieldStrength()/(E*1000); // f_yk [MPa], E [GPa], 1 GPa = 1000MPa, return [-]
}

double StructuralSteel::getCharacteristicUltimateStrainInTension()
{
    //return  15*getCharacteristicYieldStrainInTension(); //[-]
    return EPSILON_UK;
}
double StructuralSteel::getDesignUltimateStrainInTension()
{
    //return  15*getCharacteristicYieldStrainInTension(); //[-]
    return getCharacteristicUltimateStrainInTension(); //chen et al.
}


double StructuralSteel::tensileStrainLimit_eps_tu()
{
    //fprintf(stderr, "Structural Steel -> eps_tu: %g\n", getCharacteristicUltimateStrainInTension());
    return getCharacteristicUltimateStrainInTension(); //[-]
}

double StructuralSteel::compressiveStrainLimit_eps_cu()
{
    return (- getCharacteristicUltimateStrainInTension()); //[-]
}

double StructuralSteel::pureCompressionStrainLimit_eps_co()
{
    return compressiveStrainLimit_eps_cu(); //for structural steel we haven't eps_co
                                            //we assume it is equal to eps_cu in that
                                            //way we can get generic algorithm for
                                            //strain profile steerage
}

double StructuralSteel::sigmaEpsilon(double epsilon) {
    double e_yk = getDesignYieldStrainInTension(); //[-]
    double e_uk = getDesignUltimateStrainInTension(); //[-]

    if(epsilon > e_yk && epsilon <= e_uk) {
        return plasticSigmaEpsilon(epsilon); //[MPa]
    } else if(epsilon >= -e_yk && epsilon <= e_yk) {
          return elasticSigmaEpsilon(epsilon); //[MPa]
    } else if(epsilon >= -e_uk && epsilon < -e_yk){
        return (- plasticSigmaEpsilon(std::abs(epsilon))); //[MPa]
    } else {
        //strains are greater than ultimate strain value
        throw StrainOutOfBoundsException("Odkształcenia znajdują się po za zakresem dopuszczalnym.");
    }
}

//function returns stresses based on modul of elasticity multiplied by strains epsilon and 1000 to change units from GPa -> MPa
double StructuralSteel::elasticSigmaEpsilon(double epsilon)
{
    return epsilon*E*1000; // epsilon [-], E [GPa], return [MPa]
}

double StructuralSteel::plasticSigmaEpsilon(double epsilon)
{
    return getDesignYieldStrength(); //[MPa]
}

//function returns current material Module of elasticity
int StructuralSteel::getModuleOfElasticity_E(void)
{
    return E;
}

/**
 * @brief StructuralSteel::plasticTensileSigmaEpsilon
 * @param epsilon - strains must be in range [e_yk, e_uk]
 * @return stress corresponding to strains equal to epsilon
 * We assume here that in plastic range in tensile there are constant
 * stresses equal to yield strength f_yk
 */
double StructuralSteel::plasticTensileSigmaEpsilon(double epsilon)
{
    double e_yk = getDesignYieldStrainInTension(); //[-]
    double e_uk = getDesignUltimateStrainInTension(); //[-]
    if(epsilon < e_yk && epsilon > e_uk)
        throw StrainOutOfBoundsException("Odkształcenia znajdują się poza zakresm przewidzianym dla tej funkcji sigma-epsilon");
    return plasticSigmaEpsilon(epsilon);
}

/**
 * @brief StructuralSteel::plasticCompressiveSigmaEpsilon
 * @param epsilon - strains must be in range [-e_uk, -e_yk]
 * @return
 * We assume here that in plastic range in  compressive (negative strains)
 * there are constant stresses equal to -(yeild_strength) i.e. -f_yk
 */
double StructuralSteel::plasticCompressiveSigmaEpsilon(double epsilon)
{
    double e_yk = getDesignYieldStrainInTension(); //[-]
    double e_uk = getDesignUltimateStrainInTension(); //[-]
    if(epsilon > -e_yk && epsilon < -e_uk)
            throw StrainOutOfBoundsException("Odkształcenia znajdują się poza zakresm przewidzianym dla tej funkcji sigma-epsilon");
    return -plasticSigmaEpsilon(epsilon);
}

/**
 * @brief StructuralSteel::numberOfFunctionParts
 * @return number of function parts from which consist sigma-epsilon relation
 */
int StructuralSteel::numberOfFunctionParts(void)
{
       //sigma-epsilon for reinforcement consist with 3 function parts
       return noOfFunctionParts;
}

/**
 * @brief StructuralSteel::functionPartNumberFor
 * @param epsilon - strain value must be in rane [-eps_uk, eps_uk]
 * @return for given strain value returns corresponding function part
 * number. There are numerated from 0 to numberOfFunctionParts -1
 */
int StructuralSteel::functionPartNumberFor(double epsilon)
{
    double e_yk = getDesignYieldStrainInTension(); //[-]
    double e_uk = getDesignUltimateStrainInTension(); //[-]

    //e_yk, e_uk are positive by default -> tension
    //if they are with - sign i.e negative we have compression
    if(epsilon < - e_uk)
    {
        //epsilon < -e_uk
        throw StrainOutOfBoundsException("Odkształcenia w stali konstrukcyjnej poza zakresem dopuszczalnych odkształceń.");
    } else if(epsilon < -e_yk)
    {
        //range from (-e_yk, -e_uk] compression -> function part f2 (plastic range)
        return 2;
    } else if(epsilon <= e_yk)
    {
        //range from [e_yk, -e_yk]  where [e_yk, 0) - tension
        //and [0, -e_yk] compression -> function part f1 (elastic range)
        return 1;
    } else if(epsilon <= e_uk)
    {
        //range from [e_uk e_yk) tenstion -> function part f0 (plastic range)
        return 0;
    } else {
        //epsilon > e_uk
        throw StrainOutOfBoundsException("Odkształcenia w stali konstrukcyjnej są poza zakresem dopuszczalnych odkształceń.");
    }
}

/**
 * @brief StructuralSteel::functionPart
 * @param epsilon - strain value must be in range [-eps_uk, eps_uk]
 * @return sigma-epsilon relation function part corresponding to given epsilon strain
 * Function part is returned as pointer to StructuralSteel object's method
 */
StructuralSteel::SigmaEpsilon StructuralSteel::functionPart(double epsilon)
{
    //calling overloaded method functionPart(int) passing number generated based on epsilon
    return functionPart(functionPartNumberFor(epsilon));
}

/**
 * @brief StructuralSteel::functionPart
 * @param number - of function part that should be returned [0, n-1]
 * @return Function part is returned as pointer to Structural Steel object's method
 */
StructuralSteel::SigmaEpsilon StructuralSteel::functionPart(int number)
{
    switch(number)
    {
        case 0:
            //epsilon in range [eps_uk, eps_yk) tension, strains positive
            return &StructuralSteel::plasticTensileSigmaEpsilon;
        case 1:
            //epsilon in range [eps_yk, -eps_yk] tension/compression, elastic funcion part
            return &StructuralSteel::plasticSigmaEpsilon;
        case 2:
            //epsilon in range (-eps_yk, -eps_uk] compression, strains negative
            return &StructuralSteel::plasticCompressiveSigmaEpsilon;
        default:
            return NULL;
    }

}

/**
 * @brief StructuralSteel::functionPartStrainLimitA
 * @param number - of function part in stress-strain relation
 * @return
 * Method returns eps_lim_A, k for k = number which is number of
 * function part for which to return upper strain limit of
 * function part interval in stress-strain relation
 */
double StructuralSteel::functionPartStrainLimitA(int number)
{
    double e_yk = getDesignYieldStrainInTension(); //[-]
    double e_uk = getDesignUltimateStrainInTension(); //[-]

    switch(number)
    {
        case 0:
            //epsilon in range [eps_uk, eps_yk) tension, strains positive
            return e_uk;
        case 1:
            //epsilon in range [eps_yk, -eps_yk] tension/compression, elastic function part
            return e_yk;
        case 2:
            //epsilon in range (-eps_yk, -eps_uk] compression, strains negative
            return -e_yk;
    }
}

/**
 * @brief StructuralSteel::functionPartStrainLimitB
 * @param number - of function part in stress-strain relation
 * @return
 * Method returns eps_lim_B,k for k = number which is number of
 * function part for which to return lower strain limit of
 * function part interval in stress-strain relation
 */
double StructuralSteel::functionPartStrainLimitB(int number)
{
    double e_yk = getDesignYieldStrainInTension(); //[-]
    double e_uk = getDesignUltimateStrainInTension(); //[-]

    switch(number) {
        case 0:
            //epsilon in range [eps_uk, eps_yk) tension, strains positive
            return e_yk;
        case 1:
            //epsilon in range [eps_yk, -eps_yk] tension/compression, elastic function part
            return -e_yk;
        case 2:
            //epsilon in range (-eps_yk, -eps_uk] compression, strains negative
            return -e_uk;
    }
}

/**
 * @brief StructuralSteel::functionPartStrainLimitA
 * @param epsilon
 * @return
 * Is overloading of above method matching epsilon
 * to function part in stress-strain relation and next
 * returning strain limit eps_lim_A,k
 */
double StructuralSteel::functionPartStrainLimitA(double epsilon)
{
    double e_uk = getDesignUltimateStrainInTension(); //[-]
    if(epsilon < -e_uk || epsilon > e_uk)
        throw StrainOutOfBoundsException("Odkształcenia poza zakresem dopuszczalnych odkształceń.");
    return functionPartStrainLimitA(functionPartNumberFor(epsilon));
}


/**
 * @brief StructuralSteel::functionPartStrainLimitB
 * @param epsilon
 * @return
 * Is overloading of above method matching epsilon
 * to function part in stress-strain relation and next
 * returning strain limit eps_lim_B,k
 */
double StructuralSteel::functionPartStrainLimitB(double epsilon)
{
    double e_uk = getDesignUltimateStrainInTension(); //[-]
    if(epsilon < -e_uk || epsilon > e_uk)
        throw StrainOutOfBoundsException("Odkształcenia poza zakresem dopuszczalnych odkształceń.");
    return functionPartStrainLimitB(functionPartNumberFor(epsilon));
}

/**
 * @brief StructuralSteel::functionPartDegree
 * @param number - of function part in stress-strain relation
 * @return degree of this function part
 */
int StructuralSteel::functionPartDegree(int number)
{
    switch(number)
    {
        case 0:
           return 0; //in plastic range we are assuming constant stress value
        case 1:
           return 1; //in elastic range there is linear sigma-epsilon relaton
        case 2:
           return 0; //in plastic compressive range we are assuming constant stress value
    }
}

int StructuralSteel::functionPartDegree(double epsilon)
{
    double e_uk = getDesignUltimateStrainInTension(); //[-]
    if(epsilon < -e_uk || epsilon > e_uk)
        throw StrainOutOfBoundsException("Odkształcenia poza zakresem dopuszczalnych odkształceń.");
    return functionPartDegree(functionPartNumberFor(epsilon));
}
