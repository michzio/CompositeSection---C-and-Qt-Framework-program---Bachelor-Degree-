#include "reinforcement.h"

// Es [GPa], sigma [MPa], f_yk, .. [MPas]
// POPRAWIĆ JEDNOSTKI!!!

Reinforcement::Reinforcement(GradeOfSteel g) : Material(REINFORCEMENT)
{
    setGradeOfSteel(g);



}


void Reinforcement::setCharacteristicYieldStrength(int fyk){

    ///TEST ONLY - SECTION WITH <400MPa REINFORCEMENT
    if(fyk < 400 || fyk > 600) throw new InvalidF_ykException("Invalid yield strength f_yk, should be not less than 400MPa and not greater than 600MPa");

    f_yk = fyk;
}

int Reinforcement::getCharacteristicYieldStrength(void)
{
    return f_yk; //[MPa]
}

void Reinforcement::setGradeOfSteel(GradeOfSteel g)
{
    grade = g;
    classOfSteel = Reinforcement::classFromGradeOfSteel(grade);

    switch(g) {

        case St3SYb500:
            setCharacteristicYieldStrength(500); //[MPa]
            ///TEST ONLY - PAPANIKOLOAU SECTION C
            //setCharacteristicYieldStrength(220); //[MPa]
             ///CHEN
             //setCharacteristicYieldStrength(420); //[MPa]
             //setCharacteristicYieldStrength(460); //[MPa]
              //setCharacteristicYieldStrength(300); //[MPa]
            //setCharacteristicYieldStrength(413); //[MPa]
            ///ROSATI
             //setCharacteristicYieldStrength(375); //[MPa]
            /// TEST ONLY
            break;
        case RB500W:
            setCharacteristicYieldStrength(500); //[MPa]
            break;
        case Bst500S:
            setCharacteristicYieldStrength(500); //[MPa]
            break;
        case B500SP:
            setCharacteristicYieldStrength(500); //[MPa]
            break;
        default:
            setCharacteristicYieldStrength(400); //[MPa]
    }

    setDuctilityOfSteel(); //k Table C.1, PN-EN 1992
    setCharacteristicUltimateStrainInTension(); //e_uk, Table C.1, PN-EN 1992
}

GradeOfSteel Reinforcement::getGradeOfSteel(void) {
    return grade;
}

std::string Reinforcement::getGradeAsString(void) {

    switch(grade)
    {
        case St3SYb500:
            return "St3SYb500";
        case RB500W:
            return "RB500W";
        case B500SP:
            return "B500SP";
        case Bst500S:
            return "Bst500S";
    }
    //defualt value
    return "RB500W";
}

void Reinforcement::setDuctilityOfSteel(void) {
    if(classOfSteel == STEEL_A) {
         ductility = DUCTILITY_A; //[-]
    }  else if(classOfSteel == STEEL_B) {
         ductility = DUCTILITY_B; //[-]
    }  else if(classOfSteel == STEEL_C) {
         ductility = DUCTILITY_C; //[-]
    }
}

double Reinforcement::getDuctility(void) {
    return ductility; //[-]
}

void Reinforcement::setCharacteristicUltimateStrainInTension(void) {
   //Old code
    if(classOfSteel == STEEL_A) {
        epsilon_uk = STRAIN_e_uk_A; //[-]
    }  else if(classOfSteel == STEEL_B) {
         epsilon_uk = STRAIN_e_uk_B; //[-]
    }  else if(classOfSteel == STEEL_C) {
         epsilon_uk = STRAIN_e_uk_C; //[-]
    }

    //Overriden epsilon_uk in accordance with Papanikolaou materials
    epsilon_uk = EPSILON_UK;
    ///TEST ONLY - CHEN
    //epsilon_uk = 0.01;
}

double Reinforcement::getCharacteristicUltimateStrainInTension(void) {
    return epsilon_uk; //[-]
}

double Reinforcement::tensileStrainLimit_eps_tu()
{
    return getCharacteristicUltimateStrainInTension();
}

double Reinforcement::compressiveStrainLimit_eps_cu()
{
    return (- getCharacteristicUltimateStrainInTension());
}

double Reinforcement::pureCompressionStrainLimit_eps_co()
{
    return compressiveStrainLimit_eps_cu(); //for reinforcement we haven't define eps_co
                                            //we return eps_cu to make strain profile
                                            //steerage algoithm more generic
}

double Reinforcement::getDesignYieldStrength(void) {
    return f_yk/GAMMA_S; //[MPa]
}

double Reinforcement::getCharacteristicTensileStrength(void) {
    return ductility * f_yk; //[MPa]
}
double Reinforcement::getDesignTensileStrength(void) {
    return ductility *getDesignYieldStrength(); //[MPa]

}

/********************************************************/
double Reinforcement::sigmaEpsilon(double epsilon) {

    double epsilon_yk = (getDesignYieldStrength()/(Es*1000)); //eps [-], yieldStrength [MPa], Es [GPa], 1GPa = 1000 MPa

    if(epsilon > epsilon_yk && epsilon <= epsilon_uk) {
            return plasticSigmaEpsilon(epsilon); //[MPa]
    } else if(epsilon >= -epsilon_yk
              && epsilon <= epsilon_yk ) {
            return elasticSigmaEpsilon(epsilon); //[MPa]
    } else  if(epsilon >= -epsilon_uk && epsilon < -epsilon_yk){
            return ( - plasticSigmaEpsilon(std::abs(epsilon)) ); //[MPa]
    } else {
        if(epsilon > epsilon_uk && epsilon < (epsilon_uk + ERROR_OFFSET)) {
            return plasticSigmaEpsilon(epsilon_uk);
        }
        fprintf(stderr, "epsilon = %g\n", epsilon);
        //strains are greater than ultimate strain value
        throw StrainOutOfBoundsException("Odkształcenia znajdują się po za zakresem dopuszczalnym.");
    }
}

double Reinforcement::elasticSigmaEpsilon(double epsilon)
{
     double sigma = 0;

     sigma = epsilon*Es*1000; // 1 GPa = 1000 MPa!, eps[-], Es [GPa], sigma [MPa]

     return sigma; //[MPa]
}

double Reinforcement::plasticSigmaEpsilon(double epsilon)
{

    double sigma = 0;

    sigma = getDesignYieldStrength(); //[MPa]

    ///TEST ONLY - without hardening branch - delete if you want hardening branch
    return sigma;
    ///TEST ONLY

    //interpolation between e_yk and e_u k
    double epsilon_yk = getDesignYieldStrength()/(Es*1000); //eps [-], yield_strength[MPa], Es [GPa], 1GPa = 1000 MPa
    double f_yd = getDesignYieldStrength(); //[MPa]
    double f_td = getDesignTensileStrength(); //[MPa]

    sigma +=  (f_td - f_yd)/(epsilon_uk - epsilon_yk)
            *(epsilon - epsilon_yk); //[MPa]

    return sigma; //[MPa]
}

double Reinforcement::plasticTensileSigmaEpsilon(double epsilon)
{
    double epsilon_yk = (getDesignYieldStrength()/(Es*1000)); //eps [-], yieldStrength [MPa], Es [GPa], 1GPa = 1000 MPa
    if(epsilon < epsilon_yk || epsilon > epsilon_uk)
        throw StrainOutOfBoundsException("Odkształcenia poza przedziałem dla, którego stosuje się ten fragment funckji sigma-epsilon.");
    //we have positive epsilon value (tension)
    //and return for positive strain also positive strength value
    return plasticSigmaEpsilon(epsilon); //[MPa]
}

double Reinforcement::plasticCompressiveSigmaEpsilon(double epsilon)
{
    double epsilon_yk = (getDesignYieldStrength()/(Es*1000)); //eps [-], yieldStrength [MPa], Es [GPa], 1GPa = 1000 MPa
    if(epsilon > -epsilon_yk || epsilon < -epsilon_uk)
        throw StrainOutOfBoundsException("Odkształcenia poza przedziałem dla, którego stosuje się ten fragment funckji sigma-epsilon.");
    //we have negative epsilon value (compression)
    //and return for negative strain also negative strength value
    //we use this same pirvate method to calculate plastic compressive strength
    //as for plastic tensile strength but we pass absolute epsilon value
    //and return strength with negative sign (-)!
     return ( - plasticSigmaEpsilon(std::abs(epsilon)) ); //[MPa]
}

/*******************************************************/
ClassOfSteel Reinforcement::classFromGradeOfSteel(GradeOfSteel grade)
{
    switch (grade) {
    case St3SYb500:
        return STEEL_A;
    case RB500W:
        return STEEL_B;
    case Bst500S:
        return STEEL_B;
    case B500SP:
        return STEEL_C;
    default:
        return STEEL_A;
    }
}

int Reinforcement::getModuleOfElasticity_E(void)
{
    return Es; //[GPa]
}

/**
 * @brief Reinforcement::numberOfFunctionParts
 * @return number of function parts from which consist sigma-epsilon relation
 */
int Reinforcement::numberOfFunctionParts(void)
{
   //sigma-epsilon for reinforcement consists with 3 function parts
   return noOfFunctionParts;
}

/**
 * @brief Reinforcement::functionPartNumberFor
 * @param epsilon - strain value must be in range [-eps_uk, eps_uk]
 * @return for given strain value returns corresponding function part number
 * There are numerated starting from 0 to numberOfFunctionParts - 1
 */
int Reinforcement::functionPartNumberFor(double epsilon)
{
     double epsilon_yk = getDesignYieldStrength()/(Es*1000); //eps [-], yield_strength[MPa], Es [GPa], 1GPa = 1000 MPa

     //eps_uk, eps_yk are positive by default -> tension
     //if they are with - sign i.e. negative we have compression
     if(epsilon < -epsilon_uk)
     {
         throw StrainOutOfBoundsException("Odkształcenia w stali zbrojeniowej poza zakresem dopuszczalnych odkształceń.");
     } else if(epsilon < -epsilon_yk)
     {
        //range from (-eps_yk, -eps_uk] compression -> function part f2 (plastic)
        return 2;
     } else if(epsilon <= epsilon_yk)
     {
          //range from [eps_yk, -eps_yk] where [eps_yk, 0) - tension
          //and [0, -eps_yk] - compresion -> function part f1 (elastic)
         return 1;
     } else if(epsilon <= epsilon_uk)
     {
          //range from [eps_uk, eps_yk) tension  -> function part f0 (plastic)
          return 0;
     } else {
         //epsilon > epsilon_uk
         throw StrainOutOfBoundsException("Odkształcenia w stali zbrojeniowej poza zakresem dopuszczalnych odkształceń.");
     }

}

/**
 * @brief Reinforcement::functionPart
 * @param epsilon - strain value must be in range [-eps_uk, eps_uk]
 * @return sigma-epsilon relation function part corresponding to given epsilon strain
 * Function part is returned as pointer to Reinforcement object method
 */
Reinforcement::SigmaEpsilon Reinforcement::functionPart(double epsilon)
{
    //calling overloaded method functionPart(int) passing number generated based on epsilon
    return functionPart(functionPartNumberFor(epsilon));
}

/**
 * @brief Reinforcement::functionPart
 * @param number - of function part that should be returned [0, n-1]
 * @return Function part is returned as pointer to Reinforcement object method
 */
Reinforcement::SigmaEpsilon Reinforcement::functionPart(int number)
{
    switch(number)
    {
        case 0:
            //epsilon in range [eps_uk, eps_yk) tension, strains positive
            return &Reinforcement::plasticTensileSigmaEpsilon;
        case 1:
            //epsilon in range [eps_yk, -eps_yk] tension/compression, elastic function part
            return &Reinforcement::elasticSigmaEpsilon;
        case 2:
            //epsilon in range (-eps_yk, -eps_uk] compression, strains negative
            return &Reinforcement::plasticCompressiveSigmaEpsilon;
        default:
            return NULL;
    }
}

/**
 * @brief Reinforcement::functionPartStrainLimitA
 * @param number - of function part in stress-strain relation
 * @return
 * Method returns eps_lim_A,k for k = number which is number of
 * function part for which to return upper strain limit of
 * function part interval in stress-strain relation
 */
double Reinforcement::functionPartStrainLimitA(int number)
{
    double epsilon_yk = getDesignYieldStrength()/(Es*1000); //eps [-], yield_strength[MPa], Es [GPa], 1GPa = 1000 MPa

    switch(number) {
        case 0:
            //epsilon in range [eps_uk, eps_yk) tension, strains positive
            return epsilon_uk;
        case 1:
            //epsilon in range [eps_yk, -eps_yk] tension/compression, elastic function part
            return epsilon_yk;
        case 2:
            //epsilon in range (-eps_yk, -eps_uk] compression, strains negative
            return -epsilon_yk;
    }
}

/**
 * @brief Reinforcement::functionPartStrainLimitB
 * @param number - of function part in stress-strain relation
 * @return
 * Method returns eps_lim_B,k for k = number which is number of
 * function part for which to return lower strain limit of
 * function part interval in stress-strain relation
 */
double Reinforcement::functionPartStrainLimitB(int number)
{
    double epsilon_yk = getDesignYieldStrength()/(Es*1000); //eps [-], yield_strength[MPa], Es [GPa], 1GPa = 1000 MPa

    switch(number) {
        case 0:
            //epsilon in range [eps_uk, eps_yk) tension, strains positive
            return epsilon_yk;
        case 1:
            //epsilon in range [eps_yk, -eps_yk] tension/compression, elastic function part
            return -epsilon_yk;
        case 2:
            //epsilon in range (-eps_yk, -eps_uk] compression, strains negative
            return -epsilon_uk;
    }
}

/**
 * @brief Reinforcement::functionPartStrainLimitA
 * @param epsilon
 * @return
 * Is overloading of above method matching epsilon
 * to function part in stress-strain relation and next
 * returning strain limit  eps_lim_A,k
 */
double Reinforcement::functionPartStrainLimitA(double epsilon)
{
    if(epsilon < -epsilon_uk || epsilon > epsilon_uk)
        throw StrainOutOfBoundsException("Odkształcenia poza zakresem dopuszczalnych odkształceń.");
    return functionPartStrainLimitA(functionPartNumberFor(epsilon));
}

/**
 * @brief Reinforcement::functionPartStrainLimitB
 * @param epsilon
 * @return
 * Is overloading of above method matching epsilon
 * to function part in stress-strain relation and next
 * returining strain limit eps_lim_B,k
 */

double Reinforcement::functionPartStrainLimitB(double epsilon)
{
    if(epsilon < -epsilon_uk || epsilon > epsilon_uk)
        throw StrainOutOfBoundsException("Odkształcenia poza zakresem dopuszczalnych odkształeń.");
    return functionPartStrainLimitB(functionPartNumberFor(epsilon));
}

/**
 * @brief Reinforcement::functionPartDegree
 * @param number - of function part in stress-strain relation
 * @return degree of this function part
 */

int Reinforcement::functionPartDegree(int number)
{
    switch(number)
    {
        case 0:
            return 1; //in plastic range we are assuming linear plastic hardening
        case 1:
            return 1;
        case 2:
            return 1; //in plastic range we are assuming linear plastic hardening
    }
}

int Reinforcement::functionPartDegree(double epsilon)
{
    if(epsilon < -epsilon_uk || epsilon > epsilon_uk)
        throw StrainOutOfBoundsException("Odkształcenia poza zakresem dopuszczalnych odkształceń.");
    return functionPartDegree(functionPartNumberFor(epsilon));
}
