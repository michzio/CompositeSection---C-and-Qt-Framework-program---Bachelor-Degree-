#include "concrete.h"
#include<exception>
#include<cmath> //std::pow
#include "flags.h"

Concrete::Concrete(int fck) : Material(CONCRETE)
{
    //fck = 21; ///TEST ONLY
    setCharacteristicStrength_fck(fck);
}

void Concrete::setCharacteristicStrength_fck(int fck)
{
    if(fck > 90 || fck <12 )
        throw new InvalidF_ckException("Invalid f_ck value.");
    //setting characteristic strength of concrete
    f_ck = fck; //[MPa]

    //calculating strain threshold values depending on the strength value
    setStrainThresholdValues();
    setExponentN();
    setMeanStrength_fcm();
    setModuleOfElasticity_E();
}

int Concrete::getCharacteristicStrength_fck()
{
    //we return negative value as it is compressive strength
    return -f_ck; //[MPa]
}


double Concrete::getDesignStrength_fcd(void)
{
    //we return negative value as it is compressive strength
    //we use formula 3.15, PN-EN 1992-1-1
    return -(ALFA_CC*f_ck)/GAMMA_C; //[MPa]
}

double Concrete::getUltimateStrain_ecu(void)
{
    return epsilon_cu; //[-]
}

double Concrete::getPureCompressionStrain_eco()
{
    return epsilon_co; //[-]
}

double Concrete::tensileStrainLimit_eps_tu()
{
    return 0; //we assume that concrete doesn't bear tensile stresses
}

double Concrete::compressiveStrainLimit_eps_cu()
{
    return getUltimateStrain_ecu();
}

double Concrete::pureCompressionStrainLimit_eps_co()
{
    return getPureCompressionStrain_eco();
}

void Concrete::setStrainThresholdValues()
{
    if(f_ck > 50) { //[MPa]
        //if characteristic strength of concrete is greater than 50MPa
        //we use formules from Table 3.1 in PN-EN 1992-1-1
        epsilon_co = -(2.0 + 0.085*std::pow((f_ck - 50), 0.53))/1000;
        epsilon_cu = -(2.6 + 35*std::pow((0.01*(90-f_ck)), 4))/1000;
    } else {
        epsilon_co = EPSILON_CO_LESS_THAN_50;
        epsilon_cu = EPSILON_CU_LESS_THAN_50;
    }
}

double Concrete::getExponentN(void)
{
    return n; //[-]
}

void Concrete::setExponentN(void) {
    if(f_ck > 50) { //[MPa]
        n = 1.4 + 23.4*std::pow(0.01*(90-f_ck), 4); //[-]
    } else {
        n = 2.0; //[-]
    }
}

void Concrete::setModuleOfElasticity_E(void) {
     //method sets module of elasticity using formula in Table 3.1, PN-EN 1992-1-1
     Ecm = (int) 22*std::pow( (0.1 *f_cm), 0.3); // f_cm in [MPa] and resultant Ecm in [GPa]

}

int Concrete::getModuleOfElasticity_E(void) {
    if(DEBUG)
        fprintf(stderr, "Concrete Modulus of Elasticity Ecm = %d\n", Ecm);
    return Ecm; //[GPa]
}

void Concrete::setMeanStrength_fcm(void) {
    //method sets mean strength of concrete using formula in Table 3.1, PN-EN 1992-1-1
    f_cm = f_ck + 8; //[MPa]
}

int Concrete::getMeanStrength_fcm(void)
{
   return -f_cm; //[MPa]
}

double Concrete::sigmaEpsilon(double epsilon) {

    ///TEST ONLY - N, Mx, My
    //return 1;
    ///TEST ONLY - N, Mx, My
    //strains under compression are negative ( - )!!!
    if(epsilon < epsilon_co) {
        return linearSigmaEpsilon(epsilon);  //[MPa]
    } else if(epsilon < 0) {
        return parabolicSigmaEpsilon(epsilon); //[MPa]
    } else {
        return 0;
    }
}

double Concrete::linearSigmaEpsilon(double epsilon)
{
    if(epsilon < epsilon_cu || epsilon > epsilon_co)
        throw StrainOutOfBoundsException("This strains are not in linear relation with stresses. WRONG FUNCTION PART.");
    return getDesignStrength_fcd(); //[MPa]
}

double Concrete::parabolicSigmaEpsilon(double epsilon) {

     if(epsilon < epsilon_co || epsilon > 0)
         throw StrainOutOfBoundsException("This strains are not in parabolic relation with stresses. WRONG FUNCTION PART.");
      double sigma = 0;
      //To calculate compressive stress in parabolic part
      //of sigma-epsilon relation we use formula 3.17, PN-EN 1992
      double f_cd = getDesignStrength_fcd(); //[MPa]

      sigma = f_cd*(1 - std::pow((1-epsilon/epsilon_co), n));  //[MPa]

      return sigma; //[MPa]
}

/**
 * @brief Concrete::numberOfFunctionParts
 * @return - number of sigma epsilon function parts
 * For parabolic-linear relation in PN-EN 1992 there are 2 function parts
 */
int Concrete::numberOfFunctionParts(void)
{
    return noOfFunctionParts;
}

/**
 * @brief Concrete::functionPartNumberFor
 * @param epsilon - strain value must be in range [eps_cu, 0]
 * @return  For given strain value returns corresponding function part number
 * There are numerated starting from 0 to numberOfFunctionParts - 1
 *
 */
int Concrete::functionPartNumberFor(double epsilon)
{
    if(epsilon < epsilon_cu)
    {
        throw StrainOutOfBoundsException("Odkształcenia przekraczając wartość graniczną.");
    } else if(epsilon < epsilon_co)
    {
        return 1; //linear function part f1
    } else if(epsilon <= 0) {
        return 0; //parabolic function part f0
    } else {
        //tension
        throw StrainOutOfBoundsException("Odkształcenia rozciągane pomijane w betonie.");
    }
}

/**
 * @brief Concrete::functionPart
 * @param epsilon - strain value must be in range [esp_cu,0]
 * @return sigma-epsilon relation function part corresponding to given epsilon strain
 * Function part is returned as pointer to Concrete object method
 */
Concrete::SigmaEpsilon Concrete::functionPart(double epsilon)
{
    return functionPart(functionPartNumberFor(epsilon));
}

/**
 * @brief Concrete::functionPart
 * @param number - of function part that should be returned [0, n-1]
 * @return  Function part is returned as pointer to Concrete object method
 */
Concrete::SigmaEpsilon Concrete::functionPart(int number)
{
    switch(number) {
        case 0:
            return &Concrete::parabolicSigmaEpsilon;
        case 1:
            return &Concrete::linearSigmaEpsilon;
        default:
            return NULL;
    }
}

/**
 * @brief Concrete::functionPartStrainLimitA
 * @param number - of function part in stress-strain relation
 * @return
 * Method returns eps_lim_A,k for k = number which is number of
 * function part for which to return upper strain limit of
 * function part interval in stress-strain relation
 */
double Concrete::functionPartStrainLimitA(int number)
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

/**
 * @brief Concrete::functionPartStrainLimitB
 * @param number - of function part in stress-strain relation
 * @return
 * Method returns esp_lim_B,k for k = number which is number of
 * function part for which to return lower strain limit of
 * function part interval in stress-strain relation
 */
double Concrete::functionPartStrainLimitB(int number)
{
    switch(number)
    {
        case 0:
            //parabolic sigma-epsilon, function part f0
            return epsilon_co;
        case 1:
            //linear sigma-epsilon, function part f1
            return epsilon_cu;
    }

}

/**
 * @brief Concrete::functionPartStrainLimitA
 * @param epsilon
 * @return
 * Is overloading of above method matching epsilon
 * to function part in stress-strain relation and next
 * returning strain limit eps_lim_A,k
 */
double Concrete::functionPartStrainLimitA(double epsilon)
{
    if(epsilon < epsilon_cu || epsilon > 0)
        throw StrainOutOfBoundsException("Odkształcenia poza zakresem dopuszczalnych odkształceń");
    return functionPartStrainLimitA(
                functionPartNumberFor(epsilon));
}

/**
 * @brief Concrete::functionPartStrainLimitB
 * @param epsilon
 * @return
 * Is overloading of above method matching epsilon
 * to function part in stress-strain relation and next
 * returning strain limit eps_lim_B,k
 */
double Concrete::functionPartStrainLimitB(double epsilon)
{
    if(epsilon < epsilon_cu || epsilon > 0)
        throw StrainOutOfBoundsException("Odkształcenia poza zakresem dopuszczalnych odkształceń");
    return functionPartStrainLimitB(
                functionPartNumberFor(epsilon));
}

/**
 * @brief Concrete::functionPartDegree
 * @param number - of function part in stress-strain relation
 * @return degree of this function part
 */
int Concrete::functionPartDegree(int number)
{
    switch (number) {
        case 0:
            return 2;
        case 1:
            return 0;
    }
}

int Concrete::functionPartDegree(double epsilon)
{
    if(epsilon < epsilon_cu || epsilon > 0)
        throw StrainOutOfBoundsException("Odkształcenia poza zakresem dopuszczalnych odkształceń");
    return functionPartDegree(functionPartNumberFor(epsilon));
}
