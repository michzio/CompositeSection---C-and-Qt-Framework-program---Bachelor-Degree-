#include "sectionsolver.h"
#include "strainprofile.h"
#include "flags.h"
#include "vertex3d.h"
#include <QtAlgorithms>

SectionSolver::SectionSolver() : Section()
{
    //while constructing new SectionSolver object we assume
    //that local coordinate system of section (strainprofile)
    //isn't rotated in relation to initial central (global) section coordinate system
    rotationAngle = 0.0;
}

/**
 * @brief SectionSolver::setRotationAngle
 * @param omega - new strain profile rotation angle to local coordinate system while performing iterative stress intergration process
 */
void SectionSolver::setRotationAngle(double omega)
{
    //setting current rotation angle of strain profile in relation to this section
    rotationAngle = omega;
    //after setting new angle we resset strain limits in which strain steerage is performed while performing stress integration
    compressionLimit = 0;
    yCoordinateOfCompressionLimit = 0;
    pureCompressionLimit = 0;
    yCoordinateOfPureCompressionLimit = 0;
    tensionLimit = 0;
    yCoordinateOfTensionLimit = 0;
    //above values will be update on its reading try
}

/**
 * @brief SectionSolver::getRotationAngle
 * @return currently set rotation angle of strain profile
 */
double SectionSolver::getRotationAngle(void)
{
    return rotationAngle;
}

/**
 * @brief SectionSolver::findTensileeStrainLimit
 * @param eps_tu - pointer to place where tensile strain limit will be saved
 * @param y_tu - pointer to place where y coordinate of where this tensile strain limit is achived on section height will be saved
 * Method returns by pointer parameters values of eps_tu and y_tu
 * calculated ealier and stored in properties tensionLimit and yCoordinateOfTensionLimit
 * or if thay are null calls ealier method calculateTensileAndCompressiveStrainLimits()
 *
 */
void SectionSolver::findTensileStrainLimit(double *eps_tu,
                                           double *y_tu)
{

    ///TEST ONLY
    /// *eps_tu = 0.01;
    //double fiber_y_tu = this->getFiberGroups()[0]->yMinInCoordinateSystem(this->getGeometricalCenter(), rotationAngle);
    ///fprintf(stderr, "Fibergrup y_min: %g\n", fiber_y_tu);
    ///  return;
    /// TEST ONLY

    if(tensionLimit == 0
            && compressionLimit == 0)
        calculateTensileAndCompressiveStrainLimits();

    //return tension strain limit point
    *eps_tu = tensionLimit;
    *y_tu = yCoordinateOfTensionLimit;
}

/**
 * @brief SectionSolver::findCompressiveStrainLimit
 * @param eps_cu - pointer to place where compressive strain limit will be saved
 * @param y_cu - pointer to place where y coordinate of where this compressive strain limit is achived on section height will be saved
 * Method returns by pointer parameters values of eps_cu and y_cu
 * calculated ealier and stored in properties compressionLimit and yCoordinateOf CompressionLimit
 * or if thay are null call ealier method calculateTensileAndCompressiveStrainLimits()
 */
void SectionSolver::findCompressiveStrainLimit(double *eps_cu,
                                               double *y_cu)
{
    if(tensionLimit == 0
             && compressionLimit == 0)
             calculateTensileAndCompressiveStrainLimits();
    //return ultimate compression strain limit point
    *eps_cu = compressionLimit;
    *y_cu = yCoordinateOfCompressionLimit;
}

/**
 * @brief SectionSolver::findPureCompressionStrainLimit
 * @param eps_co - pointer to place where pure compression strain limit will be saved
 * @param y_co - poiter to place where y coordinate of where this pure compression strain limit is achived on section height will be saved
 * Method returns eps_co and y_co calculated by
 * calculatePureCompressionStrainLimit() method
 */
void SectionSolver::findPureCompressionStrainLimit(double *eps_co, double *y_co)
{
    if(pureCompressionLimit == 0) {
        if(compressionLimit == 0)
            calculateTensileAndCompressiveStrainLimits();
    }

    //return pure compression strain limit point
    *eps_co = pureCompressionLimit;
    *y_co = yCoordinateOfPureCompressionLimit;
}


/**
 * @brief SectionSolver::calculateTensileAndCompressiveStrainLimits
 * 1. pass through each component of section (surface, line, fibergroup) that allows tensile strains > 0 (not equal)
 *    steel, reinforcement, etc.
 *    IMPORTAN! If there isn't any such component we assume that eps_tu = 0 for y_tu = y_MIN! and pass for this point
 *    through each compressive component of section (surface, line fibergroup) that allow compressive strains < 0
 * 2. for each tensile component pass through each compressive component (surface, line fibergroup) that allow compressive
 *    strains < 0 IMPORTANT! skip components having y_max_j <= y_min_i (we assume compression on upper edge of section and tension on lower edge of section)
 * 3. we take into account two points (eps_tu_i, y_min_i) (tensile component) and (eps_cu_j, y_max_j)
 * 4. having 2 strain limits (points) first tensile strain limit and second compressive strain limit
 *    from consecutive component we extends line passing through this points to intersect with
 *    y_MIN and y_MAX coordinate limits of section. In this way we achive allowed tensile strain limit for
 *    y_MIN (section lower limit) and compressive strain limit for y_MAX (section upper limit)
 *    In each pass we compare resultant compressive strain limits for y_MAX and  tensile strain limits for y_MIN
 *    and search this with the lowest absolute values (taking into account that which are closer to vertical axis of strain profile)
 *    IMPORTAN! we assume that section has compression at the top of section and tension at the bottom
 *    i.e. we only take into account eps_cu_MAX > 0 in y_MAX or only take into account y_max_j > y_min_i
 *         and only eps_tu_MIN <= 0
 * 5. Moreover if current (eps_cu_MAX, y_MAX) and (eps_tu_MIN, y_MIN) will be choosen as limiting values
 *    we must check (eps_co_j, y_co_j) whether isn't exceeded for this choosen limits eps_cu_MAX and eps_tu_MAX
 *    tmp_eps_co < eps_co_j -> eps_co_j is exceeded we need to make adjustment. We do this adjustment by
 *    choosing limiting values as current eps_tu and new eps_co_j! and calculatind for this two points
 *    new_eps_cu_j for y_cu_j (y_max_j of component) and replacing eps_cu_j by this new_eps_cu_j as assignment
 *    to eps_cu
 * 5. At the end we has obtained three points (eps_tu, y_tu) = (eps_tu_i, y_min_i) and (eps_cu, y_cu) = (eps_cu_j or new_eps_cu_j, y_max_j)
 *     and (eps_co, y_co) = (eps_co_j, y_co_j)
 *     First is limiting tensile strain limit in section and coressponding y coordinate where this strain limit must be passed
 *     Second is limiting compressive strain limit in section and corresponding y coordinate where this strain limit must be passed
 *     Third is limiting pure compression strain limit in section and corresponding y coordinate where this strain limit must be passed
 * 6. This 3 points (eps_tu, y_tu) and (eps_cu, y_cu) guarantee that tensile and compressive limits
 *    won't be exceeded in any of section components i.e. for such defined strain profile
 *    (eps_tu_i, y_min_i) and (eps_cu_i, y_max_i) will be laying inside this strain profiles
 */
void SectionSolver::calculateTensileAndCompressiveStrainLimits(void)
{
    double eps_tu = 0.0; //ultimate tensile strain in section limiting strain profile
    double y_tu = 0.0; //corresponding to tensila strain y coordinate in current rotated by angle omega local coordinate system
    double eps_cu = 0.0; //ultimate compressive strain in section limiting strain profile
    double y_cu = 0.0; //corresponding to compressive strain y coordinate in current rotated by angle omega local coordinate system
    double eps_co = 0.0; //corresponding pure compression strain limit
    double y_co = 0.0; //y coordinate in rotated by angle omega local coordinate system which is the position of eps_co

    //centroid of local coordinate system (with strain profile) rotation angle
    Point *centroid = this->getGeometricalCenter();
    //property double rotationAngle!
    double y_MIN = this->yMinInCoordinateSystem(centroid, rotationAngle);
    double y_MAX = this->yMaxInCoordinateSystem(centroid, rotationAngle);
    //helper strains for y_MIN and y_MAX coordinate
    double eps_tu_MIN = 0.0;
    double eps_cu_MAX = 0.0;

    //we are passing through each surface component
    for(int i = 0; i < surfaces.length(); i++)
    {
        //skip this surface component if it is opening
        if(surfaces[i]->getIsOpening()) continue;

        //foreach ith surface we get its y lower coordinate limit
        double y_min_i = surfaces[i]->yMinInCoordinateSystem(centroid, rotationAngle);
        //for this component material we get ultimate tensile strain
        double eps_tu_i  = surfaces[i]->getMaterial()->tensileStrainLimit_eps_tu();
        if(DEBUG)
            fprintf(stderr, "Next surface eps_tu_i = %g\n", eps_tu_i);
        //check whether eps_tu_i (ultimate tensile strain of ith component are > 0
        //if not skip this point (this component)
        if(eps_tu_i <= 0) continue;

        //when ultimate tensile strain is greater than 0 we have for ith component point
        //(eps_tu_i, y_min_i) fixed and are passing through each component and serach its ultimate
        //compressive strain eps_cu_j such that having this two points (eps_tu_i, y_min_i)
        //and (eps_cu_j, y_max_j) we obtain simultanously minimal (as absolute value) (eps_tu_MIN, y_MIN) and
        // (eps_cu_MAX, y_MAX)
        //we pass through each component (surface, line, fiber group) as searching (eps_cu_j, y_max_j)

        //loop through each surface
        for(int j=0; j < surfaces.length(); ++j)
        {
            //skip this j-th surface component if it is opening
            if(surfaces[j]->getIsOpening()) continue;

            //foreach jth surface we get its y upper coordinate limit
            double y_max_j = surfaces[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = surfaces[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e ith component lower edge

            //we calculate direction coeffictient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;
            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for component S%d and S%d: (%g, %g)\n", i, j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for component S%d and S%d: (%g, %g)\n", i, j, tmp_eps_cu_MAX, y_MAX);
            }
            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and comressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = surfaces[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = surfaces[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't cause
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(DEBUG)
                    fprintf(stderr, "eps_co_from_tu_cu = %g \n", eps_co_from_tu_cu);
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced form it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j); //should be negative
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }

            }
        }

        //loop through each line element

        //next we do analogues iteration for line elements
        for(int j = 0; j < lines.length(); ++j)
        {
            //foreach jth line element we get its y upper coordinate limit
            double y_max_j = lines[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = lines[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e. ith component lower edge

            //we calculate direction coefficient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;
            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for component S%d and L%d: (%g, %g)\n", i, j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for component S%d and L%d: (%g, %g)\n", i, j, tmp_eps_cu_MAX, y_MAX);
            }
            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and compressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = lines[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = lines[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't casue
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced from it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j);
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX  = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }
            }
        }

        //loop through each fiber group

        //next we do analogues iteration for fiber groups
        for(int j = 0; j < fiberGroups.length(); ++j)
        {
            //foreach jth fiber group we get its y upper coordinate limit
            double y_max_j = fiberGroups[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = fiberGroups[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e. ith component lower edge

            //we calculate direction coefficient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;

            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for component S%d and FG%d: (%g, %g)\n", i, j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for component S%d and FG%d: (%g, %g)\n", i, j, tmp_eps_cu_MAX, y_MAX);
            }

            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and comressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = fiberGroups[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = fiberGroups[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't cause
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced form it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j); //should be negative
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }

            }

        }

    }

    //analogous iterative process for line element components
    for(int i=0; i < lines.length(); i++)
    {
        //foreach ith line element we get its y lower coordinate limit
        double y_min_i = lines[i]->yMinInCoordinateSystem(centroid, rotationAngle);
        //for this component material we get ultimate tensile strain
        double eps_tu_i = lines[i]->getMaterial()->tensileStrainLimit_eps_tu();
        //check whether eps_tu_i (ultimate tensile strain of ith component are > 0
        //if not skip this point (this component)
        if(eps_tu_i <= 0) continue;

        //loop through each surface
        for(int j=0; j < surfaces.length(); ++j)
        {
            //skip if j-th surface is opening
            if(surfaces[j]->getIsOpening()) continue;

            //foreach jth surface we get its y upper coordinate limit
            double y_max_j = surfaces[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = surfaces[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e ith component lower edge

            //we calculate direction coeffictient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;
            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for component S%d and S%d: (%g, %g)\n", i, j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for component S%d and S%d: (%g, %g)\n", i, j, tmp_eps_cu_MAX, y_MAX);
            }
            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE3  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                 //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and comressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = surfaces[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = surfaces[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't cause
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced form it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j); //should be negative
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }
            }
        }

        //loop through each line element

        //next we do analogues iteration for line elements
        for(int j = 0; j < lines.length(); ++j)
        {
            //foreach jth line element we get its y upper coordinate limit
            double y_max_j = lines[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = lines[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e. ith component lower edge

            //we calculate direction coefficient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;

            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for component S%d and L%d: (%g, %g)\n", i, j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for component S%d and L%d: (%g, %g)\n", i, j, tmp_eps_cu_MAX, y_MAX);
            }
            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE3  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and comressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = lines[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = lines[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't cause
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced form it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j); //should be negative
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }
            }
        }

        //loop through each fiber group

        //next we do analogues iteration for fiber groups
        for(int j = 0; j < fiberGroups.length(); ++j)
        {
            //foreach jth fiber group we get its y upper coordinate limit
            double y_max_j = fiberGroups[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = fiberGroups[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e. ith component lower edge

            //we calculate direction coefficient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;
            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for component S%d and FG%d: (%g, %g)\n", i, j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for component S%d and FG%d: (%g, %g)\n", i, j, tmp_eps_cu_MAX, y_MAX);
            }
            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE3  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and comressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = fiberGroups[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = fiberGroups[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't cause
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced form it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j); //should be negative
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }
            }

        }

    }

    //analogous iterative process for fiber group components
    for(int i=0; i < fiberGroups.length(); i++)
    {
        //foreach ith fiber group we get its y lower coordinate limit
        double y_min_i = fiberGroups[i]->yMinInCoordinateSystem(centroid, rotationAngle);
        //for this component material we get ultimate tensile strain
        double eps_tu_i = fiberGroups[i]->getMaterial()->tensileStrainLimit_eps_tu();
        //check whether eps_tu_i (ultimate tensile strain of ith component are > 0
        //if not skip this point (this component)
        if(eps_tu_i <= 0) continue;

        //loop through each surface
        for(int j=0; j < surfaces.length(); ++j)
        {
            //skip if j-th surface is opening
            if(surfaces[j]->getIsOpening()) continue;

            //foreach jth surface we get its y upper coordinate limit
            double y_max_j = surfaces[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = surfaces[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e ith component lower edge

            //we calculate direction coeffictient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;

            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for component S%d and S%d: (%g, %g)\n", i, j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for component S%d and S%d: (%g, %g)\n", i, j, tmp_eps_cu_MAX, y_MAX);
            }
            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE2  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and comressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = surfaces[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = surfaces[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't cause
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced form it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j); //should be negative
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }
            }
        }

        //loop through each line element

        //next we do analogues iteration for line elements
        for(int j = 0; j < lines.length(); ++j)
        {
            //foreach jth line element we get its y upper coordinate limit
            double y_max_j = lines[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = lines[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e. ith component lower edge

            //we calculate direction coefficient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;

            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for component S%d and L%d: (%g, %g)\n", i, j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for component S%d and L%d: (%g, %g)\n", i, j, tmp_eps_cu_MAX, y_MAX);
            }
            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE2  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and comressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = lines[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = lines[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't cause
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced form it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j); //should be negative
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }
            }
        }

        //loop through each fiber group

        //next we do analogues iteration for fiber groups
        for(int j = 0; j < fiberGroups.length(); ++j)
        {
            //foreach jth fiber group we get its y upper coordinate limit
            double y_max_j = fiberGroups[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = fiberGroups[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e. ith component lower edge

            //we calculate direction coefficient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;

            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for component S%d and FG%d: (%g, %g)\n", i, j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for component S%d and FG%d: (%g, %g)\n", i, j, tmp_eps_cu_MAX, y_MAX);
            }

            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE2  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and comressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = fiberGroups[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = fiberGroups[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't cause
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced form it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j); //should be negative
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }
            }

        }
    }

    if(eps_tu == 0) {
        //when ultimate tensile strain are still equal to 0 it means that
        //hasn't been found material wich is bearing tensile stresses
        //we assume than that this ultimate tensile strain equals to 0
        //are allowed in all section and assume eps_tu == 0 for y_tu = y_MIN!
        eps_tu = 0.0; y_tu = y_MIN;
        eps_tu_MIN = eps_tu;

        //now we must search for minimal eps_cu_MAX for y_MAX based on
        //point (eps_tu_i=eps_tu=0, y_min_i=y_tu=y_MIN) and point (eps_cu_j, y_max_j) for each component

        double y_min_i = y_tu;
        double eps_tu_i = eps_tu;

        //loop through each surface
        for(int j=0; j < surfaces.length(); ++j)
        {
            //skip if j-th surface is openning
            if(surfaces[j]->getIsOpening()) continue;

            //foreach jth surface we get its y upper coordinate limit
            double y_max_j = surfaces[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = surfaces[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e ith component lower edge

            //we calculate direction coeffictient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;

            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for pure compression for component S%d: (%g, %g)\n",  j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for pure compression for component  S%d: (%g, %g)\n",  j, tmp_eps_cu_MAX, y_MAX);
            }
            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and comressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = surfaces[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = surfaces[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't cause
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced form it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j); //should be negative
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }
            }
        }

        //loop through each line element

        //next we do analogues iteration for line elements
        for(int j = 0; j < lines.length(); ++j)
        {
            //foreach jth line element we get its y upper coordinate limit
            double y_max_j = lines[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = lines[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e. ith component lower edge

            //we calculate direction coefficient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;

            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for pure compression for component  L%d: (%g, %g)\n", j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for pure compression for component  L%d: (%g, %g)\n", j, tmp_eps_cu_MAX, y_MAX);
            }
            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and comressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = lines[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = lines[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't cause
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced form it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j); //should be negative
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }
            }
        }

        //loop through each fiber group

        //next we do analogues iteration for fiber groups
        for(int j = 0; j < fiberGroups.length(); ++j)
        {
            //foreach jth fiber group we get its y upper coordinate limit
            double y_max_j = fiberGroups[j]->yMaxInCoordinateSystem(centroid, rotationAngle);
            //for this component material we get ultimate compressive strain
            double eps_cu_j = fiberGroups[j]->getMaterial()->compressiveStrainLimit_eps_cu();

            if(y_max_j <= y_min_i)
                continue; //skip when y_max_j i.e. jth component upper edge lays lower than y_min_i i.e. ith component lower edge

            //we calculate direction coefficient (tangens of alfa of slope of line passing through
            //2 points (eps_tu_i, y_min_i) and (eps_cu_j, y_max_j)
            double direction_coefficient = (eps_tu_i - eps_cu_j)/(y_min_i - y_max_j); //should be negative
            //calculating b in eps = direction_coefficient*y + b
            double b = eps_tu_i - direction_coefficient*y_min_i;

            //new temporary values of eps_tu_MIN and eps_cu_MAX for current (eps_cu_j, y_max_j) point and (eps_tu_i, y_min_i) point
            double tmp_eps_tu_MIN = y_MIN*direction_coefficient + b;
            double tmp_eps_cu_MAX = y_MAX*direction_coefficient + b;

            if(DEBUG) {
                fprintf(stderr, "eps_tu_MIN for pure compression for component FG%d: (%g, %g)\n", j, tmp_eps_tu_MIN, y_MIN);
                fprintf(stderr, "eps_cu_MAX for pure compression for component FG%d: (%g, %g)\n", j, tmp_eps_cu_MAX, y_MAX);
            }
            //compressive strains are negative so we assume that new value should be greater than or equal previous one
            //workround for initial case when eps_cu_MAX = 0.0 (compressive strain can't be greater than or equal 0.0)
            //tensile strain are positive so we assume that new value should be less than or equal previous one
            //workround for initial case when eps_tu_MIN = 0.0 (tensile strain can't be less than 0.0)
            if( (tmp_eps_cu_MAX - eps_cu_MAX >= STRAIN_ERROR_TOLERANCE  || eps_cu_MAX == 0.0) && (tmp_eps_tu_MIN <= eps_tu_MIN || eps_tu_MIN == 0.0) )
            {
                //setting helper values for next comparisons of limiting strains
                eps_cu_MAX = tmp_eps_cu_MAX;
                eps_tu_MIN = tmp_eps_tu_MIN;

                //setting new limiting ultimate tensile and compressive strain points
                eps_cu = eps_cu_j; y_cu = y_max_j;
                eps_tu = eps_tu_i; y_tu = y_min_i;

                //next we can calculate corresponding eps_co (pure compression strain limit)
                //if it is material like concrete it has own pure compression strain limits
                //i.e. for eps_cu 0.0035 and for eps_co 0.002
                //if it is steel we assume this pure compression strains that are equal to eps_co = eps_cu
                double eps_co_j = fiberGroups[j]->getMaterial()->pureCompressionStrainLimit_eps_co();
                //next we need to calculate y_co_j and assign it to y_co!
                //for y_max_j we have eps_cu_j for y_min_j we have eps_tu_j (concrete == 0)
                double y_co_j;

                if(eps_co_j == eps_cu_j) {
                    //if we have steel, reinforcement we have workround for such case
                    //we assume that eps_co_j = eps_cu_j and that it has y_co_j = y_max_j
                    y_co_j = y_max_j;
                } else {
                    //we have material with specifically defined pure compression strains
                    //such as concrete and we must calculate y_co_j
                    //formula from PN-EN 1992, p. 77, fig. 6.1: distance eps_co strain from upper edge of component
                    //is equal to (1 - eps_co/eps_cu)*h
                    double y_min_j = fiberGroups[j]->yMinInCoordinateSystem(centroid, rotationAngle);
                    y_co_j = y_max_j - (1 - eps_co_j/eps_cu_j)*(y_max_j - y_min_j);
                }
                //assigning pure compression point
                eps_co = eps_co_j;
                y_co = y_co_j;

                //checking whether such strain limits assumption doesn't casue
                //exceed of eps_co pure compression strain limit by value derived from
                //new eps_tu and eps_cu for y_co
                double eps_co_from_tu_cu = y_co_j*direction_coefficient + b;
                if(eps_co_from_tu_cu < eps_co_j) {
                    //we must limit ultimate compressive strain eps_cu_j and induced from it eps_cu
                    //due to eps_co_j (pure compression strain limit) so as this pure compression strain limit
                    //hasn't been exceeded in any case

                    //WE CALCULATE new_eps_cu_j as intersection of line y_max_j and line passing through points
                    //(eps_tu_i, y_min_i) and (eps_co_j, y_co_j)
                    double direction_coefficient_co = (eps_tu_i - eps_co_j)/(y_min_i - y_co_j);
                    //calculating b_co in eps = direction_coefficient_co*y + b_co
                    double b_co = eps_tu_i - direction_coefficient_co*y_min_i;

                    double new_eps_cu_j = y_max_j*direction_coefficient_co + b_co;
                    //WE CAN ALSO CALCULATE: new_eps_cu_MAX  = y_MAX*direction_coefficient_co + b_co
                    //                       new_eps_tu_MIN = y_MIN*direction_coefficient_co + b_co
                    //AND OVERWRITE CURRENT VALUES IN ORDER TO TAKE THIS AMENDMENT INTO ACCOUNT

                    eps_cu = new_eps_cu_j;
                }
            }

        }
    }

    if(DEBUG) {
        fprintf(stderr, "Finally eps_cu_MAX = %g and eps_tu_MIN = %g.\n", eps_cu_MAX, eps_tu_MIN);
        fprintf(stderr, "FInally strain limit points are: (eps_tu, y_tu) = (%g, %g) and (eps_cu, y_cu) = ( %g, %g)\n", eps_tu, y_tu, eps_cu, y_cu);
        fprintf(stderr, "Finally eps_co = %g for y_co = %g.\n", eps_co, y_co);
    }

    //setting calculated strain limits (eps_cu, y_cu), (eps_tu, y_tu), (eps_co, y_co)
    tensionLimit = eps_tu;
    yCoordinateOfTensionLimit = y_tu;
    compressionLimit = eps_cu;
    yCoordinateOfCompressionLimit = y_cu;
    pureCompressionLimit = eps_co;
    yCoordinateOfPureCompressionLimit = y_co;

}

/**
 * @brief SectionSolver::solveSection
 * This method based on calculated strain limits (eps_tu, y_tu) and (eps_cu, y_cu)
 * and additional limit under pure compression i.e. (eps_co, y_co)
 * Perform in this limits strain profile steerage between this limits.
 * 1. fixed ultimate tensila strain in point (eps_tu, y_tu) and steerage
 *    section upper edge (y_cu) limit strains between eps_tu to eps_cu
 * 2. fixed ultimate compressive strain in point (eps_cu, y_cu)
 *    and steerage section lower limit strains (y_tu) from eps_tu to eps_co_on_y_tu
 *    where eps_co_on_y_tu is defined as intersection with y_tu line passing through
 *    points (eps_cu, y_cu) and (eps_co, y_co).
 * 3. Steereage strain profile using profile curvature angle fi from
 *    last strain profile fi to 0. (straight line on eps_co)
 *     ATTANTION! positive fi angle denotes compression in upper section tendons
 */
void SectionSolver::solveSection(void)
{

       //strain profile steerage limit points A, B, C definition
       //where A = (eps_tu, y_tu), B = (eps_cu, y_cu), C = (eps_co, y_co)
       double eps_cu = 0.0, y_cu = 0.0;
       double eps_tu = 0.0, y_tu = 0.0;
       double eps_co = 0.0, y_co = 0.0;
       double eps_co_on_y_tu = 0.0;

       //helper variables to store eps_cu_MAX and eps_tu_MIN
       //which are current strain profile strain valueas as defined
       //for y_MAX and y_MIN section edges
       //double eps_cu_MAX = 0.0, eps_tu_MIN = 0.0;

   //LOOP LOCAL COORDINATE SYSTEM (AND STRAIN PROFILE DIRECTION) AROUND THE SECTION
   for(int i=0; i<360; i++) {

       //set new rotationAngle value in radians as equal to: i*PI/180
       this->setRotationAngle( (i*M_PI/180) );

       //section edges in current rotated local coordinate system
       //in which we are searchin strain profiles definitions
       Point *centroid = this->getGeometricalCenter();
       double y_MIN = this->yMinInCoordinateSystem(centroid, rotationAngle);
       double y_MAX = this->yMaxInCoordinateSystem(centroid, rotationAngle);

       //finding (eps_tu, y_tu) point
       findTensileStrainLimit(&eps_tu, &y_tu);
       //finding (eps_cu, y_cu) point
       findCompressiveStrainLimit(&eps_cu, &y_cu);
       //finding (eps_co, y_co) point (pure compression strain limit)
       findPureCompressionStrainLimit(&eps_co, &y_co);

       //finding (eps_co_on_y_tu, y_tu) point which lies on intersection of line
       //passing through (eps_cu, y_cu) and (eps_co, y_co) and line y_tu
       if(eps_cu == eps_co && y_cu == y_co) {
           //we have choosen material which limits compressive strains
           //such that hasn't got defined pure compression strain limit
           //and than this limit has been assumed as equal to eps_co = eps_cu
           eps_co_on_y_tu = eps_co; //we can steerage in 2) step on y_tu from eps_tu to eps_co_on_y_tu = eps_co = eps_cu !!!
       } else {
           //else we are calculating eps_co_on_y_tu as intersection of line y_tu
           //and line passing through points (eps_cu, y_cu) and (eps_co, y_co)
            double direction_coefficient = (eps_cu - eps_co)/(y_cu - y_co);
            //calculating b in eps = y*direction_coefficient + b
            double b = eps_cu - y_cu*direction_coefficient;
            //at last we  get eps_co_on_y_tu limit of strain profile iteration on y_tu
            eps_co_on_y_tu = y_tu*direction_coefficient + b;
       }

        //NEXT WE GENERATE STRAIN PROFILE OBJECTS FOR CURRENT LOCAL COORDINATE SYSTEM
        //ROTATED IN RELATION TO CENTRAL GLOBAL COORDINATE SYSTEM BY ANGLE rotationAngle
        //WHICH SHOULD BE INTEGRATED IN ORDER TO ACHIVE INTERNAL ACTIONS N, Mx, My

        //1. loop for fixed (eps_tu, y_tu) and steerage of y_cu strain limit between [eps_tu, eps_cu]
       for(double eps_y_cu = eps_tu; eps_y_cu > eps_cu; eps_y_cu -= Y_CU_STEP)
       {
            //we have strain profile defined by 2 points
            //(eps_tu, y_tu) and (eps_y_cu, y_cu)
            //for local coordinate system (centroid, rotationAngle)
            integrateStrainProfileFrom(eps_tu, y_tu, eps_y_cu, y_cu);
       }

        //2. loop for fixed (eps_cu, y_cu) and steerage of y_tu strain limit between [eps_tu, eps_co_on_y_tu]
       for(double eps_y_tu = eps_tu; eps_y_tu > eps_co_on_y_tu; eps_y_tu -= Y_TU_STEP)
       {
           //we have strain profile defined by 2 points
           //(eps_y_tu, y_tu) and (eps_cu, y_cu)
           //for local coordinate system (centroid, rotationAngle)

           integrateStrainProfileFrom(eps_y_tu, y_tu, eps_cu, y_cu);
       }

       //3. having defined 2 points (eps_co_on_y_tu, y_tu) and (eps_cu, y_cu)
       //   we get curvature fi of this strain profile as fi = atan[(eps_co_on_y_tu - eps_cu)/(y_cu - y_tu)]
       //   next we steerage this fi curvature from fi to 0
       //   and are getting eps_y_cu and eps_y_tu values as defined below:
       //   eps_y_cu = eps_co - (y_cu - y_co)*tan(fi)
       //   eps_y_tu = eps_co + (y_co - y_tu)*tan(fi)

       if(DEBUG)
        fprintf(stderr, "(eps_co, y_co) = (%g, %g), eps_co_on_y_tu: %g\n", eps_co, y_co, eps_co_on_y_tu);

       double fi = (eps_co_on_y_tu - eps_cu)/(y_cu - y_tu);
       for(double fi_itr = fi; fi_itr >= 0; fi_itr -= FI_STEP)
       {
           double eps_y_cu = eps_co - (y_cu - y_co)*fi_itr;
           double eps_y_tu = eps_co + (y_co - y_tu)*fi_itr;

           if(DEBUG)
            fprintf(stderr, "************** STRAIN PROFILE eps_co ROTATAION *************\n");

           //we have strain profile definde by 2 points
           //(eps_y_tu, y_tu) and (eps_y_cu, y_cu)
           //for local coordinate system (centroid, rotationAngle)

           integrateStrainProfileFrom(eps_y_tu, y_tu, eps_y_cu, y_cu);
       }
   }
}


void SectionSolver::integrateStrainProfileFrom(double eps_y_tu, double y_tu,
                                               double eps_y_cu, double y_cu)

{
    //if(DEBUG)
        fprintf(stderr, "Strain Profile: omega= %g, eps_y_tu = %g, y_tu=%g, eps_y_cu=%g, y_cu=%g\n",
            rotationAngle, eps_y_tu, y_tu, eps_y_cu, y_cu);

    Point *centroid = this->getGeometricalCenter();
    double y_MIN = this->yMinInCoordinateSystem(centroid, rotationAngle);
    double y_MAX = this->yMaxInCoordinateSystem(centroid, rotationAngle);

    //fprintf(stderr, "y_MIN = %g, y_MAX = %g\n", y_MIN, y_MAX);

    //defining linear function of strain profile passing through above 2 points
    double direction_coefficient = (eps_y_tu - eps_y_cu)/(y_tu - y_cu);
    //in function eps = y*direction_coefficient + b we determine b
    double b = eps_y_tu - y_tu*direction_coefficient;

    //calculate strain limits for current profile on y_MIN and y_MAX section edges
    double eps_y_MIN = y_MIN*direction_coefficient + b; //for y_MIN
    double eps_y_MAX = y_MAX*direction_coefficient + b; //for y_MAX

    //CONSTRUCTING STRAIN PROFILE OBJECT FOR CURRENT STRAIN LIMITS!
    StrainProfile *profile = new StrainProfile(this, rotationAngle, eps_y_MIN, eps_y_MAX);

    //having StrainProfile we can integrate it and get corresponding internal actions N, Mx, My
    double N = this->summateInternalActions_N(profile); //[N]
    double localeMx = this->summateInternalActions_Mx(profile);
    double localeMy = this->summateInternalActions_My(profile);

    delete profile;

    double Mx = transferBackFromLocalMx(localeMx, localeMy, rotationAngle); //[Nm]
    double My = transferBackFromLocalMy(localeMx, localeMy, rotationAngle); //[Nm]

    //printing resultant point
    if(DEBUG)
    fprintf(stderr, "{omega = %g, eps_y_MIN = %g, y_MIN = %g, eps_y_MAX = %g, y_MAX = %g} ==> {N = %g, Mx = %g, My = %g}\n",
            rotationAngle, eps_y_MIN, y_MIN, eps_y_MAX, y_MAX, N, Mx, My);

    //creating Vertex in 3D Mx-My-N coordinate system and pushing it to QList
    Vertex3D *vertex = new Vertex3D(N/1000,Mx/1000000,My/1000000); //[kN, kNm, kNm]

    /// TEST ONLY - FOR US UNITS
    //Vertex3D *vertex = new Vertex3D(N,Mx/10000,My/10000);
    ///TEST ONLY
    vertices.push_back(vertex);
}

/**
 * @brief SectionSolver::getVertices
 * @return - list of Vertex3D objects which are points (Mx, My, N) in 3D coordinate
 *  system which specifies interaction surface and interaction curves
 */
QList<Vertex3D *> &SectionSolver::getVertices(void)
{
    //if there hasn't been calculated vertices of interaction surface solid
    //then we need to calculate them previously before we return list of vertices
    //if(vertices.length() == 0)
        //solveSection(); - omitting calculation user can load vertices from file!

    //returns const reference to QList with pointers to Vertex3D objects (Mx, My, N)
    return vertices;
}

/**
 * @brief SectionSolver::sortVerticesByN
 * this method is sorting vertices Vertex3D in QList
 * using comparator Vertex3DComparatorByN()
 * while sorting vertices are ordered firstly by increasing N value
 * and then secondly by increasing angle arctan(My/Mx) (from 0 degree to 360 degree)
 * in order to have ordered vertices for drawing interaction curves for given value of N (axial force)
 */
void SectionSolver::sortVerticesByN(void)
{
  //we are using qSort() function from QAlgorithms library
  //passing interator pointing to QList begining and ending
  //and functional object Vertex3DComparatorByN() based on
  //which sorting process is executed
  qSort(vertices.begin(), vertices.end(), Vertex3DComparatorByN() );
}

/**
 * @brief SectionSolver::sortVerticesByMx
 * this method is sorting vertices Vertex3D in QList
 * using comparator Vertex3DComparatorByMx()
 * while sorting vertices are ordered firstly by increasing Mx value
 * and then secondly by increasing angle arctan(N/My) (from 0 degree to 360 degree)
 * in order to have ordered vertices for drawing interaction curves for given value of Mx (bending moment in relation to x-axis)
 */
void SectionSolver::sortVerticesByMx(void)
{
    qSort(vertices.begin(), vertices.end(), Vertex3DComparatorByMx());
}

/**
 * @brief SectionSolver::sortVerticesByMy
 * this method is sorting vertices Vertex3D in QList
 * using comparator Vertex3DComparatorByMy()
 * while sorting vertices are ordered firstly by increasing My value
 * and then secondly by increasing angle arctan(N/Mx) (from 0 degree to 360 degree)
 * in order to have ordered vertices for drawing interaction curves for given value of My (bending moment in relation to y-axis)
 */
void SectionSolver::sortVerticesByMy(void)
{
    qSort(vertices.begin(), vertices.end(), Vertex3DComparatorByMy());
}

void SectionSolver::sortVerticesByAlfa(void)
{
    qSort(vertices.begin(), vertices.end(), Vertex3DComparatorByAlfa());
}

