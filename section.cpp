#include "section.h"
#include "flags.h"
#include "combinationgenerator.h"
#include <cmath>

Section::Section()
{
}

Surface** Section::getSurfaces() {
    return surfaces.data();
}

Line** Section::getLines() {
    return lines.data();
}

FiberGroup** Section::getFiberGroups() {
    return fiberGroups.data();
}

int Section::getNumberOfSurfaces() {
    return surfaces.length();
}

int Section::getNumberOfLines() {
    return lines.length();
}

int Section::getNumberOfFiberGroups() {
    return fiberGroups.length();
}

void Section::addSurface(Surface *surface)
{
    surfaces.append(surface);
}

void Section::addLine(Line *line)
{
    lines.append(line);
}

void Section::addFiberGroup(FiberGroup *fiberGroup)
{
    fiberGroups.append(fiberGroup);
}

/**
 * @brief Section::weightedArea
 * @return
 * Method calculate weighted area based on area of all surfaces in the section
 * subtructing all intersections of area (sometimes also adding if too much
 * has been subtructed). This surfeces are subtructed using combinatorics formula
 * in order to find all surfeces intersection. We use special generator of
 * number of combinations C_n_over_k = n!/(k!*(n-k)!)
 * While calculating weights of each surface, line, fiber group elements
 * reference module of elasticity is E of first surface S[0]
 * 1. add up all surfaces S_i multiplied by corresponding weight w_i
 * 2. if surface hasn't been adaptive linearized perform adaptive linearization
 * 3. for each generated combination we get numbers of surfaces which should be intersected together
 * 4. sign of intersection contribiution to all weighted area is determined by even or odd number of intersected areas
 * 5. we log this intersection and printing to stderr in order to check correctness of calculations
 * 6. intersected subracted or added area  is multiplied by corresponding weight
 * 7. next we take into account line elements and fiber groups
 * 8. adding to total area areas of all line elements multiplied by its weight
 * 9. subtracting material that lays underneath line element
 * !!! WE ASSUME THAT LINE ELEMENT CAN BE EMBEDED ONLY IN ONE SURFACE (IT SIMPLIFIES CALCULATIONS)
 * 10. adding to total area areas of all fiber groups multiplied by its weight
 * 11. subtracting material that lays underneath fiber groups
 * !!! WE ASSUME THAT FIBER GROUPS CAN BE EMBEDDED ONLY IN ONE SURFACE (IT SIMPLIFIES CALCULATIONS)
 */
double Section::weightedArea(void)
{
    double weightedArea = 0; //adding each surface area


    CombinationGenerator generator(surfaces.length());
    std::vector< std::vector<int> > combinations = generator.getAllCombinations();

    if(surfaces.length() < 1) return weightedArea;

    double E = surfaces[0]->getMaterial()->getModuleOfElasticity_E();

    //walk through each surface from background to foregound matrial
    for(int i=0; i< surfaces.length(); i++) {

        //modul of elasticity of ith area
        double E_i = surfaces[i]->getMaterial()->getModuleOfElasticity_E();
        //weight of ith area
        double n_i = E_i/E;
        //we are getting ith surface area
        if(!surfaces[i]->hasBeenAdaptiveLinearized())
            //if surface i hasn't been adaptive linearized
            surfaces[i]->performAdaptiveLinearization();
        //calculating ith surface area
        double A_i = surfaces[i]->calculatePolygonArea();

        //simple adding all surfaces areas multiplied by its weight
        weightedArea += n_i * A_i;

    }

    //now we need to make more sophisticated calculation
    //too subtract contribution from overlapping surfaces
     char tmpString[1000];
     std::string logString("");
     double intersectionArea = 0.0;

     for(auto itr = combinations.begin();
         itr < combinations.end(); itr++) {
            std::vector<int> row = *itr;
            //we determine the sign
            double sign = std::pow( -1.0, (int)(row.size()-1));
            //index necessary to calculate n_i
            int weightIdx = row[1] -1;
            sprintf(tmpString, "%dxn%d", (int) sign, weightIdx);
            logString += tmpString;
            //we calculate n_i
            //modul of elasticity of ith area
            double E_i = surfaces[weightIdx]->getMaterial()->getModuleOfElasticity_E();
            //weight of ith area
            double n_i = E_i/E;

            //we are retrieving foreground area
            //for which we will be checking intersection
            //with areas laying underneath it
            Surface *firstSurface = surfaces[row[0]-1];
            sprintf(tmpString, "xA%d", row[0]-1);
            logString += tmpString;

            intersectionArea = calculateIntersectionArea(firstSurface, ++row.begin(), row.end(), logString);
            /*
            //we start loop from second element
            for(auto elem = ++(row.begin());
                elem < row.end(); elem++ ) {
                //we make intersection in accordance with
                //generated combinations
                int areaIdx = *elem -1;
                sprintf(tmpString, "xA%d", areaIdx);
                logString += tmpString;
                QVarLengthArray<Surface*> *intersections =  area->intersection(surfaces[areaIdx]);
                //area = intersections->at(0);
            }*/

            sprintf(tmpString, " (intersection n_i= %g, A = %g)\n", n_i, intersectionArea);
            logString += tmpString;

            weightedArea += sign *n_i * intersectionArea;

     }

        fprintf(stderr, "%s", logString.c_str());
        fprintf(stderr, "Calculated weigth area (only surfaces): %g\n", weightedArea);

     //next we need to take into account line elements contribiution

     //adding to total area weightedArea areas of all line elements multiplied by its weights
     for(int i=0; i < lines.length(); i++)
     {
        //module of elasticity of area of ith line element
         double E_i = lines[i]->getMaterial()->getModuleOfElasticity_E();
        //weight of area of ith line element
         double n_i = E_i/E;
        //total area of ith line element
         double A_i = lines[i]->getTotalArea();

         //simple adding all line elements areas multiplied by its weight
         weightedArea += n_i * A_i;

         //subtructing for each line element contribiution from surface laying
         //underneath it. !!!Line element can be embeded only into single surface
         //moreover we assume that line elements are added on top of all surfaces
         //and therefor thay can not be embeded in surface fragment overlapped by
         //other more foreground surface!!!

         //in the end we can check in which surface (going from top to bottom)
         //line element is embedded by only searching surface which contains
         //line element origin!
         for(int j=surfaces.length()-1; j >= 0; --j)
         {
             //we are going from foreground surface (top)
             //to background surface (bottom)
             Point *lineOrigin = lines[i]->getOrigin();
             if(surfaces[j]->pointIsInsieThePolygon(lineOrigin))
             {
                 //module of elasticity of surface overlapped by line element
                 double E_j = surfaces[j]->getMaterial()->getModuleOfElasticity_E();
                 //weight of area of this surface
                 double n_j = E_j/E;
                 //total overlapped area assumed to be equal to line element total area
                 double A_j = A_i;

                 //simple subtracting surface fragment which is overlapped by line element
                 weightedArea -= n_j*A_j;

                 break;
             }
         }
    }


     //next we need to take into account in similar manner fiber group contribiution

     //adding to total area weightedArea areas of all fiber groups multiplied by its weights
     for(int i=0; i< fiberGroups.length(); i++)
     {
         //module of elasticity of ith fiber group
         double E_i = fiberGroups[i]->getMaterial()->getModuleOfElasticity_E();
         //weight of area of ith fiber group
         double n_i = E_i/E;
         //total area of ith fiber group element
         double A_i = fiberGroups[i]->getTotalArea();

         //simple adding all fiber group areas multiplied by its weight
         weightedArea += n_i*A_i;

         //subtructing for each fiber group contribiution from surface laying
         //underneath it. !!!Fiber group can be embeded only into single surface
         //moreover we assume that fiber groups are added on top of all surfaces
         //and therefore thay can not be embeded in surface fragment overlapped by
         //other more foreground surface!!!

         //in the end we can check in which surface (going from top to bottom)
         //whether fiber groupe is embedded by only searching surface which contains
         //fiber group origin!
         for(int j=surfaces.length()-1; j >= 0; --j)
         {
             //we are going from foreground surface (top)
             //to background surface (bottom)
             Point *fiberGroupOrigin = fiberGroups[i]->getOrigin();
             if(surfaces[j]->pointIsInsieThePolygon(fiberGroupOrigin))
             {
                 //module of elasticity of surface overlapped by fiber group element
                 double E_j = surfaces[j]->getMaterial()->getModuleOfElasticity_E();
                 //weight of area of this surface
                 double n_j = E_j/E;
                 //total overlapped area assumed to be equal to fiber group total area
                 double A_j = A_i;

                 //simple subtracting surface fragment which is overlapped by fieber group
                 weightedArea -= n_j*A_j;

                 break;
             }
         }
     }

     fprintf(stderr, "Calculated weight area (total): %g\n", weightedArea);

     return weightedArea;
}

/**
 * @brief Section::calculateIntersectionArea
 * @param surface - surface which will be intersected with remaining
 *                  surfaces indicated by indices in vector<int> list
 *                  starting from beginItr to endItr
 * @param beginItr - indicates index of first surface to intersect with surface (first argument)
 * @param endItr - last surface to intersect with surface (first argument)
 * @return - cumulative area of surfeces intersection
 *This method is called recursively in order to calculate area (double) of intersection
 * surface (passed as first argument) with all other surfaces indicated by vector<int> of indices
 * 1. intersect surface with first surface on the list (vector<int>)
 * 2. foreach intersection (potentially we can get >=0 surfaces)  do:
 *    3. recursively intersect this new surface with remaining surfaces on the list (vector<int>)
 *       calculateIntersectionArea(intersection[i], ++beginItr, endItr);
 *    4. cumulate area (double) returned from recursively called function into double intersectionArea
 *    5. return this cumulative area
 * 6. Ending condition is when beginItr == endItr then this function return just area of surface
 */
double Section::calculateIntersectionArea(Surface *surface,
                                          std::vector<int>::const_iterator beginItr,
                                          std::vector<int>::const_iterator endItr,
                                          std::string& logString) //logString only for testing
{
        //if there isn't any remaining surface to intesect with surface passed as first argument
        if(beginItr == endItr)
            //then just return area of surface
            return surface->calculatePolygonArea();
        //define variable cumulating intersections areas
        double intersectionArea = 0.0;

        //do intersection of surface (first parameter) with first surface on vector<int> list
        int areaIdx = *beginItr -1; //adjustment from {1-n} to {0 - (n-1)} range
        QVarLengthArray<Surface *> *intersections = surface->intersection(surfaces[areaIdx]);
        //logging intersection
        char tmpString[10];
        sprintf(tmpString, "xA%d", areaIdx);
        logString += tmpString;

        //for each intersection surfaces (we can get >=0 surfaces) recursively call
        //this function to intersect it with remaining surfaces in the list vector<int>
        //we need to increment beginItr
        //return areas are cumulated in intersectionArea variable
        ++beginItr;
        for(int i=0; i < intersections->length(); i++) {
            Surface *intersection = intersections->at(i);
            intersectionArea += calculateIntersectionArea(intersection, beginItr, endItr, logString);

        }

        return intersectionArea;
}

/**
 * @brief Section::weightedFirstMomentOfArea_Sx
 * @return
 */
double Section::weightedFirstMomentOfArea_Sx(void)
{
    //we are assuming default coordinate system
    Point *referencePoint = new Point(0,0);
    return weightedFirstMomentOfArea_Sx(referencePoint);

}

/**
 * @brief Section::weightedFirstMomentOfArea_Sx
 * @param p - reference point in relation to the first moment of area
 *            is calculated
 * @return
 * Function calculates weighted first moment of area Sx
 * in relation to Point p.
 * It takes into account contribiution form each surfaces, subtructing/addding contribiution
 * from intersection of surfaces. Additionaly it takes into account contribiution form
 * line elements and fiber groups also subtructing contribiution form overlapped surface fragments
 * Note that line element or fiber group can be embedded only in one surface material so
 * to find which overlapped surface fragment contribiution should be subtracted we go in loop
 * form the most top surface to the bottom surface and check wether line origin or fiber group origin
 * belongs to this surface.
 * 1. add up first moments of area Sx_i of all surfaces S_i multiplied by corresponding weight coefficients w_i
 * 2. if surface hasn't been adaptive linearized perform adaptive linearization
 * 3. for each generate combination we get numbers of surfaces which should be intersected together
 * 4. sign of intersection contribiution to total weighted first moment of area Sx is determined by even or odd
 *    number of intersected areas
 * 5. we log this intersections and print them to stderr in order to check correctness of calculation
 * 6. for each interected area we calculate its first moment of area Sx and subtract or add it depending on sign and multiply by weight coefficient
 * 7. next we take into account line elements and fiber groups
 * 8. adding to total weighted Sx contribiution form line elements (Sx_i*w_i) and from fiber groups (Sx_i*w_i)
 * 9. subtructing contribiutions from material (surface fragment) that lays underneath line elements or fiber groups
 * !!! WE ASSUME THAT LINE ELEMENT OR FIBER GROUP CAN BE EMBEDDED ONLY IN ONE SURFACE
 * 10. we are passing through each surface form top to bottom and check whether it contains origin of line element or fiber group
 *    and subtract Sx contribiution from this fragment as (- n_overlapped_surface * Sx_line_element)
 *    or (- n_overlapped_surface * Sx_fiber_group)
 */
double Section::weightedFirstMomentOfArea_Sx(Point *p)
{
    double weightedSx = 0; //cumulate Sx values

    CombinationGenerator generator(surfaces.length());
    std::vector< std::vector<int> > combinations = generator.getAllCombinations();

    if(surfaces.length() < 1) return weightedSx;

    double E = surfaces[0]->getMaterial()->getModuleOfElasticity_E();

    //walk through each surface from background to foreground material
    for(int i=0; i< surfaces.length(); i++) {

        //modul of elasticity of ith area
        double E_i = surfaces[i]->getMaterial()->getModuleOfElasticity_E();
        //wieght of ith area
        double n_i = E_i/E;
        //we are getting ith surface area
        if(!surfaces[i]->hasBeenAdaptiveLinearized())
            //if surface i hasn't been adaptive linearized
            surfaces[i]->performAdaptiveLinearization();
        //calculate ith surface first moment of are
        double Sx_i = surfaces[i]->firstMomentOfArea_Sx(p);

        //simple adding all surfaces first moment of area Sx multiplied by its weight
        weightedSx += n_i *Sx_i;
    }

    //now we need to make more sophisticated calculation
    //in order to subtract contribution from overlapping surfaces
    char tmpString[1000];
    std::string logString("");
    double intersectionSx = 0.0;

    for(auto itr = combinations.begin();
        itr < combinations.end(); itr++) {

        std::vector<int> row = *itr;
        //we determine the sign
        double sign = std::pow(-1.0, (int) (row.size()-1));
        //index necessery to calculate n_i
        int weightIdx = row[1]-1;
        sprintf(tmpString, "%dxn%d", (int) sign, weightIdx);
        logString += tmpString;
        //we calculate n_i
        //modul of elasticity of ith surface
        double E_i = surfaces[weightIdx]->getMaterial()->getModuleOfElasticity_E();
        //weight of ith surface
        double n_i = E_i/E;

        //we are retrieving foreground surface
        //for which we will be checking intersections
        //with surfaces laying underneath it
        Surface *firstSurface = surfaces[row[0]-1];
        sprintf(tmpString, "xSx%d", row[0]-1);
        logString +=tmpString;

        intersectionSx = calculateIntersectionSx(firstSurface,
                                                 ++row.begin(),
                                                 row.end(),
                                                 p, //origin of coordinate system
                                                 logString);

        sprintf(tmpString, " (intersection n_i = %g, Sx = %g)\n", n_i, intersectionSx);
        logString += tmpString;

        weightedSx += sign*n_i*intersectionSx;
    }

    fprintf(stderr, "%s", logString.c_str());

    fprintf(stderr, "Calculated weight Sx (only surfaces): %g\n", weightedSx);

    //next we need to take into account line elements

    //firstly we add contribiution from each line element as n_i *Sx_i
    //where n_i = E_i/Eref ; E_i is module of elasticity of ith line element
    //Sx_i is first moment of area in relation to x axis of ith line element
    for(int i=0; i < lines.length(); i++)
    {
        //module of elasticity of ith line element
        double E_i = lines[i]->getMaterial()->getModuleOfElasticity_E();
        //weight of Sx_i of ith line element
        double n_i = E_i/E; //where E i reference module of elasticity
        //first moment of area Sx_i of ith line element
        double Sx_i = lines[i]->firstMomentOfArea_Sx(p);

        //simple adding all line elements first moment of area Sx_i contribiutions multiplied by its weight
        weightedSx += n_i*Sx_i;

        //next we need to subtruct contribiution to weighted Sx from surface fragment overlapped by line element
        //in order to do that we pass through each surface from top to bottom and check whether it contains origin of line elemnt
        //if we find such surface we subtract its contribiution and break out from loop
        //IMPORTANT! WE ASSUME THAT LINE ELEMENT CAN BE EMBEDDED ONLU IN ONE SURFACE Si

        for(int j= surfaces.length()-1;  j>=0; --j)
        {
            Point *lineOrigin = lines[i]->getOrigin();
            if(surfaces[j]->pointIsInsieThePolygon(lineOrigin))
            {
                //module of elasticity of surface overlapped by line element
                double E_j = surfaces[j]->getMaterial()->getModuleOfElasticity_E();
                //weight of Sx_j of this surface
                double n_j = E_j/E;
                //first moment of area of this overlapped surface is equal to Sx_i (of ith line element)
                double Sx_j = Sx_i;

                //simple subtructing surface fragment contribiution to Sx which is overlapped by ith line element
                weightedSx -= n_j*Sx_j;
                break;
            }
        }
    }

    //next we need to take into account fiber groups
    //firstly we add contribiution from each fiber group as n_i *Sx_i
    //where n_i = E_i/Eref ; E_i is module of elasticity of ith fiber group
    //Sx_i is first moment of area in relation to x axis of ith fiber group
    for(int i=0; i < fiberGroups.length(); i++)
    {
        //module of elasticity of ith fiber group
        double E_i = fiberGroups[i]->getMaterial()->getModuleOfElasticity_E();
        //weight of Sx_i of ith fiber group
        double n_i = E_i/E; //where E i reference module of elasticity
        //first moment of area Sx_i of ith fiber group
        double Sx_i = fiberGroups[i]->firstMomentOfArea_Sx(p);

        //simple adding all fiber groups first moment of area Sx_i contribiutions multiplied by its weight
        weightedSx += n_i*Sx_i;

        //next we need to subtruct contribiution to weighted Sx from surface fragment overlapped by fiber group
        //in order to do that we pass through each surface from top to bottom and check whether it contains origin of fiber group
        //if we find such surface we subtract its contribiution and break out from loop
        //IMPORTANT! WE ASSUME THAT FIBER GROUP CAN BE EMBEDDED ONLY IN ONE SURFACE Si
        for(int j= surfaces.length()-1;  j>=0; --j)
        {
            Point *fiberGroupOrigin = fiberGroups[i]->getOrigin();
            if(surfaces[j]->pointIsInsieThePolygon(fiberGroupOrigin))
            {
                //module of elasticity of surface overlapped by ith fiber group
                double E_j = surfaces[j]->getMaterial()->getModuleOfElasticity_E();
                //weight of Sx_j of this surface
                double n_j = E_j/E;
                //first moment of area of this overlapped surface is equal to Sx_i (of ith fiber group)
                double Sx_j = Sx_i;

                //simple subtructing surface fragment contribiution to Sx which is overlapped by ith line element
                weightedSx -= n_j*Sx_j;

                break;
            }
        }
    }

    fprintf(stderr, "Calculated weight Sx (total): %g\n", weightedSx);

    return weightedSx;

}


/**
 * @brief Section::calculateIntersectionSx
 * @param surface - surface which will be intersected with remaining
 *                  surfaces indicated by indices in vector<int> list
 *                  starting from beginItr to endItr
 * @param beginItr - indicates index of first surface to intersect with surface (first argument)
 * @param endItr - last surface to intersect with surface (first argument)
 * @param referencePoint = point in relation to we calculate Sx
 * @param logString - only testing intersection combinations list
 * @return - cumulative Sx of surface intersection
 * This method is called recursively in order to calculate Sx (first moment of area)
 * of intersection surface (passed as first arg) with all other surfaces indicated
 * by vector<int>
 *1. intersect surface with first surface on the list (vector<int>)
 *2. foreach intersection (potentially we can get >=0 surfaces)  do:
 *    3. recursively intersect this new surface with remaining surfaces on the list (vector<int>)
 *       calculateIntersectionSx(intersection[i], ++beginItr, endItr);
 *    4. cumulate Sx (double - first moment of area) returned from recursively called function into double intersectionSx
 *    5. return this cumulative Sx
 * 6. Ending condition is when beginItr == endItr then this function return just Sx of surface
 */
double Section::calculateIntersectionSx(Surface *surface,
                               std::vector<int>::const_iterator beginItr,
                               std::vector<int>::const_iterator endItr,
                               Point *referencePoint,
                               std::string& logString)
{

    //if there isn't any remaining surface to intersect with surface passed as first argument
    if(beginItr == endItr)
        //then just return Sx of surface
        return surface->firstMomentOfArea_Sx(referencePoint);

    //define variable cumulating intersections Sx
    double intersectionSx = 0.0;

    //do intersection of surfaces (first parameter) with the first surface on vector<int> list
    int idx = *beginItr -1;//adjustment from {1-n} to {0-(n-1)} range

    QVarLengthArray<Surface *> *intersections = surface->intersection(surfaces[idx]);
    //logging intersection
    char tmpString[10];
    sprintf(tmpString, "xSx%d", idx);
    logString += tmpString;

    //foreach intersection surfaces (we can get >=0 surfaces) recursively call
    //this function to intersect it with remaining surfaces in the list vector<int>
    //we need to increment beginItr
    //return Sx values are cumulated in intersectionSx variable
    ++beginItr;

    for(int i=0; i < intersections->length(); i++) {
        Surface *intersection = intersections->at(i);
        intersectionSx += calculateIntersectionSx(intersection, beginItr, endItr, referencePoint, logString);

    }

    return intersectionSx;
}


/**
 * @brief Section::weightedFirstMomentOfArea_Sy
 * @return
 */
double Section::weightedFirstMomentOfArea_Sy(void)
{
    //we are assuming default coordinate system
    Point *referencePoint = new Point(0,0);
    return weightedFirstMomentOfArea_Sy(referencePoint);

}

/**
 * @brief Section::weightedFirstMomentOfArea_Sy
 * @param p - reference point in relation to the first moment of area
 *            is calculated
 * @return
 * Function calculates weighted first moment of area Sy
 * in relation to Point p.
 * It takes into account contribiution from each surfaces, subtructing/addding contribiution
 * from intersection of surfaces. Additionaly it takes into account contribiution form
 * line elements and fiber groups also subtructing contribiution form overlapped surface fragments
 * Note that line element or fiber group can be embedded only in one surface material so
 * to find which overlapped surface fragment contribiution should be subtracted we go in loop
 * form the most top surface to the bottom surface and check wether line origin or fiber group origin
 * belongs to this surface.
 * 1. add up first moments of area Sy_i of all surfaces S_i multiplied by corresponding weight coefficients w_i
 * 2. if surface hasn't been adaptive linearized perform adaptive linearization
 * 3. for each generate combination we get numbers of surfaces which should be intersected together
 * 4. sign of intersection contribiution to total weighted first moment of area Sy is determined by even or odd
 *    number of intersected areas
 * 5. we log this intersections and print them to stderr in order to check correctness of calculation
 * 6. for each interected area we calculate its first moment of area Sy and subtract or add it depending on sign and multiply by weight coefficient
 * 7. next we take into account line elements and fiber groups
 * 8. adding to total weighted Sy contribiution form line elements (Sx=y_i*w_i) and from fiber groups (Sy_i*w_i)
 * 9. subtructing contribiutions from material (surface fragment) that lays underneath line elements or fiber groups
 * !!! WE ASSUME THAT LINE ELEMENT OR FIBER GROUP CAN BE EMBEDDED ONLY IN ONE SURFACE
 * 10. we are passing through each surface form top to bottom and check whether it contains origin of line element or fiber group
 *    and subtract Sy contribiution from this fragment as (- n_overlapped_surface * Sy_line_element)
 *    or (- n_overlapped_surface * Sy_fiber_group)
 */
double Section::weightedFirstMomentOfArea_Sy(Point *p)
{
    double weightedSy = 0; //cumulate Sy values

    CombinationGenerator generator(surfaces.length());
    std::vector< std::vector<int> > combinations = generator.getAllCombinations();

    if(surfaces.length() < 1) return weightedSy;

    double E = surfaces[0]->getMaterial()->getModuleOfElasticity_E();

    //walk through each surface from background to foreground material
    for(int i=0; i< surfaces.length(); i++) {

        //modul of elasticity of ith area
        double E_i = surfaces[i]->getMaterial()->getModuleOfElasticity_E();
        //wieght of ith area
        double n_i = E_i/E;
        //we are getting ith surface area
        if(!surfaces[i]->hasBeenAdaptiveLinearized())
            //if surface i hasn't been adaptive linearized
            surfaces[i]->performAdaptiveLinearization();
        //calculate ith surface first moment of are
        double Sy_i = surfaces[i]->firstMomentOfArea_Sy(p);

        //simple adding all surfaces first moment of area Sy multiplied by its weight
        weightedSy += n_i *Sy_i;
    }

    //now we need to make more sophisticated calculation
    //in order to subtract contribution from overlapping surfaces
    char tmpString[1000];
    std::string logString("");
    double intersectionSy = 0.0;

    for(auto itr = combinations.begin();
        itr < combinations.end(); itr++) {

        std::vector<int> row = *itr;
        //we determine the sign
        double sign = std::pow(-1.0, (int) (row.size()-1));
        //index necessery to calculate n_i
        int weightIdx = row[1]-1;
        sprintf(tmpString, "%dxn%d", (int) sign, weightIdx);
        logString += tmpString;
        //we calculate n_i
        //modul of elasticity of ith surface
        double E_i = surfaces[weightIdx]->getMaterial()->getModuleOfElasticity_E();
        //weight of ith surface
        double n_i = E_i/E;

        //we are retrieving foreground surface
        //for which we will be checking intersections
        //with surfaces laying underneath it
        Surface *firstSurface = surfaces[row[0]-1];
        sprintf(tmpString, "xSx%d", row[0]-1);
        logString +=tmpString;

        intersectionSy = calculateIntersectionSy(firstSurface,
                                                 ++row.begin(),
                                                 row.end(),
                                                 p, //origin of coordinate system
                                                 logString);

        sprintf(tmpString, " (intersection n_i = %g, Sy = %g\n",  n_i, intersectionSy);
        logString += tmpString;

        weightedSy += sign*n_i*intersectionSy;
    }

    fprintf(stderr, "%s", logString.c_str());

    fprintf(stderr, "Calculated weight Sy (only surfaces): %g\n", weightedSy);

    //next we need to take into account line elements

    //firstly we add contribiution from each line element as n_i *Sy_i
    //where n_i = E_i/Eref ; E_i is module of elasticity of ith line element
    //Sy_i is first moment of area in relation to y axis of ith line element
    for(int i=0; i < lines.length(); i++)
    {
        //module of elasticity of ith line element
        double E_i = lines[i]->getMaterial()->getModuleOfElasticity_E();
        //weight of Sy_i of ith line element
        double n_i = E_i/E; //where E i reference module of elasticity
        //first moment of area Sy_i of ith line element
        double Sy_i = lines[i]->firstMomentOfArea_Sy(p);

        //simple adding all line elements first moment of area Sx_i contribiutions multiplied by its weight
        weightedSy += n_i*Sy_i;

        //next we need to subtruct contribiution to weighted Sy from surface fragment overlapped by line element
        //in order to do that we pass through each surface from top to bottom and check whether it contains origin of line elemnt
        //if we find such surface we subtract its contribiution and break out from loop
        //IMPORTANT! WE ASSUME THAT LINE ELEMENT CAN BE EMBEDDED ONLU IN ONE SURFACE Si

        for(int j= surfaces.length()-1;  j>=0; --j)
        {
            Point *lineOrigin = lines[i]->getOrigin();
            if(surfaces[j]->pointIsInsieThePolygon(lineOrigin))
            {
                //module of elasticity of surface overlapped by line element
                double E_j = surfaces[j]->getMaterial()->getModuleOfElasticity_E();
                //weight of Sy_j of this surface
                double n_j = E_j/E;
                //first moment of area of this overlapped surface is equal to Sx_i (of ith line element)
                double Sy_j = Sy_i;

                //simple subtructing surface fragment contribiution to Sx which is overlapped by ith line element
                weightedSy -= n_j*Sy_j;
                break;
            }
        }

    }

    //next we need to take into account fiber groups
    //firstly we add contribiution from each fiber group as n_i *Sy_i
    //where n_i = E_i/Eref ; E_i is module of elasticity of ith fiber group
    //Sy_i is first moment of area in relation to x axis of ith fiber group
    for(int i=0; i < fiberGroups.length(); i++)
    {
        //module of elasticity of ith fiber group
        double E_i = fiberGroups[i]->getMaterial()->getModuleOfElasticity_E();
        //weight of Sy_i of ith fiber group
        double n_i = E_i/E; //where E i reference module of elasticity
        //first moment of area Sy_i of ith fiber group
        double Sy_i = fiberGroups[i]->firstMomentOfArea_Sy(p);

        //simple adding all fiber groups first moment of area Sy_i contribiutions multiplied by its weight
        weightedSy += n_i*Sy_i;

        //next we need to subtruct contribiution to weighted Sy from surface fragment overlapped by fiber group
        //in order to do that we pass through each surface from top to bottom and check whether it contains origin of fiber group
        //if we find such surface we subtract its contribiution and break out from loop
        //IMPORTANT! WE ASSUME THAT FIBER GROUP CAN BE EMBEDDED ONLY IN ONE SURFACE Si
        for(int j= surfaces.length()-1;  j>=0; --j)
        {
            Point *fiberGroupOrigin = fiberGroups[i]->getOrigin();
            if(surfaces[j]->pointIsInsieThePolygon(fiberGroupOrigin))
            {
                //module of elasticity of surface overlapped by ith fiber group
                double E_j = surfaces[j]->getMaterial()->getModuleOfElasticity_E();
                //weight of Sy_j of this surface
                double n_j = E_j/E;
                //first moment of area of this overlapped surface is equal to Sy_i (of ith fiber group)
                double Sy_j = Sy_i;

                //simple subtructing surface fragment contribiution to Sx which is overlapped by ith line element
                weightedSy -= n_j*Sy_j;

                break;
            }
        }
    }

    fprintf(stderr, "Calculated weight Sy (total): %g\n", weightedSy);

    return weightedSy;

}


/**
 * @brief Section::calculateIntersectionSy
 * @param surface - surface which will be intersected with remaining
 *                  surfaces indicated by indices in vector<int> list
 *                  starting from beginItr to endItr
 * @param beginItr - indicates index of first surface to intersect with surface (first argument)
 * @param endItr - last surface to intersect with surface (first argument)
 * @param referencePoint = point in relation to we calculate Sy
 * @param logString - only testing intersection combinations list
 * @return - cumulative Sy of surface intersection
 * This method is called recursively in order to calculate Sy (first moment of area)
 * of intersection surface (passed as first arg) with all other surfaces indicated
 * by vector<int>
 *1. intersect surface with first surface on the list (vector<int>)
 *2. foreach intersection (potentially we can get >=0 surfaces)  do:
 *    3. recursively intersect this new surface with remaining surfaces on the list (vector<int>)
 *       calculateIntersectionSy(intersection[i], ++beginItr, endItr);
 *    4. cumulate Sy (double - first moment of area) returned from recursively called function into double intersectionSy
 *    5. return this cumulative Sy
 * 6. Ending condition is when beginItr == endItr then this function return just Sy of surface
 */
double Section::calculateIntersectionSy(Surface *surface,
                               std::vector<int>::const_iterator beginItr,
                               std::vector<int>::const_iterator endItr,
                               Point *referencePoint,
                               std::string& logString)
{

    //if there isn't any remaining surface to intersect with surface passed as first argument
    if(beginItr == endItr)
        //then just return Sx of surface
        return surface->firstMomentOfArea_Sy(referencePoint);

    //define variable cumulating intersections Sy
    double intersectionSy = 0.0;

    //do intersection of surfaces (first parameter) with the first surface on vector<int> list
    int idx = *beginItr -1;//adjustment from {1-n} to {0-(n-1)} range

    QVarLengthArray<Surface *> *intersections = surface->intersection(surfaces[idx]);
    //logging intersection
    char tmpString[10];
    sprintf(tmpString, "xSx%d", idx);
    logString += tmpString;

    //foreach intersection surfaces (we can get >=0 surfaces) recursively call
    //this function to intersect it with remaining surfaces in the list vector<int>
    //we need to increment beginItr
    //return Sy values are cumulated in intersectionSy variable
    ++beginItr;

    for(int i=0; i < intersections->length(); i++) {
        Surface *intersection = intersections->at(i);
        intersectionSy += calculateIntersectionSy(intersection, beginItr, endItr, referencePoint, logString);

    }

    return intersectionSy;
}

/**
 * @brief Section::calculateGeometricalCenterX
 * @return
 * Function calculates centroid x coordinate using formula:
 * Cx = Sy* /A* where
 * Sy* - is weighted first moment of area of section
 * A* - is weighted area of section
 */
double Section::calculateGeometricalCenterX(void)
{
        double Cx = 0.0;

        Cx = weightedFirstMomentOfArea_Sy()/weightedArea();

        return Cx;
}

/**
 * @brief Section::calculateGeometricalCenterY
 * @return
 * Function calculates centroid y coordinate using formula:
 * Cy = Sx* /A* where
 * Sy* - is weighted first moment of area of section
 * A* - is weighted area of section
 */
double Section::calculateGeometricalCenterY(void)
{
    double Cy = 0.0;

    Cy = weightedFirstMomentOfArea_Sx()/weightedArea();

    return Cy;
}

/**
 * @brief Section::calculateGeometricalCenter
 * @return
 * Function gather Centroid Cx, Cy into Point object and returns it
 */
Point *Section::calculateGeometricalCenter()
{
    double Cx = calculateGeometricalCenterX();
    double Cy = calculateGeometricalCenterY();

    return new Point(Cx, Cy);
}

/**
 * @brief Section::getGeometricalCenterX
 * @return
 * Lazy calculation of centroid x coordinate
 */
double Section::getGeometricalCenterX(void)
{
     return getGeometricalCenter()->getX();

}

/**
 * @brief Section::getGeometricalCenterY
 * @return
 *Lazy calculation of centroid y coordinate
 */
double Section::getGeometricalCenterY(void)
{
    return getGeometricalCenter()->getY();

}

/**
 * @brief Section::getGeometricalCenter
 * @return
 * Lazy calculaton of centroid point
 */
Point *Section::getGeometricalCenter(void)
{
    if(centroid == NULL)
        centroid = calculateGeometricalCenter();

     return centroid;
}


QGraphicsEllipseItem *Section::drawCentoidOnScene(QGraphicsScene *scene)
{
    double rad = 1;

    Point *point = getGeometricalCenter();

    return scene->addEllipse(point->getX()-rad, -(point->getY()-rad), rad*5.0, rad*5.0,
                QPen(QColor(Qt::red)), QBrush(Qt::SolidPattern));
}

/**
 * @brief Section::yMinInCoordinateSystem
 * @param origin
 * @param rotationAngle
 * @return
 * Method is searching minimal y_MIN for the section in local coordinate
 * system moved to origin point and rotated by angle rotationAngle.
 * this value will be used to define STRAIN PROFILE in this rotated coordinate
 * system. For y_MIN there will be assumed for example limit strain corresponding to
 * ultimate tensile strength eps_tu.
 */
double Section::yMinInCoordinateSystem(Point *origin, double rotationAngle)
{
    double y_min;

    if(surfaces.length() < 1)
        throw NoSurfacesDefinedException("Nie zdefiniowanow przekoju żadnej powierzchni.");

    y_min = surfaces[0]->yMinInCoordinateSystem(origin, rotationAngle);

    for(int i=0; i<surfaces.length(); i++)
    {
        double y_min_i = surfaces[i]->yMinInCoordinateSystem(origin, rotationAngle);
        if( y_min_i < y_min)
                y_min = y_min_i;
    }

    return y_min;
}

/**
 * @brief Section::yMaxInCoordinateSystem
 * @param origin
 * @param rotationAngle
 * @return
 * Method is searching maximal y_MAX for the section in local coordinate system
 * moved to origin point and rotated by angle rotationAngle.
 * This value will be used to define STRAIN PROFILE in this rotated coordinate
 * system. For y_MAX there will be assumed for example limit strain corresponding to
 * ultimate compressive strength eps_cu
 */
double Section::yMaxInCoordinateSystem(Point *origin, double rotationAngle)
{
    double y_max;

    if(surfaces.length() < 1)
        throw NoSurfacesDefinedException("Nie zdefiniowanow przekoju żadnej powierzchni.");

    y_max = surfaces[0]->yMaxInCoordinateSystem(origin, rotationAngle);

    for(int i=0; i<surfaces.length(); i++)
    {
        double y_max_i = surfaces[i]->yMaxInCoordinateSystem(origin, rotationAngle);
        if( y_max_i > y_max)
                y_max = y_max_i;
    }

    return y_max;
}

/*************************************************************************************/
/* SUMMATION INTERNAL ACTIONS
 * Section is composed of a specific number of components i.e.
 * nS - number of surfaces
 * nL - number of lines
 * nFG - number of fibergroups
 * Resultant actins R (N, Mx, My) are the summation of individual component contribiution
 * R = SUM(i=1, nS, sign(Si)*Rs,i) + SUM(i=1, nL, sign(Li)*R_L,i) + SUM(i=1, nFG, sign(FGi)*R_fg,i)
 *ex. N
 * N = SUM(i=1, nS, sign(Si)*Ns,i) + SUM(i=1, nL, sign(Li)*N_l,i) + SUM(i=1, nFG, sign(FGi)*N_fg,i)
 * etc. for Mx, My
 **/

/**
 * @brief Section::summateInternalActions_N
 * @param profile - StrainProfile for which function summate calculated
 * internal actions N from each surface, line and fiber group componnent
 * reducing contributions from intersecting (overlapping) surfaces, lines, fibergroups
 * @return - double internal action N calculated as contribiution from each component
 * Method has as its first parameter StrianProfile which contains information about
 * lower and upper limit strains and about strain profile configuration in relation to section
 * i.e. rotation angle, y_min, y_max, etc. Based on this information method can summate
 * internal actions N from each component (surfaces, line elements, fiber groups)
 * calling its function internalAction_N() - which executes stress integration on given subcomponent
 * While summating internal actions we must take into account contributions from overlapping surfaces
 * so we use techinque such as in case of weighted area, weighted first moment of area
 * i.e. we must use combinatorics to generate all intersections and take its internal action N
 * as positively or negatively contributing to summateInternalActions_N
 */
double Section::summateInternalActions_N(StrainProfile *profile)
{
    double N = 0.0; //cumulating N resulting internal action

    if(surfaces.length() < 1) return N;

    //pass through each surface from background to foreground material
    for(int i=0; i< surfaces.length(); i++)
    {
        //adding up internal action N contribiution from each surface component
        double surfaceN = surfaces[i]->internalAction_N(profile);
        N += surfaceN;

        if(DEBUG)
            fprintf(stderr, "Adding S%d surface influance - N = %g\n", i, surfaceN);
    }

    //taking into account positive or negative contribiution to internal action N
    //from surfaces intersections
    //N +/- = ...

    //creating 2dimensional matrix containing all possible surfaces intersection combinations
    CombinationGenerator generator(surfaces.length());
    std::vector< std::vector<int> > combinations = generator.getAllCombinations();

    //now we need to make more sophisticated calculation
    //in order to subtract/add contribiution form overlapping surfaces
    char tmpString[1000]; //storing logging info about generated intersections
    std::string logString(""); //also for logging
    double intersectionN = 0.0; //variable storing contribiution form each intersection on internal action N

    for(auto itr = combinations.begin();
        itr < combinations.end(); itr++)
    {
        std::vector<int> row = *itr;
        //we determine the sign of calculations from current intersection combination
        //it depends on number of intersected surfaces
        double sign = std::pow(-1.0, (int) (row.size()-1));
        //in the case of internal actions we don't need weight coefficient
        //but we need to have material definition which will be assigned to
        //resultant intersection surface:
        int materialIdx = row[1]-1;
        Material *material = surfaces[materialIdx]->getMaterial();
        sprintf(tmpString, "%dx%s", (int) sign, material->getMaterialName().c_str());
        logString += tmpString;

        //we are retrieving foreground surface for which we will be
        //checking intersections with surfaces laying underneath it
        Surface *firstSurface = surfaces[row[0]-1];
        sprintf(tmpString, "xS%d", row[0]-1);
        logString +=tmpString;

        intersectionN = calculateIntersectionN(firstSurface, //foreground Surface
                                               ++row.begin(), //index of start surface to intersect with
                                               row.end(), //index of end surface to intesect with
                                               profile, //strain profile information
                                               material, //material of resulting surfaces intersection
                                               logString);

        sprintf(tmpString, " (intersection sign= %g, N = %g)\n", sign, intersectionN);
        logString += tmpString;

        N += sign*intersectionN;
    }
    if(DEBUG) {
        fprintf(stderr, "%s", logString.c_str());
        fprintf(stderr, "Calculated section total N (only surfaces): %g\n", N);
    }
    //contribiutions from line element

    //firstly we add contribution form each line element Li
    //and simultanously we are searching for overlapped materail
    //and subtruct contribution of this material by constructing
    //fictious line element of this material and subtructing
    //its internal action N contribution
    for(int i=0; i< lines.length(); i++)
    {
        N += lines[i]->internalAction_N(profile);

        //next we need to subtract contribiution from the surface fragment
        //overlapped by this line element. To do this we pass through
        //each surface element from top to bottom and checking whether
        //line origin belongs to jth surface. If so  we construct
        //fictious line element based on current ith line element
        //and jth surface material and substract from total N
        //contribution from this fictious line element
        //IMPORTANT! WE ASSUME THAT LINE ELEMENT CAN BE EMBEDDED ONLY IN ONE SURFACE ELEMENT
        for(int j= surfaces.length()-1; j>=0; --j)
        {
            Point *lineOrigin = lines[i]->getOrigin();
            if(surfaces[j]->pointIsInsieThePolygon(lineOrigin))
            {
                //we construct new line element based on Sj->material()
                //and ith line element
                Line  lineIntersecting("Intersection", surfaces[j]->getMaterial());
                //copying points and areas
                Point **points = lines[i]->getPointsArray();
                double *areas = lines[i]->getAreas();
                for(int k=0; k<lines[i]->numberOfPoints(); ++k) {
                    lineIntersecting.addPointAndArea(points[k]->getX(), points[k]->getY(), areas[k]);
                }
                //subtracting contribution from surface fragment overlapped by line element
                //in this case we assume this fragment as equal to ovelapping line element
                //but with material of overlapped surface
                N -= lineIntersecting.internalAction_N(profile);

                break;
            }
        }
    }

    //contributions from fiber groups

    //firstly we add contribiution form each fiber group FGi
    //and simultanously we are searching for overlapped materail
    //and substruct contribiution of this material by constructing
    //fictious fiber group elemet of this material and subtructing
    //its internal action N contribiution
    for(int i=0; i< fiberGroups.length(); i++)
    {
        N += fiberGroups[i]->internalAction_N(profile);

        //next we need to subtract contribution from the surface fragment
        //overlapped by this fiber group. To do this we pass through
        //each surface element from top to bottom and checking whether
        //fiber group origin belongs to jth surface. If so  we construct
        //fictious fiber group element based on current ith fiber group
        //and jth surface material and substract from total N
        //contribiution from this fictious fiber group element
        //IMPORTANT! WE ASSUME THAT FIBER GROUP CAN BE EMBEDDED ONY IN ONE SURFACE ELEMENT
        for(int j= surfaces.length()-1; j>=0; --j)
        {
            Point *fiberGroupOrigin = fiberGroups[i]->getOrigin();
            if(surfaces[j]->pointIsInsieThePolygon(fiberGroupOrigin))
            {
                //we construct new fiber group element based on Sj->material()
                //and ith fiber group
                FiberGroup fiberGroupIntersecting("Intersection", surfaces[j]->getMaterial());
                //copying points and areas
                Point **points = fiberGroups[i]->getPointsArray();
                double *areas = fiberGroups[i]->getAreas();
                for(int k=0; k<fiberGroups[i]->numberOfPoints(); ++k) {
                    fiberGroupIntersecting.addPointAndArea(points[k]->getX(), points[k]->getY(), areas[k]);
                }
                //subtracting contribiution from surface fragment overlapped by fiber group element
                //in this case we assume this fragment as equal to ovelapping fiber group
                //but with material of overlapped surface
                N -= fiberGroupIntersecting.internalAction_N(profile);

                break;
            }
        }
    }
    if(DEBUG) {
        fprintf(stderr, "Calculated section total N: %g\n", N);
        fprintf(stderr, "--------------------------------------------------------\n");
    }
    return N; //[N]

}

/**
 * @brief Section::calculateIntersectionN
 * @param surface - surface which will be intersected with remaining
 *                  surfaces indcated by indices in vector<int> list
 *                  starting from beginItr to endItr
 * @param beginItr - indicates index of first surface to intersect with first surface (first argument)
 * @param endItr - last surface to intersect with first surface (first argument)
 * @param profile - strain profile
 * @param material - material definition which should be assumed for intersection surface
 * @param logString - only testing intersection combinations list
 * @return - intersection contribiution to internal action N
 */
double Section::calculateIntersectionN(Surface *surface,
                                       std::vector<int>::const_iterator beginItr,
                                       std::vector<int>::const_iterator endItr,
                                       StrainProfile *profile,
                                       Material *material,
                                       std::string &logString)
{
    //if there isn't any remaining surface to intersect with  first surface (passed as first argument)
    if(beginItr == endItr) {
            //then just return surface contribiution to internal action N
            //we must create intersection surface containing copy of
            //surface passed as first argument with modified material
            //set to material passed as parameter to function
            //we do this because of we do not want to mess materaial definition
            //in surface hierarchy
            Surface intersection("Intersection", material);
            Point **points = surface->getPointsArray();
            double *angles = surface->getAngels();
            for(int i=0; i < surface->numberOfPoints(); i++)
            {
                 intersection.addPointAndAngel(points[i]->getX(), points[i]->getY(), angles[i]);
            }
            //intersection of two polygons can only be polygon so we set that this intersection
            //surfaces is adaptive linearized in order to reduce overhead from performing unneccessery adaptive linearization process
            intersection.setIsAdaptiveLinearized(true);

           return intersection.internalAction_N(profile);

    }

    //define variable cumulating intersections contribiution to internal action N
    double intersectionN = 0.0;

    //do intersection of surfaces (first parameter) with the first surface
    //on vector<int> list
    int idx = *beginItr -1; //adjustment from {1-n} to {0-(n-1)} range

    QVarLengthArray<Surface *> *intersections = surface->intersection(surfaces[idx]);
    //logging intersection
    char tmpString[10];
    sprintf(tmpString, "xS%d", idx);
    logString += tmpString;


    //foreach intersection surfaces (we can get >=0 surfaces) recursively call
    //this function to intersect it with remaining surfaces in the list vector<int>
    //we need to increment beginItr
    //returned contribiution to internal action N are cumulated in intersectionN variable

    ++beginItr;
    for(int i=0; i < intersections->length(); i++) {
        Surface *intersection = intersections->at(i);
        intersectionN += calculateIntersectionN(intersection,
                                                beginItr, endItr,
                                                profile,
                                                material,
                                                logString);
    }

    return intersectionN;
}

/**
 * @brief Section::summateInternalActions_Mx
 * @param profile - StrainProfile for which function summate calculated
 * internal actions Mx from each surface, line and fiber group componnent
 * reducing contributions from intersecting (overlapping) surfaces, lines, fibergroups
 * @return - double internal action Mx calculated as contribiution from each component
 * Method has as its first parameter StrianProfile which contains information about
 * lower and upper limit strains and about strain profile configuration in relation to section
 * i.e. rotation angle, y_min, y_max, etc. Based on this information method can summate
 * internal actions Mx from each component (surfaces, line elements, fiber groups)
 * calling its function internalAction_Mx() - which executes stress integration on given subcomponent
 * While summating internal actions we must take into account contributions from overlapping surfaces
 * so we use techinque such as in case of weighted area, weighted first moment of area
 * i.e. we must use combinatorics to generate all intersections and take its internal action Mx
 * as positively or negatively contributing to summateInternalActions_Mx
 */
double Section::summateInternalActions_Mx(StrainProfile *profile)
{
    double Mx = 0.0; //cumulating Mx resulting internal action

    if(surfaces.length() < 1) return Mx;

    //pass through each surface from background to foreground material
    for(int i=0; i< surfaces.length(); i++)
    {
        //adding up internal action Mx contribiution from each surface component
        double surfaceMx = surfaces[i]->internalAction_Mx(profile);
        Mx += surfaceMx;

        if(DEBUG)
            fprintf(stderr, "Adding S%d surface influance - Mx = %g\n", i, surfaceMx);
    }

    //taking into account positive or negative contribiution to internal action Mx
    //from surfaces intersections
    //Mx +/- = ...

    //creating 2dimensional matrix containing all possible surfaces intersection combinations
    CombinationGenerator generator(surfaces.length());
    std::vector< std::vector<int> > combinations = generator.getAllCombinations();

    //now we need to make more sophisticated calculation
    //in order to subtruct/add contribution from overlapping surfaces
    char tmpString[1000]; //storing logging info about generated intersection combination
    std::string logString(""); //also for logging
    double intersectionMx = 0.0; //variable storing contribution from each intersection on internal action Mx

    for(auto itr = combinations.begin();
        itr < combinations.end(); itr++)
    {
        std::vector<int> row = *itr;
        //we determine the sign of calculations from current intersection combination
        //it depends on number of intersected surfaces
        double sign = std::pow(-1.0, (int) (row.size()-1));
        //in the case of internal actions we don't need weight coefficient
        //but we need to have material definition which will be assigned to
        //resultant intersection surface:
        int materialIdx = row[1]-1;
        Material *material = surfaces[materialIdx]->getMaterial();
        sprintf(tmpString, "%dx%s", (int) sign, material->getMaterialName().c_str());
        logString +=tmpString;

        //we are retrieving foreground surface for which we will be
        //checking intersections with surfaces laying underneath it
        Surface *firstSurface = surfaces[row[0]-1];
        sprintf(tmpString, "xS%d", row[0]-1);
        logString += tmpString;

        intersectionMx = calculateIntersectionMx(firstSurface, //foreground Surface
                                                 ++row.begin(), //index of start surface to intersect with
                                                 row.end(), //index of end surface to intersect with
                                                 profile, //strain profile
                                                 material, //material of resulting surfaces intersection
                                                 logString);

        sprintf(tmpString, " (intersection sign= %g, Mx = %g)\n", sign, intersectionMx);
        logString += tmpString;

        Mx += sign*intersectionMx;

    }

    if(DEBUG) {
        fprintf(stderr, "%s", logString.c_str());
        fprintf(stderr, "Calculated section total Mx (only surfaces): %g\n", Mx);
    }
    //contribiutions from line element

    //firstyl we add contribution from each line element Li
    //and simultanously we are searching for ovelapped material
    //and subtruct contribution of this material by constructing
    //fictious line element of this material and subtructing
    //its internal action Mx contribution
    for(int i=0; i < lines.length(); i++)
    {
        Mx += lines[i]->internalAction_Mx(profile);

        //next we need to subtract contribution form the surface fragment
        //overlapped by this line element. To do this we pass through
        //each surface element from top to bottom and checking whether
        //line origin belongs to jth surface. If so we construct
        //fictious line element based on current ith line element
        //and jth surface material and subtruct from total Mx
        //contribution from this fictious line element
        //IMPORTANT! WE ASSUME THAT LINE ELEMENT CAN BE EMBEDDED ONLY IN ONE SURFACE ELEMENT
        for(int j=surfaces.length()-1; j>=0; --j)
        {
            Point *lineOrigin = lines[i]->getOrigin();
            if(surfaces[j]->pointIsInsieThePolygon(lineOrigin))
            {
                //we construct new line element based on Sj->material()
                //and ith line element
                Line lineIntersecting("Intersection", surfaces[j]->getMaterial());
                //copying points and areas
                Point **points = lines[i]->getPointsArray();
                double *areas = lines[i]->getAreas();
                for(int k=0; k<lines[i]->numberOfPoints(); ++k)
                {
                    lineIntersecting.addPointAndArea(points[k]->getX(), points[k]->getY(), areas[k]);
                }
                //subtructing contribution from surface fragment overlapped by line element
                //in this case we assume this fragment as equal to overlapping line element
                //but with material of overllaped surface
                Mx -= lineIntersecting.internalAction_Mx(profile);

                break;
            }
        }
    }

    //contribution from fiber groups

    //firstly we add contribiution form each fiber group FGi
    //and simultanously we are searching for overlapped material
    //and substruct contribiution of this material by constructing
    //fictious fiber group elemet of this material and subtructing
    //its internal action Mx contribiution
    for(int i=0; i < fiberGroups.length(); i++)
    {
        Mx += fiberGroups[i]->internalAction_Mx(profile);

        //next we need to subtract contribution from the surface fragment
        //overlapped by this fiber group. To do this we pass through
        //each surface element from top to bottom and checking whether
        //fiber group origin belongs to jth surface. If so  we construct
        //fictious fiber group element based on current ith fiber group
        //and jth surface material and substract from total Mx
        //contribution from this fictious fiber group element
        //IMPORTANT! WE ASSUME THAT FIBER GROUP CAN BE EMBEDDED ONY IN ONE SURFACE ELEMENT
        for(int j= surfaces.length()-1; j>=0; --j)
        {
            Point *fiberGroupOrigin = fiberGroups[i]->getOrigin();
            if(surfaces[j]->pointIsInsieThePolygon(fiberGroupOrigin))
            {
                //we construct fiber group element based on Sj->material()
                //and ith fiber group
                FiberGroup fiberGroupIntersecting("Intersection", surfaces[j]->getMaterial());
                //copying points and areas
                Point **points = fiberGroups[i]->getPointsArray();
                double *areas = fiberGroups[i]->getAreas();
                for(int k=0; k<fiberGroups[i]->numberOfPoints(); ++k)
                {
                    fiberGroupIntersecting.addPointAndArea(points[k]->getX(), points[k]->getY(), areas[k]);
                }
                //subtracting contribution from surface fragment overlapped by fiber group element
                //in this case we assume this fragment as equal to overlapping fiber group
                //but with material of overlapped surface
                Mx -= fiberGroupIntersecting.internalAction_Mx(profile);

                break;
            }
        }
    }
    if(DEBUG) {
        fprintf(stderr, "Calculated section total Mx: %g\n", Mx);
        fprintf(stderr, "--------------------------------------------------------\n");
    }

    return Mx; //[Nmm]
}

double Section::calculateIntersectionMx(Surface *surface,
                                        std::vector<int>::const_iterator beginItr,
                                        std::vector<int>::const_iterator endItr,
                                        StrainProfile *profile,
                                        Material *material,
                                        std::string &logString)
{
    //if there isn't any remaining surface to intersect with first surfce (passed as first argument)
    if(beginItr == endItr) {
        //then just return surface contribiution to internal action Mx
        //we must create intersection surface containing copy of
        //surface passed as first argument with modified material
        //set to material passed as parameter to function
        //we do this because of we do not want to mess materaial definition
        //in surface hierarchy
        Surface intersection("Intersection", material);
        Point **points = surface->getPointsArray();
        double *angles = surface->getAngels();
        for(int i=0; i < surface->numberOfPoints(); i++)
        {
             intersection.addPointAndAngel(points[i]->getX(), points[i]->getY(), angles[i]);
        }
        //intersection of two polygons can only be polygon so we set that this intersection
        //surfaces is adaptive linearized in order to reduce overhead from performing unneccessery adaptive linearization process
        intersection.setIsAdaptiveLinearized(true);

       return  intersection.internalAction_Mx(profile);

    }

    //define variable cumulating intersections contribiution to internal action Mx
    double intersectionMx = 0.0;

    //do intersection of surfaces (first parameter) with the first surface
    //on vector<int> list
    int idx = *beginItr -1; //adjustment from {1-n} to {0-(n-1)} range

    QVarLengthArray<Surface *> *intersections = surface->intersection(surfaces[idx]);
    //logging intersection
    char tmpString[10];
    sprintf(tmpString, "xS%d", idx);
    logString += tmpString;

    //foreach intersection surfaces (we can get >=0 surfaces) recursively call
    //this function to intersect it with remaining surfaces in the list vector<int>
    //we need to increment beginItr
    //returned contribiution to internal action Mx are cumulated in intersectionMx variable
    ++beginItr;
    for(int i=0; i < intersections->length(); i++) {
        Surface *intersection = intersections->at(i);
        intersectionMx += calculateIntersectionMx(intersection,
                                                beginItr, endItr,
                                                profile,
                                                material,
                                                logString);
    }

    return intersectionMx;
}

/**
 * @brief Section::summateInternalActions_My
 * @param profile - StrainProfile
 * @return summated internal action My
 * Method works in analogical way as summateInternalAxtions_Mx()
 */
double Section::summateInternalActions_My(StrainProfile *profile)
{
    double My = 0.0; //cumulating Mx resulting internal action

     //pass through each surface from background to foreground material
    for(int i=0; i< surfaces.length(); i++)
    {
        //adding up internal action My contribiution from each surface component
        double surfaceMy = surfaces[i]->internalAction_My(profile);
        My += surfaceMy;

        if(DEBUG)
            fprintf(stderr, "Adding S%d surface influance - My = %g\n", i, surfaceMy);
    }

    //taking into account positive or negative contribiution to internal action My
    //from surfaces intersections
    //My +/- = ...

    //creating 2dimensional matrix containing all possible surfaces intersection combinations
    CombinationGenerator generator(surfaces.length());
    std::vector< std::vector<int> > combinations = generator.getAllCombinations();

    //now we need to make more sophisticated calculation
    //in order to subtruct/add contribution from overlapping surfaces
    char tmpString[1000];  //storing logging info about generated intersection combination
    std::string logString(""); //also for logging
    double intersectionMy = 0.0; //variable storing contribution from each intersection on internal action My

    for(auto itr = combinations.begin();
        itr < combinations.end(); itr++)
    {
        std::vector<int> row = *itr;
        //we determine the sign of calculations from current intersection combination
        //it depends on number of intersected surfaces
        double sign = std::pow(-1.0, (int) (row.size()-1));
        //in the case of internal actions we don't need weight coefficient
        //but we need to have material definition which will be assigned to
        //resultant intersection surface:
        int materialIdx = row[1]-1;
        Material *material = surfaces[materialIdx]->getMaterial();
        sprintf(tmpString, "%dx%s", (int) sign, material->getMaterialName().c_str());
        logString +=tmpString;

        //we are retrieving foreground surface for which we will be
        //checking intersections with surfaces laying underneath it
        Surface *firstSurface = surfaces[row[0]-1];
        sprintf(tmpString, "xS%d", row[0]-1);
        logString += tmpString;

        intersectionMy = calculateIntersectionMy(firstSurface, //foreground Surface
                                                 ++row.begin(), //index of start surface to intersect with
                                                 row.end(), //index of end surface to intersect with
                                                 profile, //strain profile
                                                 material, //material of resulting surfaces intersection
                                                 logString);

        sprintf(tmpString, " (intersection sign= %g, My = %g)\n", sign, intersectionMy);
        logString += tmpString;

        My += sign*intersectionMy;
    }
    if(DEBUG) {
        fprintf(stderr, "%s", logString.c_str());
        fprintf(stderr, "Calculated section total My (only surfaces): %g\n", My);
    }
    //contribiutions from line element

    //firstly we add contribution from each line element Li
    //and simultanously we are searching for ovelapped material
    //and subtruct contribution of this material by constructing
    //fictious line element of this material and subtructing
    //its internal action My contribution
    for(int i=0; i < lines.length(); i++)
    {
        My += lines[i]->internalAction_My(profile);

        //next we need to subtract contribution form the surface fragment
        //overlapped by this line element. To do this we pass through
        //each surface element from top to bottom and checking whether
        //line origin belongs to jth surface. If so we construct
        //fictious line element based on current ith line element
        //and jth surface material and subtruct from total My
        //contribution from this fictious line element
        //IMPORTANT! WE ASSUME THAT LINE ELEMENT CAN BE EMBEDDED ONLY IN ONE SURFACE ELEMENT
        for(int j=surfaces.length()-1; j>=0; --j)
        {
            Point *lineOrigin = lines[i]->getOrigin();
            if(surfaces[j]->pointIsInsieThePolygon(lineOrigin))
            {
                //we construct new line element based on Sj->material()
                //and ith line element
                Line lineIntersecting("Intersection", surfaces[j]->getMaterial());
                //copying points and areas
                Point **points = lines[i]->getPointsArray();
                double *areas = lines[i]->getAreas();
                for(int k=0; k<lines[i]->numberOfPoints(); ++k)
                {
                    lineIntersecting.addPointAndArea(points[k]->getX(), points[k]->getY(), areas[k]);
                }
                //subtructing contribution from surface fragment overlapped by line element
                //in this case we assume this fragment as equal to overlapping line element
                //but with material of overllaped surface
                My -= lineIntersecting.internalAction_My(profile);

                break;
            }
        }
    }

    //contribution from fiber groups

    //firstly we add contribiution form each fiber group FGi
    //and simultanously we are searching for overlapped material
    //and substruct contribiution of this material by constructing
    //fictious fiber group elemet of this material and subtructing
    //its internal action My contribiution
    for(int i=0; i < fiberGroups.length(); i++)
    {
        My += fiberGroups[i]->internalAction_My(profile);

        //next we need to subtract contribution from the surface fragment
        //overlapped by this fiber group. To do this we pass through
        //each surface element from top to bottom and checking whether
        //fiber group origin belongs to jth surface. If so  we construct
        //fictious fiber group element based on current ith fiber group
        //and jth surface material and substract from total My
        //contribution from this fictious fiber group element
        //IMPORTANT! WE ASSUME THAT FIBER GROUP CAN BE EMBEDDED ONY IN ONE SURFACE ELEMENT
        for(int j= surfaces.length()-1; j>=0; --j)
        {
            Point *fiberGroupOrigin = fiberGroups[i]->getOrigin();
            if(surfaces[j]->pointIsInsieThePolygon(fiberGroupOrigin))
            {
                //we construct fiber group element based on Sj->material()
                //and ith fiber group
                FiberGroup fiberGroupIntersecting("Intersection", surfaces[j]->getMaterial());
                //copying points and areas
                Point **points = fiberGroups[i]->getPointsArray();
                double *areas = fiberGroups[i]->getAreas();
                for(int k=0; k<fiberGroups[i]->numberOfPoints(); ++k)
                {
                    fiberGroupIntersecting.addPointAndArea(points[k]->getX(), points[k]->getY(), areas[k]);
                }
                //subtracting contribution from surface fragment overlapped by fiber group element
                //in this case we assume this fragment as equal to overlapping fiber group
                //but with material of overlapped surface
                My -= fiberGroupIntersecting.internalAction_My(profile);

                break;
            }
        }
    }
    if(DEBUG) {
        fprintf(stderr, "Calculated section total My: %g\n", My);
        fprintf(stderr, "--------------------------------------------------------\n");
    }

    return My; //[Nmm]
}

double Section::calculateIntersectionMy(Surface *surface,
                                        std::vector<int>::const_iterator beginItr,
                                        std::vector<int>::const_iterator endItr,
                                        StrainProfile *profile,
                                        Material *material,
                                        std::string &logString)
{
    //if there isn't any remaining surface to intersect with first surfce (passed as first argument)
    if(beginItr == endItr) {
        //then just return surface contribiution to internal action My
        //we must create intersection surface containing copy of
        //surface passed as first argument with modified material
        //set to material passed as parameter to function
        //we do this because of we do not want to mess materaial definition
        //in surface hierarchy
        Surface intersection("Intersection", material);
        Point **points = surface->getPointsArray();
        double *angles = surface->getAngels();
        for(int i=0; i < surface->numberOfPoints(); i++)
        {
             intersection.addPointAndAngel(points[i]->getX(), points[i]->getY(), angles[i]);
        }
        //intersection of two polygons can only be polygon so we set that this intersection
        //surfaces is adaptive linearized in order to reduce overhead from performing unneccessery adaptive linearization process
        intersection.setIsAdaptiveLinearized(true);

       return intersection.internalAction_My(profile);
    }

    //define variable cumulating intersections contribiution to internal action My
    double intersectionMy = 0.0;

    //do intersection of surfaces (first parameter) with the first surface
    //on vector<int> list
    int idx = *beginItr -1; //adjustment from {1-n} to {0-(n-1)} range

    QVarLengthArray<Surface *> *intersections = surface->intersection(surfaces[idx]);
    //logging intersection
    char tmpString[10];
    sprintf(tmpString, "xS%d", idx);
    logString += tmpString;

    //foreach intersection surfaces (we can get >=0 surfaces) recursively call
    //this function to intersect it with remaining surfaces in the list vector<int>
    //we need to increment beginItr
    //returned contribiution to internal action My are cumulated in intersectionMy variable
    ++beginItr;
    for(int i=0; i < intersections->length(); i++) {
        Surface *intersection = intersections->at(i);
        intersectionMy += calculateIntersectionMy(intersection,
                                                beginItr, endItr,
                                                profile,
                                                material,
                                                logString);

    }

    return intersectionMy; //[Nmm]
}

/**
 * @brief Section::transferBackFromLocalMx
 * @param localMx - Mx resulting internal action calculated in integration process
 *                  in local (rotated by omega angle) coordinate system
 * @param localMy - My resulting internal action calculated in integration process
 *                  in local (rotated by omega angle) coordinate system
 * @param omega - angle of coordinate system rotation from XCY to local xCy
 *                (unit: radians)
 * @return Mx in XCY coordinate system it is global but central coordinate system
 * In this function we transform back Mx internal action to XCY so we use
 * -OMEGA angle (negative) because we rotate back in reverse direction!
 * FORMULA: M_X = Mx*cos(-0) + My*sin(-0), 0 - is omega angle
 */
double Section::transferBackFromLocalMx(double localMx, double localMy, double omega)
{
    double Mx = localMx*std::cos(-omega) + localMy*std::sin(-omega);
    return Mx; //[Nmm]
}

/**
 * @brief Section::transferBackFromLocalMy
 * @param localMx - Mx resulting internal action calculated in integration process
 *                 in local (rotated by omega angle) coordinate system
 * @param localMy - My resulting internal action calculated in integration process
 *                 in local (rotated by omega angle) coordinate system
 * @param omega - angle of coordinate system rotation from XCY to local xCy
 *               (unit: radians)
 * @return  My in XCY coordinate system it is global but central coordinate system
 * In this function we transform back My internal action to XCY so we use
 * -OMEGA angle (negative) because we rotate back in reverse direction!
 * FORMULA: M_Y = -Mx*sin(-0) + My*cos(-0), 0 - is omega angle
 */
double Section::transferBackFromLocalMy(double localMx, double localMy, double omega)
{
    double My = -localMx*std::sin(-omega) + localMy*std::cos(-omega);
    return My; //[Nmm]
}
