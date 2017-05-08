#include "fibergroup.h"
#include "strainprofile.h"
#include "gausslegendre.h"
#include "flags.h"

FiberGroup::FiberGroup(std::string fibergroupName)
{
    name = fibergroupName;
    //numberOfFibers = 0;
}

FiberGroup::FiberGroup(std::string fibergroupName, Material *material)
{
    name = fibergroupName;
    this->material = material;
    //numberOfFibers = 0;
}

void FiberGroup::addPointAndArea(double x, double y, double A) {
    addPoint(x,y);
    areas.append(A);
    //numberOfFibers++;
}
void FiberGroup::setMaterial(Material *m) {
    material = m;
}
Material* FiberGroup::getMaterial(void) {
    return material;
}
double *FiberGroup::getAreas(void) {
    return areas.data();
}
void FiberGroup::setName(std::string n)
{
   name = n;

}

std::string FiberGroup::getName(void)
{
    return name;
}

/*int FiberGroup::getNumberOfFibers(void)
{
    return numberOfFibers;
}*/

void FiberGroup::clearPointsAndAreas(void) {
    points.clear();
    areas.clear();
    //numberOfFibers = 0;
}

FiberGroup::~FiberGroup() {
   //delete material;
}

Point *FiberGroup::getOrigin(void)
{
    if(points.length() > 0)
        return this->points[0];
    //gdy w fiber group nie mamy zdefiniowanych żadnych punktów
    return NULL;
}

/*****************************************************
 * getTotalArea() calculates total area of fibergroup*
 * as sum of areas of each fiber in fibergroup.      *
 * A = A1 + A2 + .. + An                             *
 *****************************************************/
double FiberGroup::getTotalArea(void) {

    double totalArea = 0;
    double *areas = getAreas();

    for(int i=0; i<numberOfPoints(); i++) {
        totalArea += areas[i];
    }

    return totalArea;
}

/*******************************************************
 * getGeometricalCenterX() - this return x_c centroid  *
 * coordinate of fibergroup. It is calculated according*
 * to formula: x_c = Sy/Atotal                         *
 * Sy - is first moment of area in relation to y axie  *
 *      which is sum of Sy1 + Sy2 + ... + Syn of each  *
 *      fiber. (x_c, and Syi are in default coordinate *
 *      system => origin = (0,0)                       *
 *Atotal - is a total area of fibergroup: A1 + .. + An *
 *******************************************************/
double FiberGroup::getGeometricalCenterX(void)
{
    return firstMomentOfArea_Sy()/getTotalArea();
}

/*******************************************************
 * getGeometricalCenterY() - this return y_c centroid  *
 * coordinate of fibergroup. It is calculated according*
 * to formula: y_c = Sx/Atotal                         *
 * Sx - is first moment of area in relation to x axie  *
 *      which is sum of Sx1 + Sx2 + ... + Sxn of each  *
 *      fiber. (y_c, and Syi are in default coordinate *
 *      system => origin = (0,0)                       *
 * Atotal - is a total area of fibergroup: A1 + .. + An*
 *******************************************************/
double FiberGroup::getGeometricalCenterY(void)
{
    return firstMomentOfArea_Sx()/getTotalArea();
}

/*******************************************************
 * getGeometricalCenter() is function which uses two   *
 * above methods getGeometricalCenterX() and ....Y()   *
 * and wrap coordinates x_c and y_c into Point object  *
 *******************************************************/
Point *FiberGroup::getGeometricalCenter(void)
{
    double x_c = getGeometricalCenterX();
    double y_c = getGeometricalCenterY();

    return new Point(x_c, y_c);
}

/*******************************************************
 * firstMomentOfArea_Sx() calculates Sx in relation to *
 * x axie based on area Ai of each fiber and yc_i      *
 * centroid coordinate of ith fiber in group.          *
 * Sx = A1*(yc_1 - y)  + ... + An*(yc_n - y)           *
 * where:                                              *
 * y is coordinate of translated coordinate system in  *
 * relation to default coordinate system. If we calcul-*
 * ate this Sx in default coordinate system origin (0,0)
 * then y = 0.                                         *
 *******************************************************/
double FiberGroup::firstMomentOfArea_Sx() {
    return firstMomentOfArea_Sx(new Point(0.0, 0.0));
}

double FiberGroup::firstMomentOfArea_Sx(Point *p) {

    double Sx = 0;
    double *areas = getAreas();

    //Calculating Sx = A1*(yc_1 - y)  + ... + An*(yc_n - y)
    for(int i =0; i < numberOfPoints(); i++) {
        Sx += areas[i]*(getYCenterOfFiber(i) - p->getY()); // Ai*(yc_i - y);
    }

    return Sx;
}

/*******************************************************
 * firstMomentOfArea_Sy() calculate Sy in relation to  *
 * y axie based on area Ai of each fiber and xc_i      *
 * centroid coordinate of ith fiber in group.          *
 * Sy = A1*(xc_1 - x) + ... + An*(xc_n - x)            *
 * where:                                              *
 * x is coordinate of translated coordinate system in  *
 * relation to default coordinate system. If we calcul-*
 * ate this Sy in default coordinate system origin (0,0)
 * then x = 0.                                         *
 *******************************************************/
double FiberGroup::firstMomentOfArea_Sy()
{
    return firstMomentOfArea_Sy(new Point(0.0, 0.0));
}

double FiberGroup::firstMomentOfArea_Sy(Point *p) {
    double Sy = 0;
    double *areas = getAreas();

    //Calculating Sy = A1(xc_1 - x) + ... + An(xc_n - x)
    for(int i = 0; i < numberOfPoints(); i++) {
        Sy += areas[i]*(getXCenterOfFiber(i) - p->getX()); //Ai*(xc_i - x);
    }

    return Sy;
}

double FiberGroup::getXCenterOfFiber(int i)
{
    return points[i]->getX();
}

double FiberGroup::getYCenterOfFiber(int i)
{
    return points[i]->getY();
}

Point *FiberGroup::getCenterOfFiber(int i) {
    return new Point( getXCenterOfFiber(i), getYCenterOfFiber(i));
}

/***********************************************************
 * STRESS INTEGRATION METHODS                              *
 ***********************************************************/

double FiberGroup::internalAction_N(StrainProfile *profile)
{
    return performStressIntegration(profile, ACTION_N);
}

double FiberGroup::internalAction_Mx(StrainProfile *profile)
{
    return performStressIntegration(profile, ACTION_Mx);
}

double FiberGroup::internalAction_My(StrainProfile *profile)
{
    return -performStressIntegration(profile, ACTION_My);
}

/**
 * @brief FiberGroup::performStressIntegration
 * @param profile - StrainProfile for which to calculate given internal action (N, Mx, My)
 *                  for this surface
 * @param actionType - type of calculated internal action N, Mx, My
 * @return - internal action (N, Mx, My) value as double
 *  The stress integration procedure for each fiber group component FGi is
 * straightforward and performed analitically. Formula:
 * Rfg,i = SUM(j=1, nFGi, Aj*xj^r*yj^s*sigma(eps_o - fi*yj) )
 */
double FiberGroup::performStressIntegration(StrainProfile *profile, InternalAction actionType)
{

    //exponents of variable xj and yj
    //depends on calculated internal action N, Mx, My
    // (r,s) = (0,0) for N, (r,s) = (0,1) for Mx, (r,s) = (1,0) for -My
    double r = 0, s = 0;
    //calculated internal action
    double internalAction = 0.0;

    if(actionType == ACTION_N)
    {
         //stress integration scheme for internal action N
        r = 0; s = 0;
    } else if(actionType == ACTION_Mx)
    {
       //stress integration scheme for internal action Mx
        r = 0; s = 1;
    } else if(actionType == ACTION_My)
    {
        //stress integration scheme for internal action -My
        r = 1; s = 0;
    }

    //strain profile basic characteristics fi and epsilon_o
    double eps_o = profile->getStrainAtOrigin_eps_o();
    double fi = profile->getCurvature_fi();

    for(int j=0; j < numberOfPoints(); ++j)
    {
        //getting each contribiuting fiber area and center coordinates
        double Aj = areas[j];
        Point *localPoint = Point::pointTransformedBy(points[j], profile->getSection()->getGeometricalCenter(),
                                                      profile->getRotationAngle());
        double xj = localPoint->getLocalX();
        double yj = localPoint->getLocalY();
        delete localPoint;
        double eps_yj = eps_o - fi*yj;

        internalAction += Aj*std::pow(xj, r)*std::pow(yj, s)*material->sigmaEpsilon(eps_yj);

    }

    //printing for testing purposes
    if(DEBUG)
        fprintf(stderr, "Result of integration for this fiber group component is: %g\n", internalAction);

    return internalAction;
}

/**
 * @brief FiberGroup::yMinInCoordinateSystem
 * @param origin - origin of local coordinate system in which we are searching for y_MIN
 * @param rotationAngle  - omega angle of rotation of local coordinate system in relation to global central coordinate system of section
 * @return minimal y coordinate of this fiber group (component)
 */
double FiberGroup::yMinInCoordinateSystem(Point *origin, double rotationAngle)
{
    double y_min; //will store y_min

    if(points.length() < 1)
        throw NoVerticesInFiberGroupException("Brak wierzchołków zdefiniowanych dla elementu fiber group");

    Point *localPoint = Point::pointTransformedBy(points[0], origin, rotationAngle);
    y_min = localPoint->getLocalY();
    delete localPoint;

    for(int i=0; i< points.length(); i++) {
        Point *localPoint = Point::pointTransformedBy(points[i], origin, rotationAngle);
        if(localPoint->getLocalY() < y_min)
            y_min = localPoint->getLocalY();
        delete localPoint;
    }
    return y_min;
}

/**
 * @brief FiberGroup::yMaxInCoordinateSystem
 * @param origin - origin of local coordinate system in which we are searching for y_MAX
 * @param rotationAngle - omega angle of rotation of local coordinate system in relation to global central coordinate system of section
 * @return maximal y coordinate of this fiber group (component)
 */
double FiberGroup::yMaxInCoordinateSystem(Point *origin, double rotationAngle)
{
    double y_max; //will store y_max

    if(points.length() < 1)
        throw NoVerticesInFiberGroupException("Brak wierzchołków zdefiniowanych dla elementu fiber group");

    Point *localPoint = Point::pointTransformedBy(points[0], origin, rotationAngle);
    y_max = localPoint->getLocalY();
    delete localPoint;

    for(int i=0; i<points.length(); ++i) {
        Point *localPoint = Point::pointTransformedBy(points[i], origin, rotationAngle);
        if(localPoint->getLocalY() > y_max)
            y_max = localPoint->getLocalY();
        delete localPoint;
    }

    return y_max;
}
