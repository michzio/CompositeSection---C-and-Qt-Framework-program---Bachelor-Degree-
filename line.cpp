#include "line.h"
#include "strainprofile.h"
#include "gausslegendre.h"
#include <cmath>
#include "flags.h"


Line::Line(std::string lineName)
{
    name = lineName;
    //numberOfPoints = 0;
}
Line::Line(std::string lineName, Material * m) {
    name = lineName;
    material = m;
    //numberOfPoints = 0;
}

void Line::addPointAndArea(double x, double y, double A)
{
    addPoint(x,y);
    areas.append(A);
    //numberOfPoints++;
}

void Line::clearPointsAndAreas(void) {
    points.clear();
    areas.clear();
    //numberOfPoints = 0;
}

/*int Line::getNumberOfPoints(void)
{
    return numberOfPoints;
}*/

double *Line::getAreas()
{
    return areas.data();
}

std::string Line::getName()
{
    return name;
}

void Line::setName(std::string n)
{
    name = n;
}

void Line::setMaterial(Material *m)
{
    material = m;
}

Point *Line::getOrigin(void)
{
    if(points.length() >0)
        return points[0];
    return NULL;
}

Line::~Line() {
   // delete material;
}

Material *Line::getMaterial(void) {
   return material;
}


/****************************************************
 * getTotalArea() - this method calculate area      *
 * of Line element by adding each line component    *
 * together. It returns area of Line element as     *
 * double.                                          *
 * Area = SUM(areas[1] + ... + areas[n-1])            *
 ****************************************************/
double Line::getTotalArea(void) {

    double *areas = getAreas();
    double totalArea = 0;

    for(int i=0; i < (numberOfPoints() - 1); i++) {
        totalArea += areas[i];

    }

    return totalArea;
}

/****************************************************
 * getGeometricalCenterX() function calculates x_c  *
 * coordinate which is first coordinate of centroid *
 * of Line element.                                 *
 * This x_c is calculated as follows:               *
 * 1) For each component of Line element we calculat*
 *    its centroid  (xi_c, yi_c)                    *
 *    (xi_c, yi_c) = [(xi + xi_1)/2, (yi + yi_1)/2] *
 * 2) Based on (xi_c, yi_c) and Ai we calculate     *
 *    Syi, Sxi in default coordinate system         *
 * 3) Based on total Area of Line element A we get  *
 *    xc = SUM(0, n, Syi)/A; yc = SUM(0,n, Sxi)     *
 ****************************************************/
double Line::getGeometricalCenterX(void) {

    return firstMomentOfArea_Sy()/getTotalArea();
}

double Line::getGeometricalCenterY(void) {

    return firstMomentOfArea_Sx()/getTotalArea();
}

/******************************************************
 * firstMomentOfArea_Sx(void) is calculated in default*
 * coordinate system.                                 *
 * Sx = SUM(0,n, Ai * (yi_c - y) )                    *
 * yi_c - is y coordinate of centroid of ith component*
 *        of Line element                             *
 * y - is client defined coordinate system translation*
 *     in this case it is assumed as default (0,0)    *
 ******************************************************/
double Line::firstMomentOfArea_Sx(void) {
    return firstMomentOfArea_Sx(new Point(0.0,0.0));
}

double Line::firstMomentOfArea_Sx(Point *p) {
    double Sx = 0;
    double *areas = getAreas();

    //Calculating Sx = A0*(y0_c - y) + A1*(y1_c - y) + ... + An-1*(yn-1_c - y)
    for(int i=0; i<(numberOfPoints() -1); i++) {
        Sx += areas[i]*( getComponentCenterY(i) - p->getY());
    }

    return Sx;
}

/******************************************************
 * firstMomentOfArea_Sy(void) is calculated in default*
 * coordinate system.                                 *
 * Sy = SUM(0,n, Ai * (xi_c - x) )                    *
 * xi_c - is y coordinate of centroid of ith component*
 *        of Line element                             *
 * x - is client defined coordinate system translation*
 *     in this case it is assumed as default (0,0)    *
 ******************************************************/
double Line::firstMomentOfArea_Sy(void) {
    return firstMomentOfArea_Sy(new Point(0.0, 0.0));
}

double Line::firstMomentOfArea_Sy(Point *p) {
     double Sy = 0;
     double *areas = getAreas();

     //Calculating Sy = A0*(x0_c - x) + A1*(x1_c - x) + ... + An-1*(xn-1_c - x)
     for(int i=0; i<(numberOfPoints() -1); i++) {
         Sy += areas[i]*( getComponentCenterX(i) - p->getX());
     }

     return Sy;
}

Point *Line::getGeometricalCenter(void) {

    double x_c = getGeometricalCenterX();
    double y_c = getGeometricalCenterY();

    return new Point(x_c, y_c);
}

/****************************************************
 * getComponentCenter(int i) function calculates    *
 * centrid of i-th component of Line element        *
 * Formula: (xi_c, yi_c) = [(xi + xi_1)/2, (yi + yi_1)/2]
 * If i == n then i+1 = 0!                          *
 ****************************************************/
Point *Line::getComponentCenter(int i) {

    double xi_c = getComponentCenterX(i);
    double yi_c = getComponentCenterY(i);
    Point *centroid = new Point(xi_c, yi_c);

    return centroid;
}

//formula xi_c = (xi + xi_1)/2;  if(i==n) no line follow last vertex;
double Line::getComponentCenterX(int i) {
    double xi, xi_1;

    if(i + 1 >= numberOfPoints()) {
        fprintf(stderr, "In method Line::getComponentCenterX(int), argument 'i' is out of bounds!");
        throw AreaIndexOutOfBounds("In method Line::getComponentCenterX(int), argument 'i' is out of bounds!");
    }

    xi = points[i]->getX();
    xi_1 = points[i+1]->getX();

    return (xi + xi_1)/2;
}

//formula yi_c = (yi + yi_1)/2; if(i==n) no line follow last vertex;
double Line::getComponentCenterY(int i) {
    double yi, yi_1;

    if(i + 1 >= numberOfPoints()) {
        fprintf(stderr, "In method Line::getComponentY(int), argument 'i' is out of bounds!");
        throw AreaIndexOutOfBounds("In method Line::getComponentCenterY(int), argument 'i' is out of bounds!");
    }

    yi = points[i]->getY();
    yi_1 = points[i+1]->getY();

    return (yi + yi_1)/2;
}


/*********************************************************************/
/* STRESS INTEGRATION METHODS                                        */
/*********************************************************************/

double Line::internalAction_N(StrainProfile *profile)
{
    return performStressIntegration(profile, ACTION_N);
}

double Line::internalAction_Mx(StrainProfile *profile)
{
    return performStressIntegration(profile, ACTION_Mx);
}

double Line::internalAction_My(StrainProfile *profile)
{
    return -performStressIntegration(profile, ACTION_My);
}

/**
  * @brief Line::performStressIntegration
  * @param profile - StrainProfile for which to calculate given internal action (N, Mx, My)
  *                  for this surface
  * @param actionType - type of calculated internal action: N, Mx, My
  * @return - internal action (N, Mx, My) value as double
  * To perform stress integration we apply numerical integration scheme:
  * the Gauss-Legendre quadrature method. We use adaptive strain-mapped integration method.
  * Resultant action Rl,i contributing to the total section response can be calculated using line integral fromula:
  * Rl,i = Sl,i(x^r*y^s*t(y)*sigma(y)dy)
  */
 double Line::performStressIntegration(StrainProfile *profile,
                                          InternalAction actionType)
 {
     double internalAction = 0.0;

     //exponents which value depends on resultant internal action:
     // (r,s) = (0,0) for N, (r,s) = (0,1) for Mx, (r,s) = (1,0) for -My
     int r = 0, s = 0;

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


     //Line integral along the segments of line element Li of nli segments
     //can be expressed as sum of integrals of each individual segment lj
     //Rl,i = SL,i x^r*y^s*t(y)*sigma(y)dy -> SUM(j=1, nli, tj * Sl,j x^r*y^s*sigma(y)*dy)
     //where tj * Sl,j x^r*y^s*sigma(y)*dy is Ij  (elementary integral for each individual segment)

     //loop through each line element segment
     //and adding up resultant action contribution fromm each integrated segment (Ij)
     //SUM(j=1, nli, tj *Sl,j x^r*y^s*sigma(y)*dy)
     double xi,yi, xi_1, yi_1;
     int numOfPoints = this->numberOfPoints();
     Point *sectionCentroid = profile->getSection()->getGeometricalCenter();
     double omega = profile->getRotationAngle();

     for(int i=0; i< numOfPoints; i++) {

         //first line element's segment end coordinates
         xi = points[i]->getX();
         yi = points[i]->getY();

        //second line element's segment end coordinates
        if( i+1 < numOfPoints) {
            xi_1 = points[i+1]->getX();
            yi_1 = points[i+1]->getY();
        } else {
            xi_1 = points[0]->getX();
            yi_1 = points[0]->getY();
        }

        //creating segment object lj with end's points transformed to local strain profile coordinate system
        Point pi(xi,yi), pi_1(xi_1, yi_1);
        pi.setOriginAndRotation(sectionCentroid->getX(),sectionCentroid->getY(), omega);
        pi_1.setOriginAndRotation(sectionCentroid->getX(), sectionCentroid->getY(), omega);
        Segment lj(pi,pi_1);

        //adding up each elementary integral Ij for each line element's segment
        //we pass current segment object lj, strain profile and
        //exponents r and s which depends on calculated internal action N, Mx, My
        internalAction += segmentIntegral_Ij(lj, profile, r, s, i);
     }

     return internalAction;
}

 /**
 * @brief Line::segmentIntegral_Ij
 * @param lj - jth line element's segment on which stress integration is performed
 * @param profile - strain profile for the current line element for which we get
 *                  stress values to integrate sigma(y)
 * @param r - exponent for x coordinate
 * @param s - exponent for y coordinate
 *        r,s exponent are used to differentiate between N, Mx, My calculation
 * @param lineIdx - is index of line element's segment used in order to get area of this segment or thickness
 * @return resulting internal action from line element's segment integration
 */
double Line::segmentIntegral_Ij(const Segment& lj, StrainProfile* profile,
                          double r, double s, int lineIdx)
{
        double internalAction = 0.0; //elementary Integral Ij for line element's segment

        //elementary integral of each segment Ij is expressed as
        //Ij = tj*Sl,j (aj*y + bj)^r*y^s*sigma(y)dy
        //where x has been substituted by x = aj*y + bj
        //we can get this from y = aj'*x + bj' -> aj'*x = y - bj' -> x = 1/aj' *y - bj'/aj'
        //aj = 1/aj' and bj = - bj'/aj'
        //where:
        //aj' = (yi_1 - yi)/(xi_1 - xi)
        //bj' =  yi - aj'*xi or bj' = yi_1 - aj'*xi_1
        //we can also directly calculate aj and bj in x = aj*y + bj
        //assuming that bj intercept x axis and traditional x axis is now y and y is now x
        //aj = (xi_1 - xi)/(yi_1 - yi)   which is in accordance with above indirect formulation
        //bj = xi - aj*yi

        double xi = lj.getStartPoint().getLocalX();
        double yi = lj.getStartPoint().getLocalY();
        double xi_1 = lj.getEndPoint().getLocalX();
        double yi_1 = lj.getEndPoint().getLocalY();

        double aj = (xi_1 - xi)/(yi_1 - yi);
        double bj = xi - aj*yi;

        //under integral we have nonlinear function of single variable y
        //we can numericaly calculate this linear integral using Gausse quadrature
        //we use formula:
        //Sl,j Fj(y)dy = SUM(k=1, nf, 1/2*(y_lim_Bk - y_lim_Ak)*SUM(m=1, nGk, wm*Fj(ym)))
        //where nf is number of material sigma-epsilon function parts
        //y_lim_Bk, y_lim_Ak are strain-to-coordinate mapped limits of function parts
        //nGk is number of Gauss points for kth function part
        //wm, ym depends on number of Gauss points
        //Fj(y) is function Fj(y) = (ajy + bj)^r * y^s *sigma(y)

        //we can see that in elementary integral this first integral is preceeded by tj - thickness
        //taking it into accoun we get tj*1/2*(y_lim_Bk - y_lim_Ak) == Ak/2
        //where Ak can be also expressed as Ak = Aj * |(y_lim_Bk - y_lim_Ak)/(yBj = yAj)|
        //where Aj is area of jth segement of line element Li

        double eps_o = profile->getStrainAtOrigin_eps_o();
        double fi = profile->getCurvature_fi();
        //jth segment area
        double Aj = areas[lineIdx];

        //IMPORTANT! contrary to surface segment integration in case of line elements
        //horizontal segments are now treated specially as they contribute to total area

        if( std::abs(yi - yi_1) < 0.000001)
        {

            //stress in all jth segment
            double eps_yi = eps_o - fi*yi;
            double sig_yj = material->sigmaEpsilon(eps_yi);

            //SEGMENT IS HORIZONTAL
            if(r == 0 && s==0) {
                //internal action calculated is N
                double N = Aj*sig_yj;
                internalAction = N;

            } else if(r == 0 && s == 1)
            {   //internal action calculated is Mx
                double Mx = Aj*yi*sig_yj;
                internalAction = Mx;
            } else if(r == 1 && s==0)
            {
                //internal action calculated is -My
                double negativeMy = Aj*((xi + xi_1)/2)*sig_yj;
                internalAction = negativeMy;
            }

            return internalAction;
        }

        //all remaining cases
        //1. we need segement y_Aj, y_Bj which are y coordinate segment min and max limits
        double y_Aj = std::min(yi, yi_1);
        double y_Bj = std::max(yi, yi_1);

        //2. loop through each material constitutive law function parts
        //   in order to integrate each contributing part
        for(int k=0; k < material->numberOfFunctionParts(); k++)
        {
            if(DEBUG)
                 fprintf(stderr, "Line %d, Material definition function part: %d.\n", lineIdx, k);
            //reading strain limits for each function part in material definition
            double eps_limA_k = material->functionPartStrainLimitA(k); //is more tensile
            double eps_limB_k = material->functionPartStrainLimitB(k); //is more compressive
            //calculating y_limA_k and y_limB_k coordinate limits for given eps_limA_k and eps_limB_k values
            //by simple strain-to-coordinate transformation
            double y_limA_k = 0.0;
            double y_limB_k = 0.0;
            if(fi == 0) {
                //we have strain profile curvature equal to 0
                //so strains are constant on height of strain profile (section)
                //then we have workround:
                // if eps_o is contained in [eps_limA_k, eps_limB_k]
                // we assume y_limA_k = y_Aj and y_limB_k = y_Bj
                // so eps_o <= eps_limA_k than y_limA_k = y_Aj
                // eps_o >= eps_limB_k than y_limB_k = y_Bj
                if(eps_o <= eps_limA_k && eps_o >= eps_limB_k)
                {
                    //y_limA_k = y_Aj;
                    //y_limB_k = y_Bj; //y_Aj always< y_Bj isn't here correct?
                    //we should rather take this values in accordance with segments traversal direction
                    y_limA_k = yi;
                    y_limB_k = yi_1;

                }
            } else {
                if( (yi_1 - yi) > 0) {  //modification 17 Dec 2013 (test)
                    y_limA_k = (eps_o - eps_limA_k)/fi;
                    y_limB_k = (eps_o - eps_limB_k)/fi;
                } else {
                    y_limA_k = (eps_o - eps_limB_k)/fi;
                    y_limB_k = (eps_o - eps_limA_k)/fi;
                }
            }

            //where eps0 is strain-at-origin and fi is strain profile curvature
            //!!!! PROBLEM WHETHER y_limA_k should be < y_limB_k always to get !!!!
            //!!!! CORRECT INTERNAL ACTION CONTRIBIUTION OR IT CAN BE IN ANY ORDER? !!!!
            //testing obtained y_limA_k, y_limB_k coordinate limits whether thay fall outside
            //a segment limit (y_Aj, y_Bj) if so they are replaced by the segment limit itself
            if(y_limA_k < y_Aj)
                y_limA_k = y_Aj; //when y_limA_k below segment y_min
            if(y_limB_k < y_Aj)
                y_limB_k = y_Aj; //when y_limB_k below segment y_min
            if(y_limA_k > y_Bj)
                y_limA_k = y_Bj; //when y_limA_k above segment y_max
            if(y_limB_k > y_Bj)
                y_limB_k = y_Bj; //when y_limB_k above segment y_max

            //if both y_limA_k and y_limB_k coordinate limits fall outside
            //the same segment limit they have gotten the same y_Aj or y_Bj segment limit
            if(y_limA_k == y_limB_k)
                continue; //we skip integration of this function part f,k as its corresponding
                          //strain range fall outside the segment

            //now we have to integrate fragment of segment between y_limA_k and y_limB_k
            //using function part f,k to calculate stresses in corresponging points
            //based on y contained in [y_limA_k, y_limB_k] y -> epsilon -> f,k(epsilon) -> stress to integrate
            //now Gaussian sampling is applied for all contributing function parts f,k of the considered segment lj

            //Ij_k - elementary integral kth component
            //Ak area of segment fragment which is within limits [y_limA_k, y_limB_k]
            double Ak = Aj*std::abs( (y_limB_k - y_limA_k)/(y_Bj - y_Aj) );
            double Ij_k = Ak/2;

            //calculation of elementary integral Ij_k kth component depends on number of utilized Gauss points nG (quadrature order)
            int fkOrder = material->functionPartDegree(k);
            int nG = (fkOrder + r + s)/2 + 1; //number of required Gaussian points to get exact results for polynomial integration of order fkOrder!

            double summation = 0.0;
            //we generate gauss object from which we will be obtaining weights w_m and y_m coordinates
            GaussLegendre gauss(nG);
            for(int m = 0; m < nG; ++m)
            {
                double ym = 0.5*(y_limB_k + y_limA_k) + gauss.lambda(m)*0.5*(y_limB_k - y_limA_k);
                double eps_ym = eps_o - fi*ym;
                double Fj_ym = std::pow(aj*ym + bj, r)*std::pow(ym, s)*material->sigmaEpsilon(eps_ym);
                //adding up consecutive parts of Gaussian quadrature
                summation += gauss.weight(m)*Fj_ym;
            }

            Ij_k *= summation;

            //cumulating elementary integral components in internalAction variable
            internalAction += Ij_k;

        }

        if(DEBUG)
             fprintf(stderr, "Wartość całkowania dla j-tego segmentu Ij = %g.\n", internalAction);

        return internalAction;
}

/**
 * @brief Line::yMinInCoordinateSystem
 * @param origin - origin of local coordinate system in which we are searching for y_MIN
 * @param rotationAngle - omega angle of rotation of local coordinate system in relation to global central coordinate system of section
 * @return minimal y coordinate of this line element (component)
 */
double Line::yMinInCoordinateSystem(Point *origin, double rotationAngle)
{
    double y_min; //will store y_min

    if(points.length() < 1)
        throw NoVerticesInLineException("Brak wierzchołków zdefiniowanych dla elementu liniowego");

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
 * @brief Line::yMaxInCoordinateSystem
 * @param origin - origin of local coordinate system in which we are searching for y_MAX
 * @param rotationAngle - omega angle of rotation of local coordinate system in relation to global central coordinate system of section
 * @return maximal y coordinate of this line element (component)
 */
double Line::yMaxInCoordinateSystem(Point *origin, double rotationAngle)
{

     double y_max; //will store y_max

     if(points.length() < 1)
         throw NoVerticesInLineException("Brak wierzchołków zdefiniowanych dla elementu liniowego");

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
