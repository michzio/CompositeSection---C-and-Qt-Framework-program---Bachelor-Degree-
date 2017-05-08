#include "surface.h"
#include "polygonintersection.h"
#include "strainprofile.h"
#include "flags.h"
#include "gausslegendre.h"
#include <memory>

Surface::Surface(std::string surfaceName) : name(surfaceName)
{
}

Surface::Surface(std::string surfaceName, Material *m)  {
    name = surfaceName;
    material = m;
}

//copy constructor
Surface::Surface(Surface &other) {

    name = other.getName();
    //copy material definition
    Point **points = other.getPointsArray();
    double *angles = other.getAngels();
    for(int i=0; i < other.numberOfPoints(); i++)
    {
       addPointAndAngel(points[i]->getX(), points[i]->getY(), angles[i]);
    }
}

Surface::~Surface()
{
   //delete material;
}

void Surface::setName(std::string n) {
    name = n;
}

void Surface::setIsAdaptiveLinearized(bool flag)
{
    isAdaptiveLinearized = flag;
}

void Surface::addPointAndAngel(double x, double y, double a) {
    addPoint(x,y);
    angles.append(a);
    //verticesNumber++;
}

void Surface::clearPointsAndAngles()
{
    this->points.clear();
    this->angles.clear();
    //verticesNumber = 0;
}

double *Surface::getAngels()
{
    return angles.data();
}

void Surface::setMaterial(Material * m)
{
   material = m;
}

void Surface::setIsOpening(bool flag) {
    isOpening = flag;
}

bool Surface::getIsOpening() {
    return isOpening;
}

std::string Surface::getName() {
    return name;
}

/*int Surface::getVerticesNumber() {
    return verticesNumber;
}*/

Material *Surface::getMaterial()
{
  return material;
}

void Surface::performAdaptiveLinearization()
{
    //area of the current polygon surface
    double Atot_sweep1 = this->calculatePolygonAreaIncludingArcs();
    double Atot_sweep2 = this->calculatePolygonArea();

    while ( std::abs((Atot_sweep2 - Atot_sweep1)/Atot_sweep1) > A_TOL)  {

        int pos; //position of longest arc, ex. arc between vertices 0-1 is 0, 5-6 is 5 etc.
        int next_pos;

        //select longest arc
        Arc *arc = getArcWithLargestSegmentArea(&pos);
        fprintf(stderr, "New longest arc (%d) starts at: (%g, %g) and end at: (%g, %g).\n", pos, arc->getStartPoint()->getX(),
                arc->getStartPoint()->getY(), arc->getEndPoint()->getX(), arc->getEndPoint()->getY());
        //divide this arc equaly by inserting arc midpoint
        Point *arcMidpoint = this->arcMidPoint(arc->getStartPoint(), arc->getEndPoint(), arc->getAngle());
        fprintf(stderr, "New midpoint has coordinates: (%g, %g).\n", arcMidpoint->getX(), arcMidpoint->getY());
        //insert obtained arc midpoint at suitable position
        next_pos = (pos+1); //%numberOfPoints();
        this->points.insert(next_pos, arcMidpoint);
        this->angles.insert(next_pos, arc->getAngle()/2);
        this->angles.replace(pos, arc->getAngle()/2);
        //calculate polygon area with new point inserted
        //Atot_sweep1 = Atot_sweep2;
        Atot_sweep2 = this->calculatePolygonArea();
    }

    if(DEBUG)
        fprintf(stderr, "Approximated area: %g.",Atot_sweep2);
    isAdaptiveLinearized = true;
    drawSurface(material->materialColor());

}


/****************************************************
 * The area of polygon is the measurement of the 2D *
 * region enclosed by the polygon. The formula for  *
 * non-self-interecting (simple) polygon with n     *
 * vertices is:                                     *
 * A = 1/2*SUM(i=0, n-1, xi*y(i+1) - x(i+1)yi))     *
 ****************************************************/
double Surface::calculatePolygonArea() {

    Point **points = this->getPointsArray();
    int numOfPoints = this->numberOfPoints();
    double area = 0.0;
    double sumOfDeterminants = 0.0;
    double xi = 0.0;
    double yi = 0.0;
    double yi_1 = 0.0;
    double xi_1 = 0.0;

    for(int i=0; i<numOfPoints; i++) {
         xi = points[i]->getX();
         yi = points[i]->getY();

        if( i+1 < numOfPoints) {
            xi_1 = points[i+1]->getX();
            yi_1 = points[i+1]->getY();
        } else {
            xi_1 = points[0]->getX();
            yi_1 = points[0]->getY();
        }

        sumOfDeterminants += (xi*yi_1 - xi_1*yi);
    }
    area = 0.5*sumOfDeterminants;

    fprintf(stderr, "Area: %g.\n", area);
    return area;
}

/*******************************************************
 * This funcion is similar to calculatePolygonArea     *
 * However it also takes into account circular arc     *
 * segments. It subtructs arc area if angle is negative*
 *and it adds arc area if angle is positive.           *
 *******************************************************/
double Surface::calculatePolygonAreaIncludingArcs() {

    Point **points = this->getPointsArray();
    int numOfPoints = this->numberOfPoints();
    double *angles = this->getAngels();
    double area = 0.0;

    area = this->calculatePolygonArea();
    fprintf(stderr, "Number of vertices: %d.\n", numOfPoints);

    for(int i=0, j=0; i<numOfPoints; i++) {
        j = ((i+1) == numOfPoints ) ? 0 : (i+1);
        fprintf(stderr, "Angel at side: %d-%d has value: %g.\n", i, j, angles[i]);
        if(angles[i] > 0) {
            area += arcSegmentArea(points[i], points[j], angles[i]);
        } else if(angles[i] < 0) {
            //arcSegmentArea if radius < 0 is also <0 so to
            //subtruct segment area we can just add it
            area += arcSegmentArea(points[i], points[j], angles[i]);
        }
    }

    fprintf(stderr, "Polygon area including arcs: %g.\n", area);

    return area;
}

/****************************************************
 * The area of the shape limited by the arc and a   *
 * straight line between the two end points is:     *
 * 1/2*r^2*(0 - sin(0))                             *
 * http://en.wikipedia.org/wiki/Arc_(geometry)#Arc_segment_area *
 ****************************************************/
double Surface::arcSegmentArea(const Point* start, const Point* end, double angle) {
      double area = 0.0;
      double lengthOfChord = lengthOfVector(start,end);

      if(DEBUG)
      fprintf(stderr, "Length of chord between point: (%g, %g) and point: (%g, %g) is: %g.\n",
              start->getX(), start->getY(), end->getX(), end->getY(), lengthOfChord);

      double radius = radiusOfArc(lengthOfChord, angle);
      if(DEBUG)
      fprintf(stderr, "Radius: %g.\n", radius);

      area = 0.5*radius*radius*(angle - std::sin(angle));

      if(DEBUG)
      fprintf(stderr, "Arc segment has area: %g.\n", area);

      return area;
}

/*****************************************************
 * This method calculates the length of vector for   *
 * given starting and ending point (coordinates)     *
 * The formula: sqrt(pow(x1-x0) + pow(y1-y0))        *
 *****************************************************/
double Surface::lengthOfVector(const Point* start, const Point* end) {
    double length;
    double delX = end->getX() - start->getX();
    double delY = end->getY() - start->getY();
    length = std::sqrt(std::pow(delX,2) + std::pow(delY,2));
    return length;
}

/*****************************************************
 * This method calculates the radius of arc for given*
 * distance between two vercities and angle of arc   *
 * The Formula is: r = W/(2*sin(0/2))  where         *
 * W - is distance between start and end point and   *
 * 0 - is omega angle of arc                         *
 *****************************************************/
double Surface::radiusOfArc(double lengthOfChord, double angle)
{
    double radius = lengthOfChord/(2*std::sin(angle/2));
    return radius;
}

/*****************************************************
 * Find middle point (xm, ym) on the line between    *
 * two points (x0, y0) and (x1, y1).                 *
 * The Formula: xm = (x0 + x1)/2, ym = (y0 + y1)/2   *
 *****************************************************/
Point *Surface::middlePointOfSegment(const Point *p0, const Point *p1)
{
   Point *middlePoint = new Point();
   middlePoint->setX((p0->getX() + p1->getX())/2);
   middlePoint->setY((p0->getY() + p1->getY())/2);
   return middlePoint;
}


/*****************************************************
 * Center of circle given two points and radius      *
 * 1) (x0, y0) - start point (x1, y1) - end point    *
 *    of arc (point on the circle)                   *
 *    (xm, ym) - middle point of the chord           *
 *    (xc, yc) - center of the circle                *
 * 2) distance between point (x1,y1) and (xc,yc) is  *
 *    sqrt((xc-x1)^2 + (yc-y1)^2)                    *
 * 3) we have set of equations:                      *
 *    (xc-x0)^2 + (yc-y0)^2 = r^2                    *
 *    (xc-x1)^2 + (yc-y1)^2 = r^2                    *
 * 4) or simpler method using "vector theorem":      *
 * 5) find middle point (xm,ym) coordinates          *
 * 6) find the direction of the line between start   *
 *    point 0 and end point 1: [x1-x0, y1-y0]        *
 * 7) find the direction of "mirror line" which is   *
 *    perpendicular to line between point 0 and 1    *
 *    ex. [x1-x0, y1-y0] ---> [y1-y0, x0-x1]         *
 * 8) next we normalize direction vector of mirror   *
 *    line. v/||v|| where ||v|| = sqrt(x^2 + y^2)    *
 * 9) the two circle centers will both be on         *
 *    the mirror line, and we can use geometry to    *
 *    find how far they are from the point (xm,ym)   *
 *10) distance to move along the mirror line is:     *
 *    sqrt(r^2-(d/2)^2) where r - radius,            *
 *    d - distance from start point to end point     *
 *11) we can move up or down the mirror line using   *
 *    its direction vector multiplied by 10)         *
 * One answer will be:                               *
 * x = xm + sqrt(r^2-(d/2)^2)*(y1-y0)/q              *
 * y = ym + sqrt(r^2-(d/2)^2)*(x0-x1)/q              *
 * Second answer will be:                            *
 * x = xm - sqrt(r^2-(d/2)^2)*(y1-y0)/q              *
 * y = ym -  sqrt(r^2-(d/2)^2)*(x0-x1)/q             *
 * Finally:                                          *
 * function returns two points of possible centers   *
 *****************************************************/
Point **Surface::centersOfCircleGiven2PointsAndRadius(const Point *start, const Point *end, double radius)
{
        Point *middlePoint = this->middlePointOfSegment(start, end);
        Vector *segmentDirectionVector = new Vector(const_cast<Point*>(start), const_cast<Point*>(end));
        double d = segmentDirectionVector->normOfVector();
        fprintf(stderr, "Norm of segment direction vector: %g.\n", d);
        Vector *mirrorDirectionVector = segmentDirectionVector->perpendicularVector();
        mirrorDirectionVector->normalizeVector();

        double distanceFromMiddlePointToCenter = std::sqrt(radius*radius - (d/2)*(d/2));
        fprintf(stderr, "Distance From Middle Point (%g, %g) to Center: %g.\n",
                middlePoint->getX(), middlePoint->getY(), distanceFromMiddlePointToCenter);
        double dx = distanceFromMiddlePointToCenter * mirrorDirectionVector->getX();
        double dy = distanceFromMiddlePointToCenter * mirrorDirectionVector->getY();

        Point *center1 = new Point(middlePoint->getX() + dx, middlePoint->getY() + dy);
        Point *center2 = new Point(middlePoint->getX() - dx, middlePoint->getY() - dy);
        fprintf(stderr, "Center 1 coordinates: (%g, %g), Center 2 coordinates(%g, %g) \n",
                center1->getX(), center1->getY(), center2->getX(), center2->getY());

        Point **centers;
        centers = new Point*[2];
        centers[0] = center1;
        centers[1] = center2;

        return centers;
}

/*****************************************************
 * Method selecting proper center of circle based on *
 * sign (+/-) of arc's angle and whether point       *
 * lies on the left or right of the segment i.e.     *
 *  - sign means subtructing arc segment area so     *
 *     center of arc (circle) should on the right of *
 *     segment or on the segment                     *
 * + sign means adding arc segment area so center of *
 *   arc (circle) should on the left of segment or   *
 *   on the segment                                  *
 *****************************************************/
 Point *Surface::selectProperCenterOfCircle(Point **centers, Point *start, Point *end, double angle)
 {
     if(angle > 0) {
        //arc segment area is added to polygon area
        //center of circle (arc) should be to the left of segment or on the segment
         if(!centers[0]->liesToTheRightOfSegment(*start, *end))
             return centers[0];
         else if(!centers[1]->liesToTheRightOfSegment(*start, *end))
             return centers[1];
         else
             return NULL;
     } else if(angle < 0) {
         //arc segment area is subtructed from polygon area
         //center of circle (arc) should be to the right of segment or on the segment
         if(!centers[0]->liesToTheLeftOfSegment(*start, *end))
             return centers[0];
         else if(!centers[1]->liesToTheLeftOfSegment(*start, *end))
             return centers[1];
         else
             return NULL;
     } else {
         //case when the segment was straight line
        return NULL;
     }
 }

 /****************************************************
  * Determining if Point lies inside the polygon.    *
  * The solution is to compute the sum of angles     *
  * mage between the test point and each pair of     *
  * points making up the polygon. If this sum is     *
  * 2*PI then the point is an interior point else if *
  * it is 0 then the point is an exterior point.     *
  * http://bbs.dartmouth.edu/~fangq/MATH/download/source/Determining%20if%20a%20point%20lies%20on%20the%20interior%20of%20a%20polygon.htm
  ****************************************************/
 bool Surface::pointIsInsieThePolygon(const Point *p)
 {
     Vector v1,v2;
     double x = 0.0, y = 0.0;
     double angle = 0.0;

     Point **points = getPointsArray();

     for(int i=0; i<numberOfPoints(); i++) {

         v1.setStartPoint(new Point(p->getX(),p->getY()));
         x = points[i]->getX();
         y = points[i]->getY();
         v1.setEndPoint(new Point(x,y));

         v2.setStartPoint(new Point(p->getX(), p->getY()));
         //remainder of division was used to change last counter N to 0
         //when we calculate angle between vertices N-1 and 0
         x = points[(i+1)%numberOfPoints()]->getX();
         y = points[(i+1)%numberOfPoints()]->getY();
         v2.setEndPoint(new Point(x,y));

         angle += calculateAngle2D(&v1, &v2);
     }

     if(std::abs(angle) < M_PI) return false;
     else
     return true;

 }

 /****************************************************
  * This method calculates the angle between two     *
  * vectors on a plane. The angle is from v1 to v2,  *
  * positive anticlockwise, between (-PI, PI)        *
  * When assesing outside/inside position of point P *
  * v1 - vector from P to first vertex of polygon    *
  * v2 - vector from P to second vertex of polygon   *
  ****************************************************/
 double Surface::calculateAngle2D(Vector *v1, Vector *v2)
 {
       double alfa_diff, alfa1, alfa2;

       alfa1 = std::atan2(v1->getY(), v1->getX());
       alfa2 = std::atan2(v2->getY(), v2->getX());
       alfa_diff = alfa2 - alfa1;

       while(alfa_diff > M_PI)
           alfa_diff -= 2*M_PI;
       while(alfa_diff < -M_PI)
           alfa_diff += 2*M_PI;

       return (alfa_diff);
 }

/*******************************************************
 * Middle Point of Arc can be calculated using dot     *
 * product. Formula for coordinates of arc midpoint is *
 * C - center, A - startPoint, B - end point, M - mid  *
 * point. Coordinates of each point are:               *
 * c, 2a, 2b, m respectively, angle is 2*Alfa  then    *
 * m = c + (a + b - c)*sec(Alfa)                       *
 *******************************************************/
 Point *Surface::arcMidPoint(Point *start, Point *end, double angle)
 {
      //calculate coordinates of Center of Circle (arc)
     double lengthOfChord = this->lengthOfVector(start, end);
     double radius = this->radiusOfArc(lengthOfChord, angle);
     Point **centers = this->centersOfCircleGiven2PointsAndRadius(start, end, radius);
     //select proper center between two possible
     Point *center = this->selectProperCenterOfCircle(centers, start, end, angle);
     fprintf(stderr, "Center of circle has coordinates: (%g, %g).\n", center->getX(), center->getY());

             //calculate coordinates of arc midpoint
     double xm, ym, xa, ya, xb, yb, xc, yc, alfa;

     xc = center->getX(); yc = center->getY();
     xa = start->getX()/2; ya = start->getY()/2;
     xb = end->getX()/2; yb = end->getY()/2;
     alfa = angle/2;

     if(angle == M_PI)  {
         if(xa == xb) { //180degree, vertical direction
               ym = yc + (ya + yb - yc)*secans(alfa);
               xm = xa + (yb - ya);
         } else if(ya == yb) { //180degree, horizontal direction
               xm = xc + (xa + xb - xc)*secans(alfa);
               ym = ya - (xb-xa);
         }
     } else if(angle == -M_PI) {
         if(xa == xb) { //-180degree, vertical direction
               ym = yc + (ya + yb - yc)*secans(alfa);
               xm = xa - (yb - ya);
         } else if(ya == yb) { //-180degree, horizontal direction
               xm = xc + (xa + xb - xc)*secans(alfa);
               ym = ya + (xb-xa);
         }
     } else {
        xm = xc + (xa + xb - xc)*secans(alfa);
        ym = yc + (ya + yb - yc)*secans(alfa);
    }

     return new Point(xm, ym);
 }

 double Surface::secans(double x) {
     return 1/std::cos(x);
 }

 Arc *Surface::getArcWithLargestSegmentArea(int *pos)
 {
     Point **points = this->getPointsArray();
     double distance, longest = 0.0, largestAngle = 0.0;
     double largestArcSegmentArea = 0.0;
     int currPos = -1;

     for(int i = 0; i < numberOfPoints(); i++) {
         double currentArcSegmentArea = arcSegmentArea(points[i], points[(i+1)%numberOfPoints()], angles[i]);
         if( currentArcSegmentArea > largestArcSegmentArea)
         {
             largestArcSegmentArea = currentArcSegmentArea;
             currPos = i;
         }

         /* OLD VERSION OF FINDING LONGEST ARC
           if(angles[i] != 0 && angles[i] > largestAngle) {
             largestAngle = angles[i];
             distance = this->lengthOfVector(points[i], points[(i+1)%numberOfPoints()]);
             if(distance > longest) {
                 longest = distance;
                 currPos = i;
             }
         }*/
     }

     if(currPos > -1) {
         *pos = currPos;
         return new Arc(points[currPos], points[(currPos+1)%numberOfPoints()], angles[currPos]);
     }

     return NULL;
 }


 /***************************************************************
  * getGeometricalCenterX(), getGeometricalCenterY() methods    *
  * return x_c and y_c coordinate of geometrical center of      *
  * current surface. This Geometrical center is calculated      *
  * based on vertices of polygon (after adaptive linearization):*
  * Formula:                                                    *
  *  x_c = 1/6A * SUM(i=0 to n-1, (x[i] + x[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
  *  y_c = 1/6A * SUM(i=0 to n-1, (y[i] + y[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
  * Source of formula: http://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
  ***************************************************************/
 double Surface::getGeometricalCenterX(void)
 {
        double x_c = 0;
        double xi, yi, xi_1, yi_1;
        int numOfPoints = this->numberOfPoints();

        for(int i=0; i< numOfPoints; i++) {

            xi = points[i]->getX();
            yi = points[i]->getY();

           if( i+1 < numOfPoints) {
               xi_1 = points[i+1]->getX();
               yi_1 = points[i+1]->getY();
           } else {
               xi_1 = points[0]->getX();
               yi_1 = points[0]->getY();
           }

           x_c +=(xi + xi_1) * (xi*yi_1 - xi_1*yi);
        }

        x_c /= 6* this->calculatePolygonArea();

        return x_c;
 }

 double Surface::getGeometricalCenterY(void)
 {
     double y_c = 0;
     double xi, yi, xi_1, yi_1;
     int numOfPoints = this->numberOfPoints();

     for(int i=0; i< numOfPoints; i++) {

         xi = points[i]->getX();
         yi = points[i]->getY();

        if( i+1 < numOfPoints) {
            xi_1 = points[i+1]->getX();
            yi_1 = points[i+1]->getY();
        } else {
            xi_1 = points[0]->getX();
            yi_1 = points[0]->getY();
        }

        y_c +=(yi + yi_1) * (xi*yi_1 - xi_1*yi);
     }

     y_c /= 6* this->calculatePolygonArea();

     return y_c;

 }
//this version of function returns geometrical center as Point *object!
 Point *Surface::getGeometricalCenter(void)
 {
     if(this->centroid == NULL) {
        double x_c, y_c;
        x_c = getGeometricalCenterX();
         y_c = getGeometricalCenterY();
         Point *point = new Point(x_c, y_c);
         centroid = point;
         return point;
      }
     return this->centroid;
 }

 /*****************************************************
  * firstMomentOfArea_Sx() and firstMomentOfArea_Sy() *
  * is calculating this values based on user indicated*
  * point of reference and geometrical center of      *
  * this surface by multiplying (y_c - y)*A and       *
  * (x_c - x)*A.                                      *
  *****************************************************/
 double Surface::firstMomentOfArea_Sx(void)
 {
     //we use default reference point (0,0)
     //in which user has defined surface
     Point *referencePoint = new Point(0,0);
     return firstMomentOfArea_Sx(referencePoint);
 }

 double Surface::firstMomentOfArea_Sx(Point *p)
 {
     double Sx = 0;
     if(!isAdaptiveLinearized) performAdaptiveLinearization();
     //Calculating Sx = (yc - y) * A, where yc is y coordinate of surface centroid
     Sx = (getGeometricalCenterY() - p->getY()) * calculatePolygonArea();
     //in above formula we use polygon area after adaptive linearization

     if(isOpening)  Sx = 0; //if surface is opening it won't contribute to Sx

     return Sx;
 }

 double Surface::firstMomentOfArea_Sy(void)
 {
     //we use default reference point (0,0)
     //in which user has defined surface
     Point *referencePoint = new Point(0,0);
     return firstMomentOfArea_Sy(referencePoint);
 }

 double Surface::firstMomentOfArea_Sy(Point *p)
 {
    double Sy = 0;

    if(!isAdaptiveLinearized) performAdaptiveLinearization();
    //Calculating Sy = (xc - x) * A where xc is x coordinate of surface centroid
    Sy = (getGeometricalCenterX() - p->getX()) * calculatePolygonArea();
    //in above formula we use polygon area after adaptive linearization

    if(isOpening) Sy = 0; //if surface is opening it won't contribute

    return Sy;
 }

 /******************************************************
  * intersection(Surface *) - function finds and return*
  * Surface which is intersection of two polygons. One *
  * polygon representing THIS surface and polygon which*
  * represents Surface passed as argument.             *
  * Functin uses API gpc (General Polygon Clipping)    *
  * to do this work. This API is wrraped by helper     *
  * PolygonIntersection class.                         *
  ******************************************************/
 QVarLengthArray<Surface *> *Surface::intersection(Surface *otherSurface)
 {
     //using helper class to get surfaces as a result
     //of intersecting two surfaces (polygons)
     return PolygonIntersection::intersect(this, otherSurface);
 }

 bool Surface::hasBeenAdaptiveLinearized(void)
 {
     return isAdaptiveLinearized;
 }

 /**
  * @brief yMinInCoordinateSystem
  * @param origin - origin of coordinate system (usually section centroid)
  * @param rotationAngle - angle by which coordinate system has been rotated
  * Method finds y_min for current surface. y_MIN = min{y_min1, y_min2, ... y_minn}
  * will be the coordinate in local coordinate system (central coordinate system rotated by angle 0 (OMEGA))
  * where we will have bottom boundary of STRAIN PROFILE ex. eps_tu in this coordinate i.e. ultimate tensile strength
  * @return
  */
 double Surface::yMinInCoordinateSystem(Point *origin, double rotationAngle)
 {
        double y_min; //will store y_min

        if(!isAdaptiveLinearized)
            performAdaptiveLinearization();

        if(points.length() <1)
            throw NoVerticesInSurfaceException("Brak wierzchołków zdefiniowanych dla powierzchni");

        Point *localPoint = Point::pointTransformedBy(points[0], origin, rotationAngle);
        y_min = localPoint->getLocalY();
        delete localPoint;

        for(int i=0; i<points.length(); i++) {
            Point *localPoint = Point::pointTransformedBy(points[i], origin, rotationAngle);
            if(localPoint->getLocalY() < y_min)
                y_min = localPoint->getLocalY();
            delete localPoint;
        }

        return y_min;
 }

 /**
  * @brief Surface::yMaxInCoordinateSystem
  * @param origin - origin of coordinate system (usually section centroid)
  * @param rotationAngle = angle by which coordinate system has been rotated
  * @return
  * Method finds y_max for current surface. y_MAX = max{ y_max1, y_max2, .., y_maxn}
  * will be the coordinate in local coordinate system (central coordinate system rotated by angle 0 (OMEGA))
  * where we will have top boundary of STRAIN PROFILE ex. eps_cu in this coordinate i.e. ultimate compression strength
  */
 double Surface::yMaxInCoordinateSystem(Point *origin, double rotationAngle)
 {
     double y_max; //will store y_max

     if(!isAdaptiveLinearized)
         performAdaptiveLinearization();

     if(points.length() <1)
         throw NoVerticesInSurfaceException("Brak wierzchołków zdefiniowanych dla powierzchni");

     Point *localPoint = Point::pointTransformedBy(points[0], origin, rotationAngle);
     y_max = localPoint->getLocalY();
     delete localPoint;

     for(int i=0; i<points.length(); i++) {
         Point *localPoint = Point::pointTransformedBy(points[i], origin, rotationAngle);
         if(localPoint->getLocalY() > y_max)
             y_max = localPoint->getLocalY();
         delete localPoint;
     }

     return y_max;
 }

 double Surface::internalAction_N(StrainProfile *profile)
 {
     if(isOpening)
         return 0;
     return performStressIntegration(profile, ACTION_N);
 }

 double Surface::internalAction_Mx(StrainProfile *profile)
 {
     if(isOpening)
         return 0;
     return performStressIntegration(profile, ACTION_Mx);
 }

 double Surface::internalAction_My(StrainProfile *profile)
 {
     if(isOpening)
         return 0;
     ///SHOULDN't BE -My here !?
     return -performStressIntegration(profile, ACTION_My);
 }

 /**
  * @brief Surface::performStressIntegration
  * @param profile - StrainProfile for which to calculate given internal action (N, Mx, My)
  *                  for this surface
  * @param actionType - type of calculated internal action: N, Mx, My
  * @return - internal action (N, Mx, My) value as double
  * To perform stress integration we apply numerical integration scheme:
  * the Gauss-Legendre quadrature method. Area integral of general form:
  * Rs,i = SSs,i(x^r*y^s*sigma(y)dxdy)
  */
 double Surface::performStressIntegration(StrainProfile *profile,
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

     //By applying Green's theorem, the area integral is transformed
     //into a line integral along the closed boundary Li that enclosed the area Si
     // Rs,i = SSs,i(x^r*y^s*sigma(y)dxdy) -> 1/(r+1) * SL,i x^(r+1)*y^s*sigma(y)dy

     //Line integral along the boundary Li which is a closed polygon of nli segments
     //can be expressed as sum of integrals of each individual polygon segment lj
     //Rs,i = 1/(r+1) * SL,i x^(r+1)*y^s*sigma(y)dy -> 1/(r+1) * SUM(j=1, nli, Sl,j x^(r+1)*y^s*sigma(y)*dy)
     //where Sl,j x^(r+1)*y^s*sigma(y)*dy is Ij  (elementary integral for each individual segment)

     //loop through each closed polygon segment
     //and adding up action contribution fromm each integrated polygon segment (Ij)
     //SUM(j=1, nli, Sl,j x^(r+1)*y^s*sigma(y)*dy)
     double xi,yi, xi_1, yi_1;
     int numOfPoints = this->numberOfPoints();
     Point *sectionCentroid = profile->getSection()->getGeometricalCenter();
     double omega = profile->getRotationAngle();

     for(int i=0; i< numOfPoints; i++) {

         //first polygon segment end coordinates
         xi = points[i]->getX();
         yi = points[i]->getY();

        //second polygon segment end coordinates
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
        Segment lj(pi,pi_1, std::to_string(i));

        //adding up each elementary integral Ij for each polygon segment
        //we pass current segment object lj, strain profile and
        //exponents r and s which depends on calculated internal action N, Mx, My
        internalAction += segmentIntegral_Ij(lj, profile, r, s);
     }

     //multiplying resultant internal action by suitable coefficient
     internalAction *= (double) 1/(r+1);

     return internalAction;
 }

 /**
  * @brief Surface::segmentIntegral_Ij
  * @param lj - j polygon segment on which stress integration is performed
  * @param profile - strain profile for the current surface from which we get
  *                  stress values to integrate sigma(y)
  * @param r - exponent for x coordinate
  * @param s - exponent for y coordinate
  *            r,s exponents use to differentiate between N, Mx, My calculation
  * @return resultng internal action from polygon segment integration
  */
 double Surface::segmentIntegral_Ij(const Segment &lj, StrainProfile *profile, double r, double s)
 {
        double internalAction = 0.0; //elementary Integral Ij

        //line segment may by expressed as lj -> x = aj*y + bj
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

        //initialization for horizontal line
        double aj = 0.0;
        double bj = (xi + xi_1)/2; //midpoint coordinate of horizontal segment?
        if((yi_1 - yi) != 0) {
            //if linear funtion isn't horizontal!
            aj = (xi_1 - xi)/(yi_1 - yi);
            bj = xi - aj*yi;
        }


        //introducing  x = aj*y + bj in elementary integral Ij makes it further simplified
        //into a line integral  of a single-variable (y) nonlinear function Fj(y):
        //Ij = Sl,j (aj*y + bj)^(r+1) *y^s *sigma(y)dy = Sl,j Fj(y)dy

        //now Gausse-Legendre integration scheme for Ij (elementary integral) is applied
        //nG = n\2 + 1 - number of Gaussian points utilized for polynomial of order n
        //it is highly probable that multiple function parts (f1,f2..) in sigma-epsilon
        //will be mapped on the same polygon segment leading to non-continuity in derivative
        //SAMPLING PROCCESS should be applied per function part f,k within each segment
        //in order to guarantee exact solutions for the widely used piecewise-polynomial
        //stress-strain relationships

        double eps_o = profile->getStrainAtOrigin_eps_o();
        double fi = profile->getCurvature_fi();


        //setting segment limits y_Aj and y_Bj in such way that
        //y_Aj is minimal y coordinate in segment and
        //y_Bj is maximal y_coordinate in segment
        double y_Aj = std::min(yi, yi_1); //where yi and yi_1 are starting and ending points y coordinates of the segment
        double y_Bj = std::max(yi, yi_1); //like in above comment

        for(int k=0; k < material->numberOfFunctionParts(); k++)
        {
            //reading strain limits for each function part in material definition
            double eps_limA_k = material->functionPartStrainLimitA(k); //is more tensile limit
            double eps_limB_k = material->functionPartStrainLimitB(k); //is more compressive limit
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
            //!!!! PROBLEM WHETHER y_limA_k should be < y_limB_k always to get      !!!!
            //!!!! CORRECT INTERNAL ACTION CONTRIBIUTION OR IT CAN BE IN ANY ORDER? !!!!
            //testing obtained y_limA_k, y_limB_k coordinate limits whether they fall outside
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
            //the same segment limit they have gotten the same y_Aj or yBj segment limit
            if(y_limA_k == y_limB_k)
                continue; //we skip integration of this function part f,k as its corresponding strain range fall outside the segment

             //now we have to integrate fragment of segment between y_limA_k and y_limB_k
             //using function part f,k  to calculate stresses in corresponding points
             //based on y contained in [y_limA_k, y_limB_k] y -> epsilon -> f,k(epsilon) -> stress to integrate
             //now Gaussian sampling is applied for all contribiuting function parts f,k of the considered segment lj
             //Ij_k - elementary integral kth component
            double Ij_k = 0.5*(y_limB_k - y_limA_k);  //!!! shouldn't be y_limB_k and y_limA_k adjusted in such way that y_limB_k > y_limA_k ? !!!
            /* test only
            if(s==1) { //calculating Mx we take into account absoulte value od (y_limB_k - y_limA_k) ?
                Ij_k = 0.5*std::abs(y_limB_k - y_limA_k);
            }*/

            //calculation of elementary integral Ij_k kth component depends on number of utilized Gauss points nG (quadrature order)
            int fkOrder = material->functionPartDegree(k);
            int nG = (fkOrder + r + 1 + s)/2 + 1; //number of required Gaussian points to get exact results for polynomial integration of order fkOrder!

            double summation = 0.0;
            //we generate gauss object from which we will be obtaining weights w_m and y_m coordinates
            GaussLegendre gauss(nG);
            for(int m = 0; m < nG; ++m)
            {
                double ym = 0.5*(y_limB_k + y_limA_k) + gauss.lambda(m)*0.5*(y_limB_k - y_limA_k);
                double eps_ym = eps_o - fi*ym;

                if(DEBUG)
                     fprintf(stderr, "(eps_ym, ym) = (%g, %g)\n", eps_ym, ym);

                double Fj_ym  = std::pow(aj*ym + bj, r+1)*std::pow(ym,s)*material->sigmaEpsilon(eps_ym);
                //adding up consecutive parts of Gaussian quadrature
                summation += gauss.weight(m)*Fj_ym; //N -> [N/mm] Mx, My -> [N]
            }
            Ij_k *= summation; //N -> [N] Mx, My -> [Nmm]

            //cumulating elementary integral components in internalAction variable
            internalAction += Ij_k; //N-> [N] Mx, My -> [Nmm]

            //loging function part integration
            if(DEBUG)
                fprintf(stderr, "Material Definition Function Part f%d Integration Result: %g.\n", k, Ij_k);
        }
        /*
        switch(material->getMaterialType())
        {
            case CONCRETE:
                break;
            case STRUCTURAL_STEEL:
                break;
            case REINFORCEMENT:
                break;
            case REINFORCEMENT_BARS:
                break;
        }*/

        if(DEBUG)
            fprintf(stderr, "Wartość całkowania dla %s-tego segmentu I%s = %g.\n", lj.getName().c_str(), lj.getName().c_str(), internalAction);

        return internalAction; //elementary integral Ij for segment lj

 }
