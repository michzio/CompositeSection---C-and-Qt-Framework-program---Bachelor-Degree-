#include "point.h"
#include<cmath>

//constructor
Point::Point(double x, double y)
{
    this->x = x;
    this->y = y;
}

//this method sets x coordinate of the point
void Point::setX(double x)
{
    this->x = x;
}

//this method sets y coordinate of the point
void Point::setY(double y)
{
    this->y = y;
}

//this method returns x coordinate of the point
double Point::getX() const
{
  return x;
}
//this method returns y coordinate of the point
double Point::getY() const
{
  return y;
}

/**
 * @brief Point::setOrigin
 * @param origin
 * Method sets coordinates of origin point of coordinate system
 */
void Point::setOrigin(Point *origin)
{
    originX = origin->getX();
    originY = origin->getY();
}

/**
 * @brief Point::setOrigin
 * @param originX - displacement of coordinate system in x direction
 * @param originY - displacement of coordinate system in y direction
 * Method sets coordinates of origin point of coordinate system
 */
void Point::setOrigin(double originX, double originY)
{
    this->originX = originX;
    this->originY = originY;
}

/**
 * @brief Point::setRotation
 * @param angle - rotation angle in radians
 * Method sets rotationAngle of new coordinate system to which
 * returned localX(), localY() will be transfered from default
 * cooraginatr system (angle = 0.0)
 */
void Point::setRotation(double angle)
{
    rotationAngle = angle;
}

/**
 * @brief Point::setOriginAndRotation
 * @param originX
 * @param originY
 * @param angle
 * Method is jointing two above methods
 */
void Point::setOriginAndRotation(double originX, double originY, double angle)
{
    this->originX = originX;
    this->originY = originY;
    rotationAngle = angle;
}

/**
 * @brief Point::getOrigin
 * @return Point containing coordinates of origin of current coordinate system
 */
Point *Point::getOrigin(void)
{
    return new Point(originX, originY);
}

/**
 * @brief Point::getOriginX
 * @return
 */
double Point::getOriginX(void)
{
    return originX;
}

/**
 * @brief Point::getOriginY
 * @return
 */
double Point::getOriginY(void)
{
    return originY;
}

/**
 * @brief Point::getRotation
 * @return
 */
double Point::getRotation(void)
{
    return rotationAngle;
}


/**
 * @brief Point::getLocalX
 * @return
 * Method is doing calculations on coordinate x, y in order to
 * provide x' coordinate in local translated and roteted coordinate system
 * typical pattern in which this function will be used for X0Y -> XCY -> xCy
 * 0 (OMEGA) - angle by which coordinate system is roteted is positive while
 * rotation is counterclockwise (CCW) i.e. 0 = 90degree if X axis is rotating onto Y axis around point C
 *
 */
double Point::getLocalX(void) const
{
     //x is in global, default coordinate system
     //not rotated not moved!
     //it is coordinate in user input coordinate system

     //this point can have set coordinate system new origin
     //originX, originY e.g. geometrical center of section
     //i.e. central coordinate system

     //1) TRANSLATION OF COORDINATE
     // x' = x - originX
     // y' = y - originX
     double localX = x - originX;
     double localY = y - originY;

     //2) ROTATION BY ANGLE OMEGA
     // transformation matrix
     // | x'| = | cos0    sin0 | | x |
     // | y'| = | -sin0   cos0 | | y |
     double rotetedX = localX*std::cos(rotationAngle) + localY*std::sin(rotationAngle);
     //double rotetedY = - localX*std::sin(rotationAngle) + localY*std::cos(rotationAngle);

     return rotetedX;
}

double Point::getLocalY(void) const
{
    //description as in getLocalX()
    double localX = x - originX;
    double localY = y - originY;

    double rotetedY = - localX*std::sin(rotationAngle) + localY*std::cos(rotationAngle);

    return rotetedY;
}


/**
 * @brief getXTransformedTo
 * @param origin
 * @param angle
 * @return
 * Function transforms this Point to coordinate system defined by
 * params origin and rotation angle
 */
double Point::getXTransformedTo(Point origin, double angle)
{
      double localX = x - origin.getX();
      double localY = y - origin.getY();

      double rotetedX = localX*std::cos(angle) + localY*std::sin(angle);

      return rotetedX;
}

double Point::getYTransformedTo(Point origin, double angle)
{
    double localX = x - origin.getX();
    double localY = y - origin.getY();

    double rotetedY = -localX*std::cos(angle) + localY*std::sin(angle);
    return rotetedY;
}


//STATIC METHODS
Point *Point::pointTransformedBy(Point *point, Point *origin, double rotation)
{
    Point *transformedPoint = new Point(point->getX(), point->getY());
    transformedPoint->setOriginAndRotation(origin->getLocalX(), origin->getY(), rotation);

    return transformedPoint;
}
