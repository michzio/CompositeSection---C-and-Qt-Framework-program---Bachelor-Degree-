#include "vector.h"

void Vector::setStartPoint(Point *start)
{
    this->startPoint = start;
}

void Vector::setEndPoint(Point *end) {
    this->endPoint = end;
}

Point *Vector::getStartPoint(void)
{
     return startPoint;
}

Point *Vector::getEndPoint(void) {
    return endPoint;
}

double Vector::getX(void) {

    return (endPoint->getX() - startPoint->getX());
}
double Vector::getY(void) {
    return (endPoint->getY() - startPoint->getY());
}

/****************************************************
 * Find vector perpendicular to given using formula *
 * [x, y] ------> [y, -x].                          *
 * [x1-x0, y1-y0] ------> [y1-y0, x0-x1]            *
 * start(x0, y0) ------> start(y0,x1)               *
 * end(x1, y1) ------> end(y1, x0)                  *
 ****************************************************/
Vector *Vector::perpendicularVector()
{
    //start(y0,x1)
    Point *newStart = new Point(startPoint->getY(), endPoint->getX());
    //end(y1, x0)
    Point *newEnd = new Point(endPoint->getY(), startPoint->getX());

    return (new Vector(newStart, newEnd));
}

/*****************************************************
 * Function normalize vector. It is done by dividing *
 * x and y vector coordinates by norm ||v|| where    *
 * ||v|| = sqrt(x^2 + y^2) next we add this new      *
 * vector coordinates x' = x/||v||, y' = y/||v|| to  *
 * cooridinates of startPoint and assigned them to   *
 * endPoint object coordinates.                      *
 *****************************************************/
void Vector::normalizeVector()
{
    double x = this->getX();
    double y = this->getY();
    double normOfVector = this->normOfVector();

    double newX = x/normOfVector;
    double newY = y/normOfVector;

    double endX = this->startPoint->getX() + newX;
    double endY = this->startPoint->getY() + newY;

    this->endPoint->setX(endX);
    this->endPoint->setY(endY);
}

double Vector::normOfVector()
{
    double x = this->getX();
    double y = this->getY();
    return std::sqrt(x*x + y*y);
}
Vector::~Vector() {
    delete startPoint;
    delete endPoint;
}
