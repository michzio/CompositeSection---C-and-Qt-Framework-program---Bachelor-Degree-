#ifndef ARC_H
#define ARC_H

#include "point.h"
#include<cstdlib>

class Arc
{
    Point *startPoint;
    Point *endPoint;
    double angle;

public:
    Arc(Point *start = NULL, Point *end = NULL, double angle = 0.0) : startPoint(start), endPoint(end), angle(angle) { }
    void setStartPoint(Point *);
    void setEndPoint(Point *);
    void setAngle(double angle);
    Point *getStartPoint();
    Point *getEndPoint();
    double getAngle();
};

#endif // ARC_H
