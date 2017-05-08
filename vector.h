#ifndef VECTOR_H
#define VECTOR_H

#include "point.h"
#include<cmath>

class Vector
{
    Point *startPoint;
    Point *endPoint;

public:
    Vector(Point *start = NULL, Point *end = NULL) : startPoint(start), endPoint(end) {}
    void setStartPoint(Point *start);
    void setEndPoint(Point *end);
    Point *getStartPoint(void);
    Point *getEndPoint(void);
    double getX(void);
    double getY(void);
    Vector *perpendicularVector();
    void normalizeVector();
    double normOfVector();

    ~Vector();
};

#endif // VECTOR_H
