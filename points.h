#ifndef POINTS_H
#define POINTS_H

#include "point.h"
#include <QVarLengthArray>
#include<exception>
#include<string>

class StrainProfile;

class Points
{
protected:
    enum InternalAction { ACTION_N , ACTION_Mx, ACTION_My };
    //set of coordinates (xj, yj)
    QVarLengthArray<Point *> points;

public:
    Points();  
    void addPoint(double x, double y);
    Point **getPointsArray(void);
    int numberOfPoints();
    //geometry manipulation methods
};

struct AreaIndexOutOfBounds : public std::exception
{
    std::string s;

    AreaIndexOutOfBounds(std::string ss) : s(ss) {}
    ~AreaIndexOutOfBounds() throw () {}

    const char* what() const throw() { return s.c_str(); }
};

#endif // POINTS_H
