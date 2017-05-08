#include "points.h"

Points::Points()
{
}

void Points::addPoint(double x, double y)
{
    points.append(new Point(x,y));
}

Point **Points::getPointsArray()
{
    Point **data = points.data();
    return data;
}

int Points::numberOfPoints()
{
  return points.count();
}

