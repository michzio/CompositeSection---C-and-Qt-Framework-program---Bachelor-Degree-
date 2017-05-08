#include "segment.h"
#include<cmath>

void Segment::setStartPoint(const Point& p)
{
    start = p;
}

void Segment::setEndPoint(const Point& p)
{
    end = p;
}

void Segment::setName(std::string n)
{
    name = n;
}

const Point& Segment::getStartPoint(void) const
{
    return start;
}
const Point& Segment::getEndPoint(void) const
{
    return end;
}

std::string Segment::getName(void) const
{
    return name;
}

double Segment::length(void) const
{
    double len = 0.0;
    //length of segment is calculated in point's local coordinate system
    double x = end.getLocalX() - start.getLocalX();
    double y = end.getLocalY() - start.getLocalY();

    len = std::sqrt( x*x + y*y);

    return len;
}
