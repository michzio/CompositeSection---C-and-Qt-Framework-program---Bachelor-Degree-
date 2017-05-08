#ifndef SEGMENT_H
#define SEGMENT_H

#include "point.h"
#include<string>

class Segment
{
    Point start;
    Point end;
    std::string name;
public:
    //Points should have set the same local coordinate system
    Segment(Point s, Point e) : start(s), end(e) { }
    Segment(Point s, Point e, std::string n) : start(s), end(e), name(n) { }

    void setStartPoint(const Point& p);
    void setEndPoint(const Point& p);
    void setName(std::string n);

    const Point& getStartPoint(void) const;
    const Point& getEndPoint(void) const;
    std::string getName(void) const;

    double length(void) const;
};

#endif // SEGMENT_H
