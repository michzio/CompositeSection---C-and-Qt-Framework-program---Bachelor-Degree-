#ifndef POINT_H
#define POINT_H
/***************************************************************
 * This class is used to create Point objects, which are used  *
 * to define set of vertices (coordinates) of each Surface,    *
 * Line, FiberGroup object                                     *
 * *************************************************************/

class Point
{
    //default coordinate system constants
    constexpr static const double defaultOriginX = 0.0;
    constexpr static const double defaultOriginY = 0.0;
    //definition of coordinates in default (global coordinate system)
    double x;
    double y;
    //definition of coordinates of origin point of current coordinate system
    double originX = defaultOriginX;
    double originY = defaultOriginY;
    //definiton of angle of rotation of current coordinate system
    double rotationAngle = 0.0;
public:
    //Point object constructor
    Point(double x = 0, double y = 0);
    //setter and getter methods for object properties
    //manipulation of Point's coordinate in default coordinate system
    void setX(double x);
    void setY(double y);
    double getX() const;
    double getY() const;


    //setter and getter methods to set up displacement and rotation
    //of coordinate system in relation to default coordinate system
    void setOrigin(Point *origin);
    void setOrigin(double originX, double originY);
    void setRotation(double angle);
    void setOriginAndRotation(double originX, double originY, double angle);
    Point *getOrigin(void);
    double getOriginX(void);
    double getOriginY(void);
    double getRotation(void);

    //getter methods in relation to current coordinate system
    //(transform and rotated coordinate system)
    double getLocalX() const;
    double getLocalY() const;

    //return x,y  transformed to coordinate system defined by params
    double getXTransformedTo(Point origin, double angle);
    double getYTransformedTo(Point origin, double angle);

    //determines whether point lies to the left, right or on the segment
    //using cross product
    //| x2-x1  x3-x1 |
    //| y2-y1  y3-y1 |
    //if determinant > 0 (x3,y3) lies to the left of segment(line) [(x1,y1),(x2,y2)]
    //if determinant < 0 (x3,y3) lies to the right of segment(line) [(x1,y1),(x2,y2)]
    //if determinant = 0 (x3,y3) lies on the segment (line) [(x1,y1),(x2,y2)]
    bool liesToTheLeftOfSegment(Point start, Point end) {
         return ((end.x - start.x)*(this->y - start.y) - (end.y - start.y)*(this->x - start.x)) > 0;
    }
    bool liesToTheRightOfSegment(Point start, Point end) {
         return ((end.x - start.x)*(this->y - start.y) - (end.y - start.y)*(this->x - start.x)) < 0;
    }

    bool liesOnTheSegment(Point start, Point end) {
         return ((end.x - start.x)*(this->y - start.y) - (end.y - start.y)*(this->x - start.x)) == 0;
    }



    //static methods to perform point transformation
    static Point *pointTransformedBy(Point *point, Point *origin, double rotation);
};

#endif // POINT_H
