#ifndef SURFACE_H
#define SURFACE_H

#include "points.h"
#include <QVarLengthArray>
#include "concrete.h"
#include "structuralsteel.h"
#include "fafitisconcrete.h"
#include "hollow.h"
#include<cmath>
#include "vector.h"
#include "arc.h"
#include<string>
#include <QGraphicsView>
#include <QGraphicsPolygonItem>
#include <QGraphicsScene>
#include "segment.h"

const double A_TOL = 0.01; // 1% tolerance of area approximation

class Surface : public Points
{

    //class inherits from Points and adds additional properties to coordinates e.g. angles
    QVarLengthArray<double> angles; //set of angel values of the succeeding circular arc segment
                                    //ignored (set to 0) for straight segments
                                    //negative angel value leads to subtraction of the corresponding circular
                                    //sector from the surface
    Material *material; //instance of the material object
    std::string name;
    bool isOpening = false;
    bool isAdaptiveLinearized = false;
    Point *centroid = NULL;
    //int verticesNumber = 0;
public:
    Surface(std::string);
    Surface(std::string, Material *);
    Surface(Surface& other); //copy-constructor

    void setName(std::string);
    void setMaterial(Material *);
    void setIsOpening(bool);
    void addPointAndAngel(double, double, double);
    void clearPointsAndAngles(void);
    double *getAngels(void);
    Material *getMaterial(void);
    bool getIsOpening(void);
    std::string getName(void);
    //int getVerticesNumber(void);
    bool hasBeenAdaptiveLinearized(void);
    void setIsAdaptiveLinearized(bool flag);

    //method to perform Adaptive Linearization
    void performAdaptiveLinearization(void);
    double calculatePolygonArea();
    double calculatePolygonAreaIncludingArcs();
    double arcSegmentArea(const Point* start, const Point* end, double angle);
    double lengthOfVector(const Point* start, const Point* end);
    double radiusOfArc(double lengthOfChord, double angle);
    Point *middlePointOfSegment(const Point *p0, const Point *p1);
    Point **centersOfCircleGiven2PointsAndRadius(const Point *start, const Point *end, double radius);
    Point *selectProperCenterOfCircle(Point **centers, Point *start, Point *end, double angle);
    bool pointIsInsieThePolygon(const Point *p);
    double calculateAngle2D(Vector *v1, Vector *v2);
    Point *arcMidPoint(Point *startPoint, Point *endPoint, double angle);
    double secans(double x);
    Arc *getArcWithLargestSegmentArea(int *pos);

    //method to calculate geometrical center of surface
    double getGeometricalCenterX(void);
    double getGeometricalCenterY(void);
    Point* getGeometricalCenter(void);

    //calculating first moment of area in relation to coordinate system in point P
    double firstMomentOfArea_Sx();
    double firstMomentOfArea_Sx(Point *p); //y*A
    double firstMomentOfArea_Sy();
    double firstMomentOfArea_Sy(Point *p); //x*A

    //function finds intersection of two polygons
    //this polygon (after adaptive linearization
    //and polygon passed to function as first parameter
    QVarLengthArray<Surface *> *intersection(Surface *);

    //functions searching y_min and y_max values amongst surface
    //vertices for coordinate system defined by parameters
    double yMaxInCoordinateSystem(Point *origin, double rotationAngle);
    double yMinInCoordinateSystem(Point *origin, double rotationAngle);

    //method calling performStressIntegration and getting each internal action N, Mx, My for this surface
    //and for given StrainProfile
    double internalAction_N(StrainProfile *);
    double internalAction_Mx(StrainProfile *);
    double internalAction_My(StrainProfile *);

    //method to perform Stress Integration of this surface
    double performStressIntegration(StrainProfile *, InternalAction);
    //method to perform Stress Integration of each polygon segment lj
    double segmentIntegral_Ij(const Segment& lj, StrainProfile* profile,
                              double r, double s);

    ~Surface();

    //drawing surface on QGraphicsScene
private:
    QGraphicsPolygonItem *graphicsItem = NULL;
    QGraphicsScene *scene = NULL;
public:
    void setGraphicsItem(QGraphicsPolygonItem *gi)
    {
        graphicsItem = gi;
    }

    QGraphicsPolygonItem *getGraphicsItem(void) {
        return graphicsItem;
    }

    void setGraphicsScene(QGraphicsScene *s) {
        scene = s;
    }
    //Surface drawing itself on QGraphicsScene passed to this object
    void drawSurface(Qt::GlobalColor fillColor = Qt::red) {
        //declaring QPolygonF object which will be
        //added to QGraphicsScene
        QPolygonF polygon;
        Point **point = getPointsArray();

        if(getGraphicsItem() != NULL)
            scene->removeItem(getGraphicsItem());

        fprintf(stderr, "Liczba punktow: %d\n", numberOfPoints());
        //filling QPolygonF with QPointF
        for(int i=0; i<numberOfPoints(); i++)
        {
            fprintf(stderr, "{%g,%g}\n",  point[i]->getX(),
                                          point[i]->getY());
            polygon.append(QPointF( static_cast<int>(point[i]->getX() + 0.5),
                                    static_cast<int>(-point[i]->getY() + 0.5)));
        }
        QGraphicsPolygonItem* polygonItem =
                scene->addPolygon(polygon, QPen(Qt::black), QBrush(fillColor));
        setGraphicsItem(polygonItem);
    }

    void drawGeometricalCenter(void) {
        double rad = 1;

        Point *point = getGeometricalCenter();
        scene->addEllipse(point->getX()-rad, -point->getY() -rad, rad*5.0, rad*5.0,
                    QPen(), QBrush(Qt::SolidPattern));
    }
};

//Exception thrown when there are performed calculations on
//surface which doesn't have any vertices
struct NoVerticesInSurfaceException : public std::exception
{
   std::string s;
   NoVerticesInSurfaceException(std::string ss) : s(ss) {}
   ~NoVerticesInSurfaceException() throw () {} // Updated
   const char* what() const throw() { return s.c_str(); }
};

#endif // SURFACE_H
