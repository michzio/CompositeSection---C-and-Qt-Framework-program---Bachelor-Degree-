#ifndef LINE_H
#define LINE_H

#include "points.h"
#include "reinforcement.h"
#include "frp.h"
#include "qvarlengtharray.h"
#include "segment.h"
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QPainter>

/************************************************************
 * Lines are defined as a series of vertex coordinates      *
 * (multi-segment), including the contribution area of the  *
 * succeeding segment Aj.                                   *
 ************************************************************/
class Line : public Points
{
     //class inherits from Points and adds additional properties to coordinates e.g. contributing areas
    Material *material;
    QVarLengthArray<double> areas; //set of contributing areas of the succeeding line segments Aj
    std::string name;
    //int numberOfPoints;

public:
    Line(std::string lineName);
    Line(std::string lineName, Material *);
    void addPointAndArea(double, double, double);
    double *getAreas(void);
    //int getNumberOfPoints(void);
    void setMaterial(Material *m);
    void clearPointsAndAreas(void);
    Material *getMaterial(void);
    std::string getName(void);
    void setName(std::string);
    Point *getOrigin(void);
    ~Line();

    //methods to find y coordinate of lower and upper edge of component
    double yMaxInCoordinateSystem(Point *origin, double rotationAngle);
    double yMinInCoordinateSystem(Point *origin, double rotationAngle);

    //methods calling performStressIntegration and getting each
    //internal action N, Mx, My for this line element and for given
    //StrainProfile
    double internalAction_N(StrainProfile *profile);
    double internalAction_Mx(StrainProfile *profile);
    double internalAction_My(StrainProfile *profile);

    //method to perform Stress Integration of this line element
    double performStressIntegration(StrainProfile *, InternalAction);
    //method to perform Stress Integration of each line element segment lj
    double segmentIntegral_Ij(const Segment& lj, StrainProfile* profile,
                              double r, double s, int lineIdx);


    //method to get cumulative area of Line element
    double getTotalArea(void);

    //method to calculate geometrical center of surface
    double getGeometricalCenterX(void);
    double getGeometricalCenterY(void);
    Point *getGeometricalCenter(void);

    //calculating first moment of area in relation to coordinate system in point P
    double firstMomentOfArea_Sx(void);
    double firstMomentOfArea_Sx(Point *p); //y*A
    double firstMomentOfArea_Sy(void);
    double firstMomentOfArea_Sy(Point *p); //x*A

    double getComponentCenterX(int i);
    double getComponentCenterY(int i);
    Point *getComponentCenter(int i);

private:
    QGraphicsScene *scene = NULL;
public:
    void setGraphicsScene(QGraphicsScene *s) {
        scene = s;
    }

    void drawLine( Qt::GlobalColor fillColor = Qt::black)
    {
        Point **point = this->getPointsArray();
        double *area = this->getAreas();

        for(int i=0; i<numberOfPoints() - 1; i++)
        {
            double length = sqrt( pow(point[i+1]->getX() - point[i]->getX(), 2.0)
                                + pow(point[i+1]->getY() - point[i]->getY(), 2.0));
            //thickness of line
            double t = area[i]/length;
            QPen pen1;
            pen1.setColor(QColor(fillColor));
            pen1.setWidth((int)t);
            pen1.setCosmetic(false);
            scene->addLine(point[i]->getX(), -point[i]->getY(),
                           point[i+1]->getX(), -point[i+1]->getY(),
                           pen1);
        }
    }
};

//Exception thrown when there are performed calculations on
//line element which doesn't have any vertices
struct NoVerticesInLineException : public std::exception
{
   std::string s;
   NoVerticesInLineException(std::string ss) : s(ss) {}
   ~NoVerticesInLineException() throw () {} // Updated
   const char* what() const throw() { return s.c_str(); }
};

#endif // LINE_H
