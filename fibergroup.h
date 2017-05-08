#ifndef FIBERGROUP_H
#define FIBERGROUP_H

#include "points.h"
#include "reinforcement.h"
#include "fafitisreinforcement.h"
#include <QVarLengthArray>
#include <string>
#include <QGraphicsScene>
#include <QGraphicsView>

class FiberGroup : public Points
{
     //class inherits from Points and adds additional properties to coordinates e.g. contributing areas
    Material *material;
    QVarLengthArray<double> areas; //set of contribution areas of the succeeding fibers in the group
    std::string name;
    //int numberOfFibers;
public:
    FiberGroup(std::string);
    FiberGroup(std::string, Material *);

    void addPointAndArea(double, double, double);
    void setMaterial(Material *);
    Material* getMaterial(void);
    double *getAreas(void);
    void setName(std::string);
    std::string getName(void);
    //int getNumberOfFibers(void);
    void clearPointsAndAreas(void);
    Point *getOrigin(void);

    //methods to find y coordinate of lower and upper edge of component
    double yMaxInCoordinateSystem(Point *origin, double rotationAngle);
    double yMinInCoordinateSystem(Point *origin, double rotationAngle);

    //methods calling performStressIntegration and getting each internal
    //action N, Mx, My for this fiber group and for given StrainProfile
    double internalAction_N(StrainProfile *profile);
    double internalAction_Mx(StrainProfile *profile);
    double internalAction_My(StrainProfile *profile);
    //method to perfrom Stress Integration
    double performStressIntegration(StrainProfile *profile, InternalAction actionType);

    //method to get cumulative area of Fibergroup element
    double getTotalArea(void);

    //method to calculate geometrical center of surface
    double getGeometricalCenterX(void);
    double getGeometricalCenterY(void);
    Point *getGeometricalCenter(void);

    //calculating first moment of area in relation to coordinate system in point P
    double firstMomentOfArea_Sx(void); //default coordinate system
    double firstMomentOfArea_Sx(Point *p); //coordinate system at point p
    double firstMomentOfArea_Sy(void); //default coordinate system
    double firstMomentOfArea_Sy(Point *p); //coordinate system at point p

    double getXCenterOfFiber(int i);
    double getYCenterOfFiber(int i);
    Point *getCenterOfFiber(int i);

    ~FiberGroup();

private:
    QGraphicsScene *scene = NULL;
public:
    void setGraphicsScene(QGraphicsScene *s) {
        scene = s;
    }

    void drawFiberGroup(Qt::GlobalColor fillColor = Qt::black)
    {
        Point **point = getPointsArray();
        double *area = this->getAreas();

        for(int i=0; i<numberOfPoints(); i++)
        {
            double rad = sqrt(area[i]/M_PI);

            scene->addEllipse(point[i]->getX()-rad, -point[i]->getY() -rad, rad*2.0, rad*2.0,
                        QPen(), QBrush(fillColor, Qt::SolidPattern));
        }
    }
};

//Exception thrown when there are performed calculations on
//fiber group which doesn't have any vertices
struct NoVerticesInFiberGroupException : public std::exception
{
   std::string s;
   NoVerticesInFiberGroupException(std::string ss) : s(ss) {}
   ~NoVerticesInFiberGroupException() throw () {} // Updated
   const char* what() const throw() { return s.c_str(); }
};

#endif // FIBERGROUP_H
