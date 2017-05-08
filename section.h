#ifndef SECTION_H
#define SECTION_H

#include "surface.h"
#include "line.h"
#include "fibergroup.h"
#include<QVarLengthArray>
#include <vector>

/**********************************************
 * Section class contains an unlimited number *
 * of section component objects (Surfeces,    *
 * Line, FiberGroup) together with methods to *
 * define centroids, determine the ULS, and   *
 * summate the internal actions from the      *
 * individually integreted component          *
 * contributions.                             *
 **********************************************/
class Section
{
 protected:
    //assumption: Lines and Fibergroups are in the foregroud of surfaces
    //            Surfaces in array are in order from background material to foreground
    //            User while defining this surfaces must put them in correct order!
    //            During calculation program will be subtracting multi-nested materials
    //            and openings.
    //            Material of first surface in "surfaces" array is referencing material
    //            While calculating weighted area (A*) and weighted first moment of area (Sx*, Sy*) its
    //            modul of elasticity (E) will be used in accordance with formula: n_i = E_i/E
    //            where E_i - is modul of elasticity ith material and n_i is weight for ith material
    QVarLengthArray<Surface*> surfaces;
    QVarLengthArray<Line*> lines;
    QVarLengthArray<FiberGroup*> fiberGroups;
private:
    //section centoid (geometrical center)
    Point *centroid = NULL;

public:
    Section();

    void addSurface(Surface *);
    void addLine(Line *);
    void addFiberGroup(FiberGroup *);

    Surface **getSurfaces();
    Line **getLines();
    FiberGroup **getFiberGroups();

    int getNumberOfSurfaces();
    int getNumberOfLines();
    int getNumberOfFiberGroups();

    //method to define centroid
    Point *defineCentroid() {
        return getGeometricalCenter();
    }

    //methods to summate the internal actions from the individually integrated component contributions
    double summateInternalActions_N(StrainProfile *profile);
    double summateInternalActions_Mx(StrainProfile *profile);
    double summateInternalActions_My(StrainProfile *profile);

    //methods to calculate wieghted area (A*) and weighted first moment of area (Sx*, Sy*)
    //which are used to get proper centroid of section consisted with different materials
    double weightedArea(void);
    double weightedFirstMomentOfArea_Sx(void); //in default coordinate system
    double weightedFirstMomentOfArea_Sx(Point *p); //in coordinate system with origin p
    double weightedFirstMomentOfArea_Sy(void); //in default coordinate system
    double weightedFirstMomentOfArea_Sy(Point *p); //in coordinate system with origin p

    //methods gets centroid from Section state or calculate it if
    //hasn't been calculated yet
    double getGeometricalCenterX(void);
    double getGeometricalCenterY(void);
    Point *getGeometricalCenter(void);

    //methods that calculate centroid
    double calculateGeometricalCenterX(void);
    double calculateGeometricalCenterY(void);
    Point *calculateGeometricalCenter(void);

    //method draws centroid on pointed QGraphicsScene
    QGraphicsEllipseItem *drawCentoidOnScene(QGraphicsScene *scene);

    //methods searching y_min and y_max values amongst section
    //surfaces vertices in order to find lower and upper limit of
    //section in local coordinate system, i.e. central coordinate system
    //roteted by angle OMEGA.
    double yMaxInCoordinateSystem(Point *origin, double rotationAngle);
    double yMinInCoordinateSystem(Point *origin, double rotationAngle);

    //Following the integration process moments Mx and My are transferred back
    //from local (rotated) coordinate system xCy to not-rotated coordine system
    //XCY
    static double transferBackFromLocalMx(double localMx, double localMy, double omega);
    static double transferBackFromLocalMy(double localMx, double localMy, double omega);

private:
    //helper function in calculating weightedArea of section
    //used to include overlapping surfaces (intersections)
    double calculateIntersectionArea(Surface *,
                                     std::vector<int>::const_iterator,
                                     std::vector<int>::const_iterator,
                                     std::string&);
    //helper function in calculating weightedSx of section
    //used to include overlapping surfaces (intersections) contribution
    double calculateIntersectionSx(Surface *,
                                   std::vector<int>::const_iterator,
                                   std::vector<int>::const_iterator,
                                   Point *,
                                   std::string&);
    //helper function in calculating weightedSy of section
    //used to include overlapping surfaces (intersections) contribution
    double calculateIntersectionSy(Surface *,
                                   std::vector<int>::const_iterator,
                                   std::vector<int>::const_iterator,
                                   Point *,
                                   std::string&);
    //helper function in calculating total N of section
    //used to include overlapping surfaces (intersections) contribution
    double calculateIntersectionN(Surface *surface,
                                   std::vector<int>::const_iterator beginItr,
                                   std::vector<int>::const_iterator endItr,
                                   StrainProfile *profile,
                                   Material *material,
                                   std::string& logString);
    //helper function in calculating total Mx of section
    //used to include overlapping surfaces (intersections) contribution
    double calculateIntersectionMx(Surface *surface,
                                   std::vector<int>::const_iterator beginItr,
                                   std::vector<int>::const_iterator endItr,
                                   StrainProfile *profile,
                                   Material *material,
                                   std::string& logString);
    //helper function in calculating total My of section
    //used to include overlapping surfaces (intersections) contribution
    double calculateIntersectionMy(Surface *surface,
                                   std::vector<int>::const_iterator beginItr,
                                   std::vector<int>::const_iterator endItr,
                                   StrainProfile *profile,
                                   Material *material,
                                   std::string& logString);
};

struct NoSurfacesDefinedException : public std::exception
{
    std::string s;
    NoSurfacesDefinedException(std::string ss) : s(ss) {}
    ~NoSurfacesDefinedException() throw() {} //Updated
    const char* what() const throw() { return s.c_str(); }
};

#endif // SECTION_H
