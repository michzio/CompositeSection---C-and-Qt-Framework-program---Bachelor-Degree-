#ifndef SECTIONSOLVER_H
#define SECTIONSOLVER_H

#include "section.h"
#include <QList>
#include "vertex3d.h"

class SectionSolver : public Section
{
    const double STRAIN_ERROR_TOLERANCE = 0.0;
    const double STRAIN_ERROR_TOLERANCE2 = -0.0000000001;
    const double STRAIN_ERROR_TOLERANCE3 = 0.0;
    //const double Y_CU_STEP = 0.0001;
    //const double Y_TU_STEP = 0.0001;
    //const double FI_STEP = M_PI/180/10000/10;
    const double Y_CU_STEP = 0.0001;
    const double Y_TU_STEP = 0.0001;
    const double FI_STEP = M_PI/180/10000/10;
    double tensionLimit; //tensile strain limit calculated for given strain profile rotation angle
    double yCoordinateOfTensionLimit; //y coordinate where tension limit occures in the section
    double compressionLimit; //compressive strain limit calculated for given strain profile rotation angle
    double yCoordinateOfCompressionLimit; //y coordinate where compression limit occures in the section
    double pureCompressionLimit; //pure compression strain limit calculated for given compressive strain limit
    double yCoordinateOfPureCompressionLimit; //y coordinate of pure compression limit in the section
    double rotationAngle; //strain profile coordinate system rotation in relation to section definition
    QList<Vertex3D *> vertices;

public:
    SectionSolver();

    //necessary methods to apply force equilibrum conditions and calculate the incremental
    //and ultimate section response (solution strategies).

    //method perform strain profile steerage between strain compressive and tensile limits
    void solveSection(void);

    //method setting new strain profile rotation angle
    void setRotationAngle(double omega);
    double getRotationAngle(void);

    //method returning vertices (Mx, My, N) constructed from calculated internal actions
    QList<Vertex3D *> &getVertices(void);
    //function sorting QList of Vertex3D vertices
    void sortVerticesByN(void);
    void sortVerticesByMx(void);
    void sortVerticesByMy(void);
    void sortVerticesByAlfa(void);

    //method searching current section tensile and compressive strain limit
    //based on components strain limits
    void findTensileStrainLimit(double *eps_tu, double *y_tu);
    void findCompressiveStrainLimit(double *eps_cu, double *y_cu);
    //method finding pure compression strain limit
    void findPureCompressionStrainLimit(double *eps_co, double *y_co);
    //calculation tensile and compressive strain limits points
    void calculateTensileAndCompressiveStrainLimits(void);
    //methods construct StrainProfile object from 2 points (eps_y_tu, y_tu) and (eps_y_cu, y_cu) and perform its integration
    //calculating resultant internal actions N, Mx, My
    void integrateStrainProfileFrom(double eps_y_tu, double y_tu, double eps_y_cu, double y_cu);
};

#endif // SECTIONSOLVER_H
