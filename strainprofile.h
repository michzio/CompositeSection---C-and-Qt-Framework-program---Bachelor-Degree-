#ifndef STRAINPROFILE_H
#define STRAINPROFILE_H

#include "section.h"
#include <cmath>

class StrainProfile
{
    //angle by which local coordinate system is rotated
    double omega;
    //section for which strain profile is generated
    Section *section;
    //strain profile vertical limits
    double y_min;
    double y_max;
    //strain at lower edge of strain profile
    double eps_lower;
    //strain at upper edge of strain profile
    double eps_upper;
    //strain at origin
    double eps_o;
    //neutral axis
    double d;
    //curvature FI
    double fi;

public:
    StrainProfile(Section *sect, double omg, double el, double eu);
    //setters - all call calculateStrainProfileProperties(void)
    //in order to update Strain Profile
    void setSection(Section *sect);
    void setRotationAngle(double omg);
    void setEpsilonLower(double el);
    void setEpsilonUpper(double eu);
    void setEpsilonLimits(double el, double eu);
    //getters methods
    Section *getSection();
    double getRotationAngle(void);
    double getEpsilonLower(void);
    double getEpsilonUpper(void);

    //methods to get calculated properties of Strain Profile
    double getNeutralAxis_d(void);
    double getCurvature_fi(void);
    double getStrainAtOrigin_eps_o(void);
    double getStrainAtPoint(Point *p);
    double getStrainAtLocalPoint(Point *p);


private:
    void calculateStrainProfileProperties(void);
};

#endif // STRAINPROFILE_H
