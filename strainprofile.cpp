#include "strainprofile.h"
#include "flags.h"

StrainProfile::StrainProfile(Section *sect,
                             double omg,
                             double el,
                             double eu) : section(sect), omega(omg), eps_lower(el), eps_upper(eu)
{
    calculateStrainProfileProperties();
}

void StrainProfile::calculateStrainProfileProperties(void)
{
    //y_min, y_max
    //setting y_min, y_max vertical limits of strain profile
    Point *centroid = section->getGeometricalCenter();
    y_min = section->yMinInCoordinateSystem(centroid, omega);
    y_max = section->yMaxInCoordinateSystem(centroid, omega);

    //CURVATURE fi
    //strain difference at lower and upper edge
    //esp_lower assume tension +, eps_upper assume compression -
    double dx = eps_lower - eps_upper;
    //vertical length of strain profile
    double dy = y_max - y_min;
    //calculating and setting curvature
    fi = std::atan(dx/dy);

    //NEUTRAL AXIS d
    double dx_cmpr = 0.0 - eps_upper;
    double dx_tens = eps_lower - 0.0;

    //calculating distance from y_max (upper) to neutral axis
    double dy_neutral_upper = dx_cmpr/std::tan(fi);
    //calculating distance form y_min (lower) to nerutral axis
    double dy_neutral_lower = dx_tens/std::tan(fi);

    d = y_max - dy_neutral_upper;

    if(d != (y_min + dy_neutral_lower)) {
        if(DEBUG)
            fprintf(stderr, "Bład podczas wyznaczania polozenia osi obojetnej w profilu odkształceń.\n 1 opcja: %g, 2 opcja: %g.\n", d, (y_min + dy_neutral_lower) );
    }

    // STRAIN AT ORIGIN
    // eps_o = d*tan(fi) for small angles fi eps_o = d*fi
    //calculating strain at origin
    if(fi == 0 && std::isinf(d) && eps_upper == eps_lower) {
        //we have square strain profile with constatnt strain value on the height of strain profile
        //d lies in infinity, eps_upper == eps_lower
        eps_o = eps_upper; //or eps_lower
    } else {
         eps_o = d*std::tan(fi);
    }

    //log all calculated informations about profile
    if(DEBUG) {
        fprintf(stderr, "\n################## CREATING NEW STRAIN PROFILE ###################\n");
        fprintf(stderr, "Strain Profile: omega = %g, centroid = {%g, %g}, y_min = %g, y_max = %g, eps_lower = %g, eps_upper = %g, eps_o = %g, d = %g, fi = %g\n",
            omega, centroid->getX(), centroid->getY(), y_min, y_max, eps_lower, eps_upper, eps_o, d, fi);
    }
}


void StrainProfile::setSection(Section *sect)
{
    section = sect;
    calculateStrainProfileProperties();
}

void StrainProfile::setRotationAngle(double omg)
{
    omega = omg;
    calculateStrainProfileProperties();
}

void StrainProfile::setEpsilonLower(double el)
{
    eps_lower = el;
    calculateStrainProfileProperties();
}

void StrainProfile::setEpsilonUpper(double eu)
{
    eps_upper = eu;
    calculateStrainProfileProperties();
}

void StrainProfile::setEpsilonLimits(double el, double eu)
{
    eps_lower = el;
    eps_upper = eu;
    calculateStrainProfileProperties();
}

Section *StrainProfile::getSection(void)
{
    return section;
}

double StrainProfile::getRotationAngle(void)
{
    return omega;
}

double StrainProfile::getEpsilonLower(void)
{
    return eps_lower;
}

double StrainProfile::getEpsilonUpper(void)
{
    return eps_upper;
}

double StrainProfile::getNeutralAxis_d(void)
{
    return d;
}

double StrainProfile::getCurvature_fi(void)
{
    return fi;
}

double StrainProfile::getStrainAtOrigin_eps_o(void)
{
    return eps_o;
}

double StrainProfile::getStrainAtPoint(Point *p)
{
    //we are getting y coordinate of point p
    //1. setting Point *p coordinate system properly
    Point *adjustedP = new Point(p->getX(), p->getY());

    adjustedP->setOrigin(section->getGeometricalCenter());
    adjustedP->setRotation(omega);

    //2. geting y coordinate as p->getLocalY();
    double y = adjustedP->getLocalY();
    delete adjustedP;

    //Calculating strain at point p
    // (d-y) is distance from neutral axis!
    double eps = (d - y)*std::tan(fi);

    return eps;
}

//now Point *p is given in current local coordinate system
double StrainProfile::getStrainAtLocalPoint(Point *p)
{
    double y = p->getLocalY();

    //calculating strain at point p
    // (d-y) is distance from neutral axis!
    double eps = (d - y)*std::tan(fi);

    return eps;
}

