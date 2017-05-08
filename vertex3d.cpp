#include "vertex3d.h"
#include "math.h"

Vertex3D::Vertex3D(double N, double Mx, double My)
{
    setN(N); setMx(Mx); setMy(My);
}

//setter methods
void Vertex3D::setN(double N)
{
    //rounding and casting to int
    z = (int) (N + 0.5);
}

void Vertex3D::setMx(double Mx)
{
    //rounding and casting to int
    x = (int) (Mx + 0.5);
}

void Vertex3D::setMy(double My)
{
    //rounding and casting to int
    y = (int) (My + 0.5);
}

//getter methods
int Vertex3D::N(void) const
{
    return z;
}

int Vertex3D::Mx(void) const
{
    return x;
}

int Vertex3D::My(void) const
{
    return y;
}

int Vertex3D::M(void) const
{
    //calculating M value
    int M = (int) (sqrt(pow(Mx(),2.0) + pow(My(),2.0)) + 0.5);
    //adjusting sign of M in new 2D M-N coordiate system
    if(My() < 0) M = -M;
    else if(My() == 0 && Mx() < 0) M = -M;
    //fprintf(stderr, "M = %d, Mx = %d, My = %d\n", M, Mx(), My());

    return M;
}

double Vertex3D::argumentMxMy(void) const
{
    return argument( Mx(), My() );
}

double Vertex3D::argumentMyN(void) const
{
    return argument( My(), N());
}

double Vertex3D::argumentMxN(void) const
{
    return argument( Mx(), N() );
}

double Vertex3D::argumentMN(void) const
{
    return argument( M(), N() );
}

/**
 * @brief Vertex3D::argument
 * @param a - abscissa
 * @param b - ordinata
 * @return value of argument calculated with below formula
 * atan(b/a) for a > 0
 * atan(b/a) + PI for a < 0 and b >= 0
 * atan(b/a) - PI for a < 0 and b < 0
 * +PI/2 for a == 0 and b > 0
 * -PI/2 for a == 0 and b < 0
 * 0 for a == 0 and b == 0
 */
double Vertex3D::argument(double a, double b) const
{

    double fi = atan2(b,a);

    //fprintf(stderr, "a: %g b: %g => fi: %g\n", a, b, fi);

    if(fi >= 0)
        return fi;

    return fi + 2*M_PI;

    /* calculated using formula with conditions
     double fi = 0.0;
    if(a>0)
        fi = atan(b/a);
    if(a < 0 && b >= 0)
        fi = atan(b/a) + M_PI;
    if(a < 0 && b < 0)
        fi = atan(b/a) - M_PI;
    if(a == 0 && b > 0)
        fi =  M_PI/2;
    if(a == 0 && b < 0)
        fi = -M_PI/2;

    if(fi >= 0)
        return fi;

     return fi + 2*M_PI;
    */
}
