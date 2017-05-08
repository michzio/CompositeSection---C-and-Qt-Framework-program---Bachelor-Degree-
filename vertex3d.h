#ifndef VERTEX3D_H
#define VERTEX3D_H

#include "stdio.h"
#include "math.h"

class Vertex3D
{
    int x;
    int y;
    int z;

public:

    Vertex3D(double N, double Mx, double My);
    //setter methods
    void setN(double);
    void setMx(double);
    void setMy(double);

    //getter methods
    int N(void) const;
    int Mx(void) const;
    int My(void) const;
    int M(void) const;

    //calculates argument (fi angle) in indicated plane of coordinate system
    double argumentMxMy(void) const;
    double argumentMyN(void) const;
    double argumentMxN(void) const;
    double argumentMN(void) const;

private:
    //calculatrs argument (fi angle) in indicated plane of coordinate system
    double argument(double a, double b) const;
};

//Vertex3D comparators used by qSort() function in sorting
//containers containing Vertex3D objects as elements
struct Vertex3DComparatorByN
{
    bool operator()(const Vertex3D *a, const Vertex3D *b) const
    {
        //if N values of vertex a and vertex b are equal
        //then we compare arguments (fi angles) of vertices
        //in plane defined by Mx and My axis
        //values with lower angle are less then values with greater angle
        //(angle values are in range 0 to 360 degrees or 0 - 2PI radians)
        if(a->N() == b->N())
        {
            //checking angle value
            return a->argumentMxMy() < b->argumentMxMy();
        }

        return a->N() < b->N();
    }
};

struct Vertex3DComparatorByMx
{
    bool operator()(const Vertex3D *a, const Vertex3D *b) const
    {
        //if Mx values of vertex a and vertex b are equal
        //then we compare arguments (fi angles) of vertices
        //in plane defined by N and My axis
        //values with lower angle are less then values with greater angle
        //(angle values are in range 0 to 360 degrees or 0 - 2PI radians)
        if(a->Mx() == b->Mx())
        {
            //checking angle value
            return a->argumentMyN() < b->argumentMyN();
        }

        return a->Mx() < b->Mx();
    }
};

struct Vertex3DComparatorByAlfa
{
    bool operator()(const Vertex3D *a, const Vertex3D *b) const
    {
        //we firstly sort vertices by angle alfa which is
        //angle of moment M inclination M = sqrt(Mx^2 + My^2)
        //then we sort them by angle i M-N coordinate system
        //i.e. Beta equal to angle arctan(N/M)

        //converstion alfa to degree value! only whole integer angles taken into account!
        //values in range [180, 360] are transformmed to [0,180]
        //i.e. vertices with angle alfa 45 degree and 225 degree are treated as laying in
        //the same plane M-N
        int aAlfa = (int)(a->argumentMxMy() *180/M_PI)%180;
        int bAlfa = (int)(b->argumentMxMy() *180/M_PI)%180;
        if(aAlfa == bAlfa)
        {

            //checking angle value
            return a->argumentMN() < b->argumentMN();
        }

        return aAlfa < bAlfa;
    }
};

struct Vertex3DComparatorByMy
{
    bool operator()(const Vertex3D *a, const Vertex3D *b) const
    {
        //if My values of vertex a and vertex b are equal
        //then we compare arguments (fi angles) of vertices
        //in plane defined by N and Mx axis
        //values with lower angle are less then values with greater angle
        //(angle values are in range 0 to 360 degrees or 0 - 2PI radians)
        if(a->My() == b->My())
        {
            //checking angle value
            return a->argumentMxN() < b->argumentMxN();
        }

        return a->My() < b->My();
    }
};

struct Vertex3DComparatorByAlfaInMxMy
{
    bool operator()(const Vertex3D *a, const Vertex3D *b) const
    {
         return a->argumentMxMy() < b->argumentMxMy();
    }
};

struct Vertex3DComparatorByAlfaInMyN
{
    bool operator()(const Vertex3D *a, const Vertex3D *b) const
    {
         return a->argumentMyN() < b->argumentMyN();
    }
};

struct Vertex3DComparatorByAlfaInMxN
{
    bool operator()(const Vertex3D *a, const Vertex3D *b) const
    {
         return a->argumentMxN() < b->argumentMxN();
    }
};

struct Vertex3DComparatorByBetaInMN
{
    bool operator()(const Vertex3D *a, const Vertex3D *b) const
    {
        return a->argumentMN() < b->argumentMN();
    }
};

#endif // VERTEX3D_H
