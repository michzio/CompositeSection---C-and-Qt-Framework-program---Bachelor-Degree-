#ifndef MATERIAL_H
#define MATERIAL_H

#include<string>
#include <QColor>

enum MaterialType { CONCRETE, STRUCTURAL_STEEL, REINFORCEMENT, REINFORCEMENT_BARS, FAFITIS_CONCRETE, FAFITIS_REINFORCEMENT, HOLLOW, FIBER_REINFORCED_PLASTIC};

class Material
{
    //material properties
    MaterialType materialType;
public:
    Material(MaterialType);
    std::string getMaterialName(void);
    MaterialType getMaterialType(void);
    //stress calculation methods
    virtual double sigmaEpsilon(double epsilon) = 0;

    //virtual methods, must be overriden in subclases of Material class
    virtual int getModuleOfElasticity_E() = 0;
    virtual int numberOfFunctionParts(void) = 0;
    virtual int functionPartNumberFor(double epsilon) = 0;
    virtual double functionPartStrainLimitA(int number) = 0;
    virtual double functionPartStrainLimitB(int number) = 0;
    virtual double functionPartStrainLimitA(double epsilon) = 0;
    virtual double functionPartStrainLimitB(double epsilon) = 0;
    virtual int functionPartDegree(int number) = 0;
    virtual int functionPartDegree(double epsilon) = 0;

    //material tensile and compressive strain limits
    virtual double tensileStrainLimit_eps_tu(void) = 0;
    virtual double compressiveStrainLimit_eps_cu(void) = 0;
    virtual double pureCompressionStrainLimit_eps_co(void) = 0;
    virtual Qt::GlobalColor materialColor(void) = 0;

};

//Exception thrown when user has choosen not existant material
struct MaterialNotFoundException : public std::exception
{
   std::string s;
   MaterialNotFoundException(std::string ss) : s(ss) {}
   ~MaterialNotFoundException() throw () {} // Updated
   const char* what() const throw() { return s.c_str(); }
};

//Exception thrown when user pass strain epsilon value out of bounds
struct StrainOutOfBoundsException : public std::exception
{
   std::string s;
   StrainOutOfBoundsException(std::string ss) : s(ss) {}
   ~StrainOutOfBoundsException() throw () {} // Updated
   const char* what() const throw() { return s.c_str(); }
};
#endif // MATERIAL_H
