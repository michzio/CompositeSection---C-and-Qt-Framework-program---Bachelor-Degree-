#ifndef GAUSSLEGENDRE_H
#define GAUSSLEGENDRE_H

#include<string>
#include<cmath>

struct NotSupportedNumberOfGaussPointsException : public std::exception
{
    std::string s;
    NotSupportedNumberOfGaussPointsException(std::string ss) : s(ss) {}
    ~NotSupportedNumberOfGaussPointsException() throw () {} // Updated
    const char* what() const throw() { return s.c_str(); }
};

struct OutOfRangeGaussPointException : public std::exception
{
    std::string s;
    OutOfRangeGaussPointException(std::string ss) : s(ss) { }
    ~OutOfRangeGaussPointException() throw () {} //Updated
    const char* what() const throw() { return s.c_str(); }
};

class GaussLegendre
{
 private:
    //weight values for 1,2,3,4 or number of Gauss points
    const double weights[5][5] = {
      { 2.0 }, //nG = 1
      {1.0, 1.0 }, //nG = 2
      {5.0/9.0, 8.0/9.0, 5.0/9.0}, //nG = 3
      { (18.0 - std::sqrt(30.0))/36, //nG = 4
        (18.0 + std::sqrt(30.0))/36,
        (18.0 + std::sqrt(30.0))/36,
        (18.0 - std::sqrt(30.0))/36
       },
       { (322.0 - 13.0*std::sqrt(70.0))/900, //nG = 5
         (322.0 + 13.0*std::sqrt(70.0))/900,
         128.0/225.0,
         (322.0 + 13.0*std::sqrt(70.0))/900,
         (322.0 - 13.0*std::sqrt(70.0))/900,
        }
    };
    //lamda values for 1,2,3,4 or 5 number of Gauss points
    const double lambdas[5][5] = {
      { 0 }, //nG = 1
      { -std::sqrt(1.0/3.0), std::sqrt(1.0/3.0) }, //nG = 2
      { -std::sqrt(3.0/5.0), 0, std::sqrt(3.0/5.0) }, //nG = 3
      { -std::sqrt((3.0 + 2.0*std::sqrt(6.0/5.0))/7), //nG = 4
        -std::sqrt((3.0 - 2.0*std::sqrt(6.0/5.0))/7),
         std::sqrt((3.0 - 2.0*std::sqrt(6.0/5.0))/7),
         std::sqrt((3.0 + 2.0*std::sqrt(6.0/5.0))/7),
       },
       {  -(1.0/3.0)*std::sqrt(5.0 + 2.0*std::sqrt(10.0/7.0)), //nG = 5
          -(1.0/3.0)*std::sqrt(5.0 - 2.0*std::sqrt(10.0/7.0)),
           0,
           (1.0/3.0)*std::sqrt(5.0 - 2.0*std::sqrt(10.0/7.0)),
           (1.0/3.0)*std::sqrt(5.0 + 2.0*std::sqrt(10.0/7.0)),
        },

    };

    int nG;

public:
    //constructor takes in number of Gaussian points
    GaussLegendre(int noOfPoints) : nG(noOfPoints) {
        if(nG>5)
            throw NotSupportedNumberOfGaussPointsException("Nie wspierana taka liczba pkt Gaussa.");
    }

    double weight(int m)
    {
        if(m >= nG)
            throw OutOfRangeGaussPointException("Próba pobrania wagi spoza bieżącego zakresu punktów Gaussa.");
        return weights[nG-1][m];
    }

    double lambda(int m)
    {
        if(m >= nG)
            throw OutOfRangeGaussPointException("Próba pobrania lambda_m spoza bieżącego zakresu punktów Gaussa.");
        return lambdas[nG-1][m];
    }
};

#endif // GAUSSLEGENDRE_H
