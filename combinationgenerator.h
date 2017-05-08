#ifndef COMBINATIONGENERATOR_H
#define COMBINATIONGENERATOR_H

#include <vector>

class CombinationGenerator
{
    //maximal number of elements i.e layers, surfaces that can intersect
    static const int MAX_ELEM = 100;
    //helper table to generete combinations
    int tab[MAX_ELEM + 1];
    //n - is a number of elements from which we choose
    //k - is a length of choosing sequences of elements
    int n, k;
    int noOfElements;

    std::vector< std::vector<int> > combinations;
public:
    CombinationGenerator(int noOfElements);
    std::vector< std::vector<int> > getAllCombinations();
private:
    void generateCombinations(void);
    void addCombination(void);
};

#endif // COMBINATIONGENERATOR_H
