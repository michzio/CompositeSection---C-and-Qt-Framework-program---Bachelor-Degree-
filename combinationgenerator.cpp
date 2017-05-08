#include "combinationgenerator.h"

CombinationGenerator::CombinationGenerator(int no) : noOfElements(no) {
    //we choose k no of elements from set of n no of elements
    for(k=1; k<noOfElements; k++) {

            for(n=k; n<noOfElements; n++) {
              // printf("A%d", n+1);
               //printf("---------\n");
               generateCombinations();
        }
       }

}

/**
 * @brief CombinationGenerator::generateCombinations
 * function is generating k elements combinations from
 * current n elements set
 */
void CombinationGenerator::generateCombinations(void)
{
    int i, r = k;
    for(i=1;i<=k;i++)
              tab[i] = i;
    addCombination();
    while(r >= 1) {
      if(tab[r] < (n - k + 1) + (r - 1)) {
          tab[r]++;
          for(i=r+1;i<=k;i++)
            tab[i] = tab[i - 1] + 1;
          addCombination();
          r = k;
      }else r--;
    }
}

/**
 * @brief CombinationGenerator::addCombination
 * adds currently generated combination to multidimensional vector
 * table.
 */
void CombinationGenerator::addCombination(void)
{
    //initialize new row with k+1 elements for adding next combination
    std::vector<int> row(k+1);
    //firstly we add n+1 as first element in each row
    row[0] = n+1;

    //add to row all generated elements in current combination
    for(int i=1;i<=k;i++) {
        row[i] = tab[i];
      //printf("%d ", tab[i]);
    }
    combinations.push_back(row);
}

/***
 * getAllCombinations() is a function that
 * generates multidimensional table (vector) with
 * all possible combinations of 1,2,3,4,....n elements from
 * the set of 1,2,3,...n elements ex  2 elements from 3 elements
 * each row in table is preceeded by number representng
 * numbers of elements in set from which sequences are genereated
 * + one i.e. (noOfElements + 1)
 * purpose of this function is to return multidimensional array
 * with indices of surfaces (layers) between which we need to
 * check overlapping while calculating A*, Sx*, Sy*
 * where * means weighted
 */
std::vector< std::vector<int> > CombinationGenerator::getAllCombinations()
{
   return combinations;
}

