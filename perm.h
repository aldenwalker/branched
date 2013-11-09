#ifndef __PERM_H__
#define __PERM_H__

#include <vector>
#include <iostream>

/******************************************************************************
 * these functions identify a range {smin, ..., -1, 1, ..., smax} with 
 * {dmin, ..., -1, 1, ..., dmax}, where all the things in the smaller 
 * range are associated with things in the larger range, but things in the 
 * larger range might not be glued to anything (they are set to 0)
 * ***************************************************************************/
struct PDPerm {
  int smin, smax, szero;
  int dmin, dmax, dzero;
  std::vector<int> map;
  std::vector<int> inverse_map;
  
  PDPerm();
  PDPerm(int smi, int sma, int dmi, int dma);
  
  //int ap(int x);
  //int inv_ap(int x);
};

std::ostream& operator<<(std::ostream& os, const PDPerm& p);

/******************************************************************************
 * This is just a normal permuation class
 * ***************************************************************************/
class Perm {
  private:
    std::vector<int> p;
    std::vector<std::vector<int> > cycles;
    
  public:
    Perm();
    Perm(int n);
    //Perm(const Perm& P);
    Perm(std::vector<int>& targets);
    Perm(int n, std::vector<int>& targets);
    int ap(int x);
    int size();
    Perm inverse();
    void find_cycles();
    int num_cycles();
    int cycle_len(int c);
    int cycle_element(int c, int i);
    //permutations act on the left, so this computes the permutation
    //which is other, followed by this
    Perm operator*(const Perm& other);
    
    friend std::ostream& operator<<(std::ostream& os, const Perm& p);
};


#endif