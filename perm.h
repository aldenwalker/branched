#ifndef __PERM_H__
#define __PERM_H__

#include <vector>
#include <iostream>

#include "dslist.h"

/******************************************************************************
 * these functions identify a range {smin, ..., -1, 1, ..., smax} with 
 * {dmin, ..., -1, 1, ..., dmax}, where all the things in the smaller 
 * range are associated with things in the larger range, but things in the 
 * larger range might not be glued to anything (they are set to 0)
 *
 * the map and inverse map are just vectors, so there's some annoyance 
 * in figuring out what index in the list stores what value
 * the map_at and imap_at functions do this
 * ***************************************************************************/
struct PDPerm {
  DSList<int> map;
  DSList<int> inverse_map;
  
  PDPerm();
  PDPerm(int smi, int sma, int dmi, int dma);
  
  int smin();
  int smax();
  int dmin();
  int dmax();
  
  int max_size();
  
  void reset(); //reset to the first pdperm
  bool next();  //advance to the next pdperm (returns true if we're at the end)
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