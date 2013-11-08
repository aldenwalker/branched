#ifndef PERM_H
#define PERM_H

#include <vector>
#include <iostream>

/******************************************************************************
 * This class is for functions which identify a set with a larger set
 * i.e. injective functions [n]->[m], where m>=n.  Or, equivalently, 
 * partially-defined injective surjective functions [m]->[n]
 * ***************************************************************************/
struct PDPerm {
  int source_size;
  int dest_size;
  std::vector<int> map;
  std::vector<int> inverse_map;
  
  PDPerm();
  PDPerm(int n, int m);
  PDPerm(int n, int m, const std::vector<int>& map, bool is_inverse_map=false);
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