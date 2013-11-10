#ifndef __DSLIST_H__
#define __DSLIST_H__

#include <vector>
#include <iostream>
#include <cstdlib>
#include <stdexcept>

template <class T>
class DSList {
  private:
    int m, M, am;
    std::vector<T> L;
  public:
    DSList() {
      m = am = M = 0;
      L.resize(0);
    }
    
    DSList(int mi, int ma) {
      m = mi;
      am = abs(mi);
      M = ma;
      L.resize(am+M);
    }
    
    DSList(int mi, int ma, const T& val) {
      m = mi;
      am = abs(mi);
      M = ma;
      L = std::vector<T>(am+M, val);
    }
    
    typename std::vector<T>::reference operator[](int x) {
      if (x == 0) {
        throw std::out_of_range("Cannot access DSList at index 0");
      }
      return L[am+x-(x>0)];
    }
    
    int min() {
      return m;
    }
    
    int max() {
      return M;
    }
};


#endif