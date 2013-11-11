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
    
    T at(int x) const {
      return L[am+x-(x>0)];
    }
    
    int min() const{
      return m;
    }
    
    int bottom() const {
      if (m == 0) {
        if (M == 0) return 0;
        return 1;
      }
      return m;
    }
    
    int max() const {
      return M;
    }
    
    int top() const {
      if (M == 0) {
        if (m==0) return 0;
        return -1;
      }
      return M;
    }
      
    
    int size() const {
      return M + am;
    }
    
    void reset(const T& val) {
      for (int i=0; i<(int)L.size(); ++i) {
        L[i] = val;
      }
    }
};

template <class T>
std::ostream& operator<<(std::ostream& os, const DSList<T>& d) {
  os << "[";
  for (int i=d.min(); i<0; ++i) {
    os << i << ":" << d.at(i);
    if (i<-1) os << ",";
  }
  os << "_";
  for (int i=1; i<=d.max(); ++i) {
    os << i << ":" << d.at(i);
    if (i<d.max()) os << ",";
  }
  os << "]";
  return os;
}

inline int add_no_zero(int a, int b) {
  int ans = a+b;
  if (ans == 0) return (b>0 ? 1 : -1);
  return ans;
}

#endif