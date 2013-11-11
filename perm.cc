#include <vector>
#include <iostream>

#include "perm.h"

/******************************************************************************
 * partially defined permuations
 * ****************************************************************************/
PDPerm::PDPerm() {
  map = DSList<int>();
  inverse_map = DSList<int>();
}

PDPerm::PDPerm(int smi, int sma, int dmi, int dma) {
  map = DSList<int>(smi, sma, 0);
  inverse_map = DSList<int>(dmi, dma, 0);
  int map_pos = map.min();
  int imap_pos = inverse_map.min();
  while (true) {
    if (map_pos == 0) ++map_pos;
    if (imap_pos == 0) ++imap_pos;
    if (map_pos > map.max() || imap_pos > inverse_map.max()) break;
    map[map_pos] = imap_pos;
    inverse_map[imap_pos] = map_pos;
    ++map_pos;
    ++imap_pos;
  } 
}

int PDPerm::smin() {
  return map.min();
}

int PDPerm::smax() {
  return map.max();
}

int PDPerm::dmin() {
  return inverse_map.min();
}

int PDPerm::dmax() {
  return inverse_map.max();
}

int PDPerm::max_size() {
  int ss = map.size();
  int sd = inverse_map.size();
  return (ss > sd ? ss : sd);
}

void PDPerm::reset() {
  map.reset(0);
  inverse_map.reset(0);
  int map_pos = map.min();
  int imap_pos = inverse_map.min();
  while (true) {
    if (map_pos == 0) ++map_pos;
    if (imap_pos == 0) ++imap_pos;
    if (map_pos > map.max() || imap_pos > inverse_map.max()) break;
    map[map_pos] = imap_pos;
    inverse_map[imap_pos] = map_pos;
    ++map_pos;
    ++imap_pos;
  } 
}


//this returns true if a >=  b, except where we consider 
//zero to be bigger than everything
bool ge_zero_big(int a, int b) {
  if (a==0) return true;
  if (b==0) return false;
  return a >= b;
}

//this does the actual next finding
bool next_dslist(DSList<int>& L) {
  int T = L.top();
  int B = L.bottom();
  //find the first time that the value decreases
  int j = T;
  int i = add_no_zero(T, -1);
  while (i >= B && ge_zero_big(L[i], L[j])) {
    j = i;
    i = add_no_zero(i, -1);
  }
  if (i < B) return true;
  //find the smallest thing forward which is bigger
  int k = add_no_zero(j, 1);
  while (k <= T && ge_zero_big(L[k], L[i])) {
    j = k;
    k = add_no_zero(k,1);
  }
  //at this point L[j] needs to be swapped with L[i]
  int temp = L[i];
  L[i] = L[j];
  L[j] = temp;
  //and finally, we need to reverse L[i+1] through the end
  j = add_no_zero(i,1);
  k=T;
  while (j<k) {
    temp = L[j];
    L[j] = L[k];
    L[k] = temp;
    j = add_no_zero(j, 1);
    k = add_no_zero(k,-1);
  }
  return false;
}

bool PDPerm::next() {
  //we'll just use lexicographic ordering, where 0 is the largest
  //it's easier to advance the larger list (which has zeros)
  if (map.size() > inverse_map.size()) {
    bool at_end = next_dslist(map);
    if (at_end) return true;
    for (int i=map.min(); i<=(int)map.max(); ++i) {
      if (i == 0 || map[i] == 0) continue;
      inverse_map[map[i]] = i;
    }
  } else {
    bool at_end = next_dslist(inverse_map);
    if (at_end) return true;
    for (int i=inverse_map.min(); i<=(int)inverse_map.max(); ++i) {
      if (i == 0 || inverse_map[i] == 0) continue;
      map[inverse_map[i]] = i;
    }
  }
  return false;
}

std::ostream& operator<<(std::ostream& os, const PDPerm& p) {
  os << "PDPerm(" << p.map << "," << p.inverse_map << ")";
  return os;
}

/*****************************************************************************
 * usual permuations
 * ***************************************************************************/
Perm::Perm() {
  p.resize(0);
  cycles.resize(0);
}

Perm::Perm(int n) {
  p.resize(n);
  cycles.resize(0);
  for (int i=0; i<n; ++i) {
    p[i] = i;
  }
}

//Perm::Perm(const Perm& P) {
//  p.resize(0);
//  cycles.resize(0);
  //p = P.p;
  //cycles = P.cycles;
//}

Perm::Perm(std::vector<int>& targets) {
  int n = (int)targets.size();
  p.resize(n);
  cycles.resize(0);
  for (int i=0; i<n; ++i) {
    p[i] = targets[i];
  }
}

Perm::Perm(int n, std::vector<int>& targets) {
  p.resize(n);
  cycles.resize(0);
  for (int i=0; i<n; ++i) {
    p[i] = targets[i];
  }
}

int Perm::ap(int x) {
  return p.at(x);
}

int Perm::size() {
  return (int)p.size();
}

Perm Perm::inverse() {
  int n = (int)p.size();
  std::vector<int> t(n);
  for (int i=0; i<n; ++i) {
    t[ p[i] ] = i;
  }
  return Perm(t);
}

void Perm::find_cycles() {
  int n = (int)p.size();
  std::vector<bool> done(n, false);
  cycles.resize(0);
  for (int i=0; i<n; ++i) {
    if (done[i]) continue;
    std::vector<int> temp_cycle(0);
    int j=i;
    do {
      temp_cycle.push_back(j);
      done[j] = true;
      j = p[j];
    } while (j != i);
    cycles.push_back(temp_cycle);
  }
}

int Perm::num_cycles() {
  if (cycles.size() != 0) {
    return (int)cycles.size();
  }
  this->find_cycles();
  return (int)cycles.size();
}

int Perm::cycle_len(int c) {
  if (cycles.size() != 0) {
    return (int)cycles[c].size();
  }
  this->find_cycles();
  return (int)cycles[c].size();
}

int Perm::cycle_element(int c, int i) {
  if (cycles.size() != 0) {
    return cycles[c][i];
  }
  this->find_cycles();
  return cycles[c][i];
}


Perm Perm::operator*(const Perm& other) {
  Perm np(p.size());
  for (int i=0; i<(int)p.size(); ++i) {
    np.p[i] = p[other.p[i]];
  }
  return np;
}
  


std::ostream& operator<<(std::ostream& os, const Perm& p) {
  if (p.p.size()==0) {
    os << "[]";
    return os;
  }
  os << "[" << p.p[0];
  for (int i=1; i<(int)p.p.size(); ++i) {
    os << " " << p.p[i];
  }
  os << "]";
  return os;
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  