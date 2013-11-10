#include <vector>
#include <iostream>

#include "perm.h"

/******************************************************************************
 * partially defined permuations
 * ****************************************************************************/
PDPerm::PDPerm() {
  smin = smax = asmin = dmin = dmax = admin = 0;
  map.resize(0);
  inverse_map.resize(0);
}

PDPerm::PDPerm(int smi, int sma, int dmi, int dma) {
  smin = smi;
  asmin = abs(smi);
  smax = sma;
  dmin = dmi;
  admin = abs(dmi);
  dmax = dma;
  map = std::vector<int>(smax-smin, 0);
  inverse_map = std::vector<int>(dmax-dmin, 0);
  int map_pos = smin;
  int imap_pos = dmin;
  while (true) {
    if (map_pos == 0) ++map_pos;
    if (imap_pos == 0) ++imap_pos;
    if (map_pos > smax || imap_pos > dmax) break;
    map_at(map_pos) = imap_pos;
    imap_at(imap_pos) = map_pos;
  } 
}

int& PDPerm::map_at(int x) {
  if (x == 0) {
    std::cout << "Error; cannot get map at index 0\n";
    return map[-1];
  }
  return map[asmin+x-(x>0)];
}

int& PDPerm::imap_at(int x) {
  if (x==0) {
    std::cout << "Error; cannot get map at index 0\n";
    return inverse_map[-1];
  }
  return inverse_map[admin+x-(x>0)];
}

int PDPerm::max_size() {
  int ss = asmin + smax;
  int sd = admin + dmax;
  return (ss > sd ? ss : sd);
}

std::ostream& operator<<(std::ostream& os, const PDPerm& p) {
  os << "PDPerm([";
  for (int i=0; i<(int)p.map.size()-1; ++i) {
    os << p.map[i] << ",";
  }
  if (p.map.size()>0) {
    os << p.map[p.map.size()-1];
  }
  os << "])";
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  