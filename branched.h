#ifndef __BRANCHED_H__
#define __BRANCHED_H__

#include <vector>
#include <string>
#include <iostream>


typedef SignedInd int;

int sgn(int x) { 
  return ((x > 0) - (x < 0));
}

// we could use a standard simplex type, but the boundaries 
// are confusing, so we'll use different types
struct Vertex {
  Vertex();
  std::vector<SignedInd> in_bd_of;
};
std::ostream& operator<<(std::ostream& os, Vertex& v);

struct Edge {
  Edge();
  std::vector<SignedInd> bd;
  std::vector<SignedInd> in_bd_of;
  SignedInd& operator[](int i);
};
std::ostream& operator<<(std::ostream& os, Edge& e);

struct Triangle {
  Triangle();
  std::vector<SignedInd> bd;
  std::vector<SignedInd> in_bd_of;
  SignedInd& operator[](int i);
};
std::ostream& operator<<(std::ostream& os, Triangle& t);


/*****************************************************************************
 * Triangulation                                                             
 *****************************************************************************/
//all lists in here are ONE-BASED, so that SignedInds are easier
struct Triangulation {
    std::vector<Simplex> triangles;
    std::vector<Simplex> edges;
    std::vector<Simplex> vertices;

    Triangulation();
    void read_file(std::string filename);
    void write_file(std::string filename);
    void print(ostream& os);
  
    Triangulation branched_surface_from_vector(std::vector<int>& weights);
    Triangulation resolve_branched_surface();
};

#endif