#ifndef __BRANCHED_H__
#define __BRANCHED_H__

#include <vector>
#include <string>
#include <iostream>


typedef int SignedInd;

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
  std::vector<SignedInd> in_bd_pos;
  std::vector<SignedInd> in_bd_neg;
  SignedInd& operator[](int i);
};
std::ostream& operator<<(std::ostream& os, Edge& e);

struct Triangle {
  Triangle();
  std::vector<SignedInd> bd;
  SignedInd& operator[](int i);
};
std::ostream& operator<<(std::ostream& os, Triangle& t);


/*****************************************************************************
 * Triangulation                                                             
 *****************************************************************************/
//all lists in here are ONE-BASED, so that SignedInds are easier
struct Triangulation {
    std::vector<Triangle> triangles;
    std::vector<Edge> edges;
    std::vector<Vertex> vertices;
    std::vector<std::vector<SignedInd> > fundamental_loops;

    Triangulation();
    int add_edge(int v0, int v1);
    int add_triangle(SignedInd e0, SignedInd e1, SignedInd e2);
    void set_closed_surface(int genus);
    void read_file(std::string filename);
    void write_file(std::string filename);
    void print(std::ostream& os);
  
    Triangulation branched_surface_from_vector(std::vector<int>& weights);
    Triangulation resolve_branched_surface();
};

#endif