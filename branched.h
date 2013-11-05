#ifndef __BRANCHED_H__
#define __BRANCHED_H__

#include <vector>
#include <string>
#include <iostream>

#include "point.h"
#include "rational.h"
#include "graphics.h"

typedef int SignedInd;

int sgn(int x);
int pos_mod(int a, int b);

// we could use a standard simplex type, but the boundaries 
// are confusing, so we'll use different types
struct Vertex {
  Vertex();
  std::vector<SignedInd> in_bd_of;
};
std::ostream& operator<<(std::ostream& os, Vertex& v);

struct Edge {
  Edge();
  int start;
  int end;
  bool two_sided;
  Point2d<Rational> start_pos;
  Point2d<Rational> end_pos;
  Point2d<Rational> start_neg;
  Point2d<Rational> end_neg;
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



struct Cell {
  std::vector<SignedInd> bd;
};
std::ostream& operator<<(std::ostream& os, Cell& c);


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
    int chi(bool resolve_vertices=true);
  
    Triangulation branched_surface_from_vector(std::vector<int>& weights);
    Triangulation resolve_branched_surface();
};

struct Cellulation {
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  std::vector<Cell> cells;
  std::vector<std::vector<SignedInd> > loops;
  
  void print(std::ostream& os);
  void draw_to_XGraphics(XGraphics& X);
};

#endif