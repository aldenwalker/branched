#ifndef __BRANCHED_H__
#define __BRANCHED_H__

#include <vector>
#include <string>
#include <iostream>

#include "point.h"
#include "rational.h"
#include "graphics.h"


//forward declaration of looparrangement
struct LoopArrangement;

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
  bool boundary_loop;
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
  int sign;
  bool contains_boundary;
  bool computed_winding_number;
  int winding_number;
  std::vector<SignedInd> bd;
  Point2d<Rational> coords;
};
std::ostream& operator<<(std::ostream& os, Cell& c);


struct Cellulation {
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  std::vector<Cell> cells;
  std::vector<std::vector<SignedInd> > loops;
  
  void print(std::ostream& os);
  SignedInd next_edge(SignedInd e);
  std::vector<SignedInd> follow_edge(SignedInd e);
  void draw_to_xgraphics(XGraphics& X);
  int verbose;
  
  //this computes the winding number of the loops around a given cell, relative 
  //to another cell
  void compute_winding_numbers(LoopArrangement& LA, int relative_to_cell=0);
  
  //this computes a naive lower bound on negative euler characteristic of a surface 
  //with the desired boundary
  int chi_upper_bound(LoopArrangement& LA);
  
  
};

#endif