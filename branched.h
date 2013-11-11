#ifndef __BRANCHED_H__
#define __BRANCHED_H__

#include <vector>
#include <string>
#include <iostream>

#include "point.h"
#include "rational.h"
#include "graphics.h"
#include "perm.h"
#include "dslist.h"


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
std::ostream& operator<<(std::ostream& os, const Edge& e);

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


struct BranchedSurface {
  const Cellulation* C;
  std::vector<std::pair<int, int> > cell_coefficients;
  bool eperms_valid;
  int verbose;
  
  //the edge pdperms record the pdperm that is applied as we go 
  //from the negative side to the positive side
  //i.e. it is the pdperm applied as we rotate around the vertex 
  //counterclockwise and the edge is pointing out
  std::vector<PDPerm> edge_pdperms;
  
  BranchedSurface(const Cellulation* const C, int verbose=1);
  BranchedSurface(const Cellulation* const C, const std::vector<std::pair<int, int> >& cc, int verbose=1);
  
  //initialize the edge perms to something arbitrary
  void init_edge_pdperms();
  
  //compute the euler characteristic of a gluing
  int chi();
  
  //compute a bound on a partially defined chi
  //it'll only compute with vertices that are determined by the
  //edges in edge_is_set
  int partially_defined_chi(const std::vector<bool>& edge_is_set);
  
  //technical function to follow stuff around and record things
  //if direction is 1, it means that we cross all edges going 
  //right to left in the cyclic order on the vertex 
  //(so if the edge is incoming, it means we apply the PDPerm inverse)
  void follow_gluing_around_vertex(const Vertex& vert, int start_edge, int start_level, 
                                   int direction, std::vector<DSList<bool> >& edges_visited, 
                                   int& boundary);
  
  //the only issue in computing the euler characteristic is computing the 
  //number of vertices.  this function computes how many vertices live 
  //over the given vertex
  int num_vertices_over_vertex(int vi, bool quiet=false);
  
  //this optimizes over all gluings to obtain the best possible
  int brute_minimal_gluing();
  
  //this tries to improve the gluing locally
  int hillclimb_minimal_gluing();
  
  //just try to guess a minimal gluing
  //this function can be randomized with the seed
  int guess_minimal_gluing(int seed=0);
  
  //print the data as text
  void print(std::ostream& os);
};

//get the sum of a pair
int sum(std::pair<int, int>& p);

//print a pair
std::ostream& operator<<(std::ostream& os, std::pair<int, int>& p);

#endif