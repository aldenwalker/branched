#include <vector>
#include <iostream>
#include <string>
#include <stringstream>

#include "branched.h"


/****************************************************************************
 * Simplex                                                                  
 ****************************************************************************/
Simplex::Simplex() {
  dim = 0;
  bd.resize(0);
  in_bd_of.resize(0);
}

Simplex::Simplex(int n) {
  dim = n;
  bd = std::vector<SignedInd>(n, 0);
  in_bd_od.resize(0);
}

Simplex::Simplex(std::vector<SignedInd>& bd) {
  dim = bd.size();
  this->bd = bd;
  in_bd_of.resize(0);
}

std::string Simplex::ihat_string(int i) {
  std::ostringstream ans;
  for (int j=0; j<dim; ++j) {
    if (j==i) continue;
    ans << j;
  }
  return ans.string();
}

int Simplex::dimension() {
  return dim;
}

SignedInd& Simplex::boundary(int i) {
  return bd[i];
}

SignedInd& Simplex::operator[](int i) {
  return bd[i];
}

SignedInd& Simplex::in_boundary_of(int i) {
  return in_bd_of[i];
}

std::ostream& operator<<(std::ostream& os, Simplex& S) {
  os << S.dim << "s{";
  for (int i=0; i<dim; ++i) {
    os << ihat_string(i) << ": " << bd[i];
    if (i<dim) os << ", ";
  }
  os << "}";
  return os;
}

/*****************************************************************************
 * Triangulation                                                             
 *****************************************************************************/
Triangulation::Triangulation() {
  triangles.resize(0);
  edges.resize(0);
  vertices.resize(0);
}

void Triangulation::read_file(std::string filename) {}

void Triangulation::write_file(std::string filename) {}

void Triangulation::print(ostream& os) {
  os << "Vertices (" << vertices.size()-1 << "):\n";
  for (int i=1; i<=vertices.size(); ++i) {
    os << i << ": " << vertices[i] << "\n";
  }
  os << "Edges (" <<edges.size()-1 << ")\n";
  for (int i=1; i<=edges.size(); ++i) {
    os << i << ": " << edges[i] << "\n";
  }
  os << "Triangles (" <<triangles.size()-1 << ")\n";
  for (int i=1; i<=triangles.size(); ++i) {
    os << i << ": " << triangles[i] << "\n";
  }
}

Triangulation Triangulation::branched_surface_from_vector(std::vector<int>& weights) {
  
}

Triangulation Triangulation::resolve_branched_surface() {

}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


