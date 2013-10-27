#include <vector>
#include <iostream>
#include <string>
#include <stringstream>
#include <cstdlib>

#include "branched.h"

Vertex::Vertex() {
  in_bd_of.resize(0);
}

std::ostream& operator<<(std::ostream& os, Vertex& v) {
  os << "V(); {";
  for (int i=0; i<(int)v.in_bd_of.size(); ++i) {
    os << v.in_bd_of[i];
    if (i<v.in_bd_of.size()-1) os << ", ";
  }
  os << "}";
  return os;
}

Edge::Edge() {
  bd.resize(2,0);
  in_bd_of.resize(0);
}

SignedInd& Edge::operator[](int i) {
  return bd[i];
}

std::ostream& operator<<(std::ostream& os, Edge& e) {
  os << "E(" << e.bd[0] << "," << e.bd[1] << "); {"
  for (int i=0; i<(int)e.in_bd_of.size(); ++i) {
    os << e.in_bd_of[i];
    if (i<e.in_bd_of.size()-1) os << ", ";
  }
  os << "}";
  return os;
}

Triangle::Triangle() {
  bd.resize(3,0);
  in_bd_of.resize(0);
}

SignedInd& Triangle::operator[](int i) {
  return bd[i];
}

std::ostream& operator<<(std::ostream& os, Triangle& t) {
  os << "T(" << t.bd[0] << "," << t.bd[1] << "," << t.bd[2] << "); {"
  for (int i=0; i<(int)t.in_bd_of.size(); ++i) {
    os << t.in_bd_of[i];
    if (i<t.in_bd_of.size()-1) os << ", ";
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
  for (int i=1; i<(int)vertices.size(); ++i) {
    os << i << ": " << vertices[i] << "\n";
  }
  os << "Edges (" << edges.size()-1 << ")\n";
  for (int i=1; i<(int)edges.size(); ++i) {
    os << i << ": " << edges[i] << "\n";
  }
  os << "Triangles (" << triangles.size()-1 << ")\n";
  for (int i=1; i<(int)triangles.size(); ++i) {
    os << i << ": " << triangles[i] << "\n";
  }
}


// Given a integral weight vector on the triangles, produce a 
// branched surface triangulation containing appropriate duplicates 
// of the triangles in the weights.  Note any unused simplices are 
// removed, so the indices of everything can be different
Triangulation Triangulation::branched_surface_from_vector(std::vector<int>& weights) {
  Triangulation new_T;

  //Determine which edges are used, and figure out what the new 
  //indices of the used edges will be
  int num_edges_used = 0;
  int num_vertices_used = 0;
  std::vector<bool> is_edge_used(edges.size(), false);
  for (int i=0; i<(int)weights.size(); ++i) {
    if (weights[i] == 0) continue;
    for (int j=0; j<triangles[i].dim; ++j) {
      is_edge_used[ abs(triangles[i][j]) ] = true;
    }
  }
  std::vector<int> edge_index_translation_table(edges.size(), 0);
  int current_new_index = 1;
  for (int i=1; i<(int)edges.size(); ++i) {
    if (is_edge_used[i]) {
      edge_index_translation_table[i] = current_new_index;
      ++current_new_index;
    }
  }
  num_edges_used = current_new_index-2;
  
  //Determine which vertices are used, and figure out what the new 
  //indices will be
  std::vector<bool> is_vertex_used(vertices.size(), false);
  for (int i=1; i<(int)edges.size(); ++i) {
    if (is_edge_used[i]) {
      is_vertex_used[ abs(edges[i][0]) ] = true;
      is_vertex_used[ abs(edges[i][1]) ] = true;
    }
  }
  std::vector<int> vertex_index_translation_table(vertices.size(), 0);
  current_new_index = 1;
  for (int i=1; i<(int)vertices.size(); ++i) {
    if (is_vertex_used[i]) {
      vertex_index_translation_table[i] = current_new_index;
      ++current_new_index;
    }
  }
  num_vertices_used = current_new_index-2;
  
  //Create the new lists of vertices and edges
  //record which edges are incident to which vertices
  new_T.vertices.resize(num_vertices_used+1, Simplex(0));
  new_T.edges.resize(num_edges_used+1, Simplex(1));
  for (int i=1; i<(int)edges.size(); ++i) {
    int ei = edge_index_translation_table[i];
    if (ei == 0) continue;
    int v0i = vertex_index_translation_table[ edges[i][0] ];
    new_T.edges[i][0] = v0i;
    new_T.vertices[v0i].in_bd_of.push_back(ei);
    int v1i = vertex_index_translation_table[ edges[i][1] ];
    new_T.edges[i][1] = v1i;
    new_T.vertices[v1i].in_bd_of.push_back(-ei);
  }
  
  //Create the new list of triangles
  new_T.triangles.resize(1);
  for (int i=0; i<(int)weights.size(); ++i) {
    Simplex temp_simp(2);
    int ei[3];
    for (int j=0; j<3; ++j) {
      ei[j] = sgn(triangles[i][j])*edge_index_translation_table[ abs(triangles[i][j]) ];
      temp_simp[j] = ei[j];
    }
    for (int j=0; j<weights[i]; ++j) {
      int this_t_ind = new_T.triangles.size();
      new_T.edges[ 
      
  
  
  
  
  
  
}

Triangulation Triangulation::resolve_branched_surface() {

}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


