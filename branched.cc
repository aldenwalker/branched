#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>

#include "branched.h"
#include "surface.h"

//return +/-1 depending on the sign
int sgn(int x) { 
  return ((x > 0) - (x < 0));
}

//this is a handy function to get the generator index from the tile position
SignedInd free_gen_from_tile_index(int ind) {
  int q = ind/4;
  int r4 = ind%4;
  int r2 = ind%2;
  int gen_index = (2*q) + r2;
  int gen_index_1_based = gen_index + 1;
  int sign = (r4 < 2 ? 1 : -1);
  return sign * gen_index_1_based;
}

//mod, but always return a nonnegative number
int pos_mod(int a, int b) {
  return (a%b + b)%b;
}


/*****************************************************************************
 * Vertex
 *****************************************************************************/
Vertex::Vertex() {
  in_bd_of.resize(0);
}

std::ostream& operator<<(std::ostream& os, Vertex& v) {
  os << "V(); {";
  for (int i=0; i<(int)v.in_bd_of.size(); ++i) {
    os << v.in_bd_of[i];
    if (i<(int)v.in_bd_of.size()-1) os << ", ";
  }
  os << "}";
  return os;
}

/*****************************************************************************
 * Edge
 *****************************************************************************/
Edge::Edge() {
  bd.resize(2,0);
  in_bd_pos.resize(0);
  in_bd_neg.resize(0);
}

SignedInd& Edge::operator[](int i) {
  return bd[i];
}

std::ostream& operator<<(std::ostream& os, Edge& e) {
  os << "E(" << e.bd[0] << "," << e.bd[1] << "); {{";
  for (int i=0; i<(int)e.in_bd_pos.size(); ++i) {
    os << e.in_bd_pos[i];
    if (i<(int)e.in_bd_pos.size()-1) os << ", ";
  }
  os << "},{";
  for (int i=0; i<(int)e.in_bd_neg.size(); ++i) {
    os << e.in_bd_neg[i];
    if (i<(int)e.in_bd_neg.size()-1) os << ", ";
  }
  os << "}}";
  return os;
}

/*****************************************************************************
 * Triangle
 *****************************************************************************/
Triangle::Triangle() {
  bd.resize(3,0);
}

SignedInd& Triangle::operator[](int i) {
  return bd[i];
}

std::ostream& operator<<(std::ostream& os, Triangle& t) {
  os << "T(" << t.bd[0] << "," << t.bd[1] << "," << t.bd[2] << ");";
  return os;
}


/*****************************************************************************
 * Triangulation                                                             
 *****************************************************************************/
Triangulation::Triangulation() {
  triangles.resize(0);
  edges.resize(0);
  vertices.resize(0);
  fundamental_loops.resize(0);
}

int Triangulation::add_edge(int v0, int v1) {
  Edge temp_edge;
  temp_edge[0] = v0;
  temp_edge[1] = v1;
  int ei = edges.size();
  edges.push_back(temp_edge);
  vertices[v0].in_bd_of.push_back(ei);
  vertices[v1].in_bd_of.push_back(-ei);
  return ei;
}
  
//the edges read around counterclockwise
int Triangulation::add_triangle(SignedInd e0, SignedInd e1, SignedInd e2) {
  Triangle temp_tri;
  temp_tri[0] = e0; temp_tri[1] = e1; temp_tri[2] = e2;
  int ti = triangles.size();
  triangles.push_back(temp_tri);
  if (e0>0) edges[e0].in_bd_pos.push_back(ti); else edges[-e0].in_bd_neg.push_back(ti);
  if (e1>0) edges[e1].in_bd_pos.push_back(ti); else edges[-e1].in_bd_neg.push_back(ti);
  if (e2>0) edges[e2].in_bd_pos.push_back(ti); else edges[-e2].in_bd_neg.push_back(ti);
  return ti;
}  
  

void Triangulation::read_file(std::string filename) {}

void Triangulation::write_file(std::string filename) {}


/*****************************************************************************
 * set the triangulation to be a closed surface of genus g
 * ***************************************************************************/
void Triangulation::set_closed_surface(int genus) {
  //zigzag triangulation
  //the middle edges zigzag up face, starting at vertex 0->2, then 2->(2-3), then (2-3)->(2-3)+4...
  //it's the quotient of this, of course
  
  //add a single vertex and clear the other lists
  vertices = std::vector<Vertex>(2);
  edges.resize(1);
  triangles.resize(1);
  fundamental_loops.resize(0); //these will list a_1, b_1, a_2, b_2, ..., where \prod_i [a_i,b_i] = 1
  //add all the generator edges, which will be edges 1 through 2g
  for (int i=0; i<genus; ++i) {
    for (int j=0; j<2; ++j) {
      int ei = add_edge(1,1);
      fundamental_loops.push_back( std::vector<int>(1,ei) );
    }
  }
  //now add the edges that go up the middle, and add the triangles
  //make the first, special case triangle
  int just_added_edge = add_edge(1,1);
  (void)add_triangle(1,2,-just_added_edge); //this is the first triangle
  
  //set up the iterative step
  int previous_tile_vertex = 1;
  int current_tile_vertex = 0;
  int next_tile_vertex=2;
  int prev_added_edge=-1;
  int outside_edge=-1;
  int tile_step = 2;
  while (true) {
    //update the vertices
    previous_tile_vertex = current_tile_vertex;
    current_tile_vertex = next_tile_vertex;
    
    //compute the next tile vertex
    ++tile_step;
    if (tile_step == 4*genus-1) break;
    next_tile_vertex = current_tile_vertex + ((tile_step%2)==0 ? tile_step : -tile_step);
    next_tile_vertex = pos_mod(next_tile_vertex, 4*genus);
    
    //update the edges;
    prev_added_edge = just_added_edge;
    just_added_edge = add_edge(1,1);
    if ((tile_step%2)==0) { 
      //we're moving right, so the counterclockwise tile vertex before the outside edge is previous_tile_vertex
      outside_edge = free_gen_from_tile_index(previous_tile_vertex);
      //add the triangle
      (void)add_triangle( outside_edge, -just_added_edge, -prev_added_edge );
      
    } else {
      //we're moving left, so the counterclockwise tile vertex before the outside edge is next_tile_vertex
      outside_edge = free_gen_from_tile_index(next_tile_vertex);
      //add the triangle
      (void)add_triangle( prev_added_edge, just_added_edge, outside_edge );
    }
  }
  //add the last triangle
  prev_added_edge = just_added_edge;
  previous_tile_vertex = current_tile_vertex;
  current_tile_vertex = next_tile_vertex;
  next_tile_vertex = current_tile_vertex + 1;
  int outside_edge_1 = free_gen_from_tile_index(current_tile_vertex);
  int outside_edge_2 = free_gen_from_tile_index(next_tile_vertex);
  (void)add_triangle( prev_added_edge, outside_edge_1, outside_edge_2 );
}


/*****************************************************************************
 * Print a triangulation
 * ***************************************************************************/
void Triangulation::print(std::ostream& os) {
  os << "Vertices (" << vertices.size()-1 << "):\n";
  for (int i=1; i<(int)vertices.size(); ++i) {
    os << i << ": " << vertices[i] << "\n";
  }
  os << "Edges (" << edges.size()-1 << "):\n";
  for (int i=1; i<(int)edges.size(); ++i) {
    os << i << ": " << edges[i] << "\n";
  }
  os << "Triangles (" << triangles.size()-1 << "):\n";
  for (int i=1; i<(int)triangles.size(); ++i) {
    os << i << ": " << triangles[i] << "\n";
  }
  os << "Fundamental loops:\n";
  for (int i=0; i<(int)fundamental_loops.size(); ++i) {
    os << i << ": ";
    for (int j=0; j<(int)fundamental_loops[i].size(); ++j) {
      os << fundamental_loops[i][j] << " ";
    }
    os << "\n";
  }
}

/*****************************************************************************
 * compute Euler characteristic.  This requires that every edge have at most 
 * one positive incident triangle and at most one negative incident triangle.
 * If resolve_vertices=false, it will assume that there are no branch points
 * and just compute v-e+f.  If resolve_vertices=true, it'll figure out how many 
 * sheets of the surface touch together at the vertex, and it'll use the sum of 
 * these counts instead of v
 * ***************************************************************************/
int Triangulation::chi(bool resolve_vertices) {
  if (!resolve_vertices) {
    return (vertices.size()-1) - (edges.size()-1) + (triangles.size()-1);
  }
  int real_v=0;
  for (int i=1; i<(int)vertices.size(); ++i) {
    int num_ie = vertices[i].in_bd_of.size();
    std::vector<bool> have_done_incident_edge(num_ie, false);
    //find an undone edge
    //int first_ei = -1;
    for (int j=0; j<num_ie; ++j) {
      if (have_done_incident_edge[j]==false) {
        
      //NOT DONE
      }
    }
  }
  return real_v - (edges.size()-1) + (triangles.size()-1);
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
  for (int i=1; i<(int)weights.size(); ++i) {
    if (weights[i] == 0) continue;
    for (int j=0; j<3; ++j) {
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
  new_T.vertices.resize(num_vertices_used+1);
  new_T.edges.resize(num_edges_used+1);
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
    Triangle temp_tri;
    for (int j=0; j<3; ++j) {
      //copy the edges, except translate to the new edge indices
      temp_tri[j] = sgn(triangles[i][j])*edge_index_translation_table[ abs(triangles[i][j]) ];
    }
    for (int j=0; j<weights[i]; ++j) {
      int this_tri_ind = new_T.triangles.size();
      for (int k=0; k<3; ++k) {
        //record in the edge that this new triangle has it as boundary
        if (sgn(temp_tri[k]) > 0) {
          new_T.edges[abs(temp_tri[k])].in_bd_pos.push_back( this_tri_ind );
        } else {
          new_T.edges[abs(temp_tri[k])].in_bd_neg.push_back( this_tri_ind );
        }
      }
      new_T.triangles.push_back(temp_tri);
    }
  }
  
  return new_T;
}



Triangulation Triangulation::resolve_branched_surface() {
  return Triangulation();
}

int main(int argc, char* argv[]) {

  if (argc < 2) {
  
    Surface S(2,0);
    S.print(std::cout);
    std::string w("aabABbBcdCDA");
    std::string wg = S.geodesic(w);
    std::cout << "Word: " << w << " geodesic: " << wg << "\n";
    w = std::string("abABcababab");
    wg = S.geodesic(w);
    std::cout << "Word: " << w << " geodesic: " << wg << "\n";
    w = std::string("dababbaBA");
    wg = S.geodesic(w);
    std::cout << "Word: " << w << " geodesic: " << wg << "\n\n";
  
    Surface S2(3,2);
    S2.print(std::cout);
    w = std::string("abABcdCD");
    wg = S2.geodesic(w);
    std::cout << "Word: " << w << " geodesic: " << wg << "\n";
    w = std::string("ddcDCbaBA");
    wg = S2.geodesic(w);
    std::cout << "Word: " << w << " geodesic: " << wg << "\n\n";
  
    Surface S3(3,0);
    std::vector<std::string> W_words(3);
    W_words[0] = std::string("efEFabbca");
    W_words[1] = std::string("ababaccdcdeaea");
    W_words[2] = std::string("aBeDfcbaEd");
    LoopArrangement LA(S3, W_words, 2);
    std::cout << "After just initializing, " << LA.count_crossings() << " crossings:\n";
    LA.print(std::cout);
    LA.show();
  
    LA.minimal_position();
    std::cout << "After minimizing, " << LA.count_crossings() << " crossings:\n";
    LA.print(std::cout);
  
    LA.show();
    return 0;
  
  } else {
    std::vector<std::string> words(0);
    int genus, nboundaries;
    int current_arg = 1;
    int verbose=1;
    if (argc < 4) {
      std::cout << "usage: ./branched -v[n] <genus> <nbounaries> <loops>\n";
      return 0;
    }
    while (argv[current_arg][0] == '-') {
      switch(argv[current_arg][1]) {
        case 'v':
          if (argv[current_arg][2] == '\0') {
            verbose = 2;
          } else {
            verbose = atoi(&argv[current_arg][2]);
          }
          break;
      }
      ++current_arg;
    }     
    genus = atoi(argv[current_arg]);
    nboundaries = atoi(argv[current_arg+1]);
    for (int i=current_arg+2; i<argc; ++i) {
      words.push_back(std::string(argv[i]));
    }
    Surface S(genus, nboundaries, verbose);
    LoopArrangement LA(S, words, verbose);
    std::cout << "After just initializing, " << LA.count_crossings() << " crossings:\n";
    LA.print(std::cout);
    LA.show();
  
    LA.minimal_position();
    std::cout << "After minimizing, " << LA.count_crossings() << " crossings:\n";
    LA.print(std::cout);
  
    LA.show();
    return 0;
  }
}





