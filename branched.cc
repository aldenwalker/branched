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

std::ostream& operator<<(std::ostream& os, const Edge& e) {
  os << "E(" << e.start << "," << e.end << "); {{";
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
  if (e.boundary_loop) os << " (boundary loop)";
  return os;
}

/*****************************************************************************
 * Cell
 * ***************************************************************************/
std::ostream& operator<<(std::ostream& os, Cell& c) {
  os << "C" << c.bd;
  if (c.contains_boundary) os << " (contains boundary)";
  if (c.computed_winding_number) os << " Winding number: " << c.winding_number;
  return os;
}

/*****************************************************************************
 * print a cellulation
 * ***************************************************************************/
void Cellulation::print(std::ostream& os) {
  os << "Vertices (" << vertices.size()-1 << "):\n";
  for (int i=1; i<(int)vertices.size(); ++i) {
    os << i << ": " << vertices[i] << "\n";
  }
  os << "Edges (" << edges.size()-1 << "):\n";
  for (int i=1; i<(int)edges.size(); ++i) {
    os << i << ": " << edges[i] << "\n";
  }
  os << "Cells (" << cells.size()-1 << "):\n";
  for (int i=1; i<(int)cells.size(); ++i) {
    os << i << ": " << cells[i] << "\n";
  }
  os << "Loops:\n";
  for (int i=0; i<(int)loops.size(); ++i) {
    os << i << ": " << loops[i] << "\n";
  }
}


/*****************************************************************************
 * both of the following functions follows as though the edges are on the 
 * boundary of a relator disk, and it's reading off the relator disk
 * this will give the INVERSE boundary as reading it off as a fatgraph
 *****************************************************************************/
/*****************************************************************************
 * get the next edge in the boundary:
 * this assumes that all vertices have cyclically ordered edges
 *****************************************************************************/
SignedInd Cellulation::next_edge(SignedInd e) {
  int v = (e>0 ? edges[abs(e)].end : edges[abs(e)].start);
  int vi = -1;
  for (int i=0; i<(int)vertices[v].in_bd_of.size(); ++i) {
    if (vertices[v].in_bd_of[i] == -e) {
      vi = i;
      break;
    }
  }
  int vi_prev = pos_mod(vi-1, vertices[v].in_bd_of.size());
  return vertices[v].in_bd_of[vi_prev];
}

/*****************************************************************************
 * follow an edge until we loop back
 * this assumes that all vertices have cyclically ordered edges
 *****************************************************************************/
std::vector<SignedInd> Cellulation::follow_edge(SignedInd e) {
  std::vector<SignedInd> ans(0);
  ans.push_back(e);
  int next_e = next_edge(e);
  while (next_e != e) {
    ans.push_back(next_e);
    next_e = next_edge(next_e);
  }
  return ans; 
}

/*****************************************************************************
 * draw the edges and the cells to the given XGraphics
 ****************************************************************************/
void Cellulation::draw_to_xgraphics(XGraphics& X) {

  bool label_edge_arrows = false;
  bool gray_from_winding = false;
  int wn_max = 0;
  if (cells[1].computed_winding_number) {
    gray_from_winding = true;
    for (int i=1; i<(int)cells.size(); ++i) {
      if (abs(cells[i].winding_number) > wn_max) {
        wn_max = abs(cells[i].winding_number);
      }
    }
  }    
  
  //draw the cells; make them random gray levels
  srand(2);
  for (int i=1; i<(int)cells.size(); ++i) {
    if (cells[i].sign < 0 || cells[i].contains_boundary) continue;
    double gray_level;
    if (gray_from_winding) {
      if (cells[i].winding_number == 0) {
        gray_level = 1;
      } else {
        gray_level = 0.5 + 0.45*(1-(double)abs(cells[i].winding_number)/(double)wn_max);
      }
    } else {
      double rand_gray_level = (double)rand()/(double)RAND_MAX;
      gray_level = rand_gray_level*0.3 + 0.5;
    }
      
    int col = X.get_rgb_color(gray_level, gray_level, gray_level);
    //std::cout << "Random gray level: " << rand_gray_level << "\n";
    //std::cout << "Returned color: " << col << "\n";
    std::vector<Point2d<float> > points(0);
    for (int j=0; j<(int)cells[i].bd.size(); ++j) {
      int e = cells[i].bd[j];
      if (e>0) {
        points.push_back( Point2d<float>(edges[e].start_pos.x.get_d(), 
                                         edges[e].start_pos.y.get_d()) );
      } else {
        if (edges[-e].two_sided) {
          points.push_back( Point2d<float>(edges[-e].end_neg.x.get_d(), 
                                           edges[-e].end_neg.y.get_d()) );
        } else {
          points.push_back( Point2d<float>(edges[-e].end_pos.x.get_d(), 
                                           edges[-e].end_pos.y.get_d()) );
        }
      }
    }
    X.draw_filled_polygon(points, col);
  }

  //draw the edges
  Point2d<float> temp1;
  Point2d<float> temp2;
  int black_color = X.get_color(std::string("black"));
  for (int i=1; i<(int)edges.size(); ++i) {
    temp1 = Point2d<float>(edges[i].start_pos.x.get_d(), edges[i].start_pos.y.get_d());
    temp2 = Point2d<float>(edges[i].end_pos.x.get_d(), edges[i].end_pos.y.get_d());
    std::stringstream edge_label_s;
    edge_label_s << i;
    std::string edge_label = (label_edge_arrows ? edge_label_s.str() : "");
    //X.draw_arrowed_labeled_line(temp1, temp2, black_color, 1, std::string(""));//edge_label);
    //X.draw_arrowed_labeled_line(temp1, temp2, black_color, 1, edge_label);
    if (edges[i].two_sided) {
      temp1 = Point2d<float>(edges[i].start_neg.x.get_d(), edges[i].start_neg.y.get_d());
      temp2 = Point2d<float>(edges[i].end_neg.x.get_d(), edges[i].end_neg.y.get_d());
      //X.draw_arrowed_labeled_line(temp1, temp2, black_color, 1, std::string(""));//edge_label);
      //X.draw_arrowed_labeled_line(temp1, temp2, black_color, 1, edge_label);
    }
  }

  //if we have a winding number, might as well print it
  for (int i=1; i<(int)cells.size(); ++i) {
    if (!cells[i].computed_winding_number) continue;
    std::stringstream edge_label_s;
    edge_label_s << cells[i].winding_number;
    std::string edge_label = edge_label_s.str();
    temp1 = Point2d<float>(cells[i].coords.x.get_d(), cells[i].coords.y.get_d());
    X.draw_text_centered(temp1, edge_label, black_color);
  }
  
}

/*****************************************************************************
 * compute the winding number of the cells, relative to one of the cells
 * this only works for cellulations coming from loop arrangements
 * ***************************************************************************/
void Cellulation::compute_winding_numbers(LoopArrangement& LA, int relative_to_cell) {
  if (relative_to_cell == 0) {
    //if the surface is closed, choose cell 1
    //if it has boundary, choose one of the cells on the boundary
    if (LA.S->nboundaries == 0) {
      relative_to_cell = 1;
    } else {
      for (int i=1; i<(int)cells.size(); ++i) {
        if (cells[i].contains_boundary && cells[i].bd.size() > 1) {
          relative_to_cell = i;
          break;
        }
      }
    }
  }
  if (verbose > 1) std::cout << "Computing winding number relative to " << relative_to_cell << "\n";
  for (int i=1; i<(int)cells.size(); ++i) {
    cells[i].winding_number = LA.algebraic_segment_intersection_number(cells[relative_to_cell].coords,
                                                                       cells[i].coords);
    cells[i].computed_winding_number = true;
  }
}
  
/*****************************************************************************
 * compute a naive (upper) bound for Euler characteristic surface with the desired boundary
 * in this case, every cell appears exactly as many times as its winding number
 * there's a choice about where to base the winding number computation 
 * 
 * In a surface with boundary, it'll always be based at the surface boundary
 * 
 * In a surface without boundary, we need to maximize over a 
 * finite range of possible values (by adding and subtracting the fundamental
 * class of the surface)
 * ***************************************************************************/
int Cellulation::chi_upper_bound(LoopArrangement& LA) {
  
  //compute winding numbers relative to an arbitrary cell
  compute_winding_numbers(LA);
  
  if (verbose > 1) std::cout << "Starting to compute chi upper bound\n";
  
  //if the surface has boundary, figure out what the right offset is 
  //to make the boundary cell have winding number 0
  //if it doesn't have boundary, find all possible values for the offset
  std::vector<int> offset_values(0);
  if (LA.S->nboundaries > 0) {
    //find the cell with boundary (and multiple sides)
    for (int i=1; i<(int)cells.size(); ++i) {
      if (cells[i].contains_boundary && cells[i].bd.size() > 0) {
        offset_values.push_back(-cells[i].winding_number);
        break;
      }
    }
  } else {
    int max_wn = 0; //we computed relative to a cell, so there is always a 0
    int min_wn = 0; //same
    for (int i=1; i<(int)cells.size(); ++i) {
      if (cells[i].winding_number > max_wn) max_wn = cells[i].winding_number;
      if (cells[i].winding_number < min_wn) min_wn = cells[i].winding_number;
    }
    if (max_wn > abs(min_wn)) { 
      //range is more positive than negative, so we can only subtract
      for (int i=0; abs(min_wn-i) <= max_wn; ++i) {
        offset_values.push_back(-i);
      }
    } else {
      //range is more negative than positive, so we can only add
      for (int i=0; max_wn+i <= abs(min_wn); ++i) {
        offset_values.push_back(i);
      }
    }
  }
  
  if (verbose > 1) std::cout << "Offset values: " << offset_values << "\n";
  
  //for each offset value, add it to the winding number for each 
  //cell, and total up the (naively computed) euler characteristic
  int largest_chi = 0;
  for (int off_i = 0; off_i<(int)offset_values.size(); ++off_i) {
    int off = offset_values[off_i];
    int cell_chi = 0;
    for (int i=1; i<(int)cells.size(); ++i) {
      if (cells[i].contains_boundary) continue;
      if (verbose > 2) std::cout << "cell " << i << " contributes " << abs(cells[i].winding_number + off) << "\n";
      cell_chi += abs(cells[i].winding_number + off);
    }
    if (verbose > 1) std::cout << "Total cells: " << cell_chi << "\n";
    int edge_chi = 0;
    for (int i=1; i<(int)edges.size(); ++i) {
      int edge_val = max( abs(cells[edges[i].in_bd_pos[0]].winding_number + off),
                          abs(cells[edges[i].in_bd_neg[0]].winding_number + off) );
      edge_chi += edge_val;
      if (verbose > 2) std::cout << "edge " << i << " contributes " << edge_val << "\n";
    }
    if (verbose > 1) std::cout << "Total edges: " << edge_chi << "\n";
    int vertex_chi = 0;
    for (int i=1; i<(int)vertices.size(); ++i) {
      int max_wn = 0;
      bool has_neg = false;
      bool has_pos = false;
      for (int j=0; j<(int)vertices[i].in_bd_of.size(); ++j) {
        int e = vertices[i].in_bd_of[j];
        int c = (e>0 ? edges[e].in_bd_neg[0] : edges[-e].in_bd_pos[0]);
        int wn = cells[c].winding_number + off;
        if (wn < 0) has_neg = true;
        if (wn > 0) has_pos = true;
        if (abs(wn) > max_wn) max_wn = abs(wn);
      }
      int vert_val;
      if (has_pos && has_neg) {
        vert_val = 2;
      } else {
        vert_val = max_wn;
      }
      if (verbose > 2) std::cout << "vertex " << i << " contributes " << vert_val << "\n";
      vertex_chi += vert_val;
    }
    if (verbose > 1) std::cout << "Total vertices: " << vertex_chi << "\n";
    int putative_chi = vertex_chi - edge_chi + cell_chi;
    if (verbose > 1) {
      std::cout << "Computed potential chi = " << putative_chi << " with offset " << off << "\n";
    }
    if (putative_chi > largest_chi || largest_chi == 0) {
      largest_chi = putative_chi;
    }
  }
  return largest_chi;
  
}




/*****************************************************************************
 * construct a branched surface
 * ***************************************************************************/
BranchedSurface::BranchedSurface(const Cellulation* const C, int verbose) {
  eperms_valid = false;
  this->C = C;
  this->verbose = verbose;
  cell_coefficients = std::vector<std::pair<int, int> >(C->cells.size(), std::make_pair(0,0));
  if (C->cells[1].computed_winding_number) { //use the winding numbers why not
    for (int i=1; i<(int)C->cells.size(); ++i) {
      int wn = C->cells[i].winding_number;
      if (wn >= 0) {
        cell_coefficients[i].first = wn;
      } else {
        cell_coefficients[i].second = -wn;
      }
    }
  }
}

BranchedSurface::BranchedSurface(const Cellulation* const C, 
                                 const std::vector<std::pair<int, int> >& cc,
                                 int verbose) {
  eperms_valid = false;
  this->C = C;
  cell_coefficients = cc;
  this->verbose = verbose;
}

int sum(std::pair<int, int>& p) {
  return p.first + p.second;
}

std::ostream& operator<<(std::ostream& os, std::pair<int, int>& p) {
  return os << "{" << p.first << "," << p.second << "}";
}

/*****************************************************************************
 * initialize the edge perms to a default
 * ***************************************************************************/
void BranchedSurface::init_edge_pdperms() {
  //set all the edge pdperms to be the uniform one (identify the low indices, 
  //and leave the top alone
  edge_pdperms.resize(C->edges.size());
  for (int i=1; i<(int)C->edges.size(); ++i) {
    int neg_cell = C->edges[i].in_bd_neg[0];
    int pos_cell = C->edges[i].in_bd_pos[0];
    if (verbose > 2) std::cout << "Edge " << i;
    if (verbose > 2) std::cout << " about to initialize a PDPerm( " << cell_coefficients[pos_cell].second << ","
                                                                   << cell_coefficients[neg_cell].first  << ","
                                                                   << cell_coefficients[neg_cell].second << ","
                                                                   << cell_coefficients[pos_cell].first << " )\n";
    edge_pdperms[i] = PDPerm(-cell_coefficients[pos_cell].second,  //cells on the right that contain me negatively
                             cell_coefficients[neg_cell].first,   //cells on the left containing positively
                             -cell_coefficients[neg_cell].second,  //cell on the right containing negatively
                             cell_coefficients[pos_cell].first);  //cells on the left containing positively
  }
}


/*****************************************************************************
 * compute the euler characteristic of the current gluing
 * ***************************************************************************/
int BranchedSurface::chi() {
  int ncells = 0;
  for (int i=1; i<(int)cell_coefficients.size(); ++i) {
    ncells += cell_coefficients[i].first + cell_coefficients[i].second;
  }
  if (verbose>1) std::cout << "I found there were " << ncells << " cells\n";
  int nedges = 0;
  for (int i=1; i<(int)C->edges.size(); ++i) {
    nedges += edge_pdperms[i].max_size();
  }
  if (verbose>1) std::cout << "I found there were " << nedges << " edges\n";
  int nverts = 0;
  for (int i=1; i<(int)C->vertices.size(); ++i) {
    nverts += num_vertices_over_vertex(i);
  }
  if (verbose>1) std::cout << "I found there were " << nverts << " vertices\n";
  return nverts - nedges + ncells;
}


/******************************************************************************
 * follow the glued edges around to figure out where the boundary is
 * direction = 1 means cross right to left in the cyclic order on the vertex
 * when we are "at" an edge, it means we're sitting immediately to the 
 * right in the cyclic order, so if we're going positively, 
 * we apply the edge pdperm as we move away
 * if we're going negatively, it means we apply the edge pdperm as 
 * we arrive at the edge
 *****************************************************************************/
void BranchedSurface::follow_gluing_around_vertex(const Vertex& vert, 
                                                  int start_edge, 
                                                  int start_level, 
                                                  int direction, 
                                                  std::vector<DSList<bool> >& edges_visited, 
                                                  int& boundary) {
  int nedges = vert.in_bd_of.size();
  int current_edge = start_edge;
  int current_level = start_level;
  int next_edge, next_level;
  int global_ei = vert.in_bd_of[current_edge];
  int next_global_ei;
  while (true) {
    if (direction > 0) {
      next_level = (global_ei > 0 
                         ? edge_pdperms[global_ei].map[current_level] 
                         : edge_pdperms[-global_ei].inverse_map[current_level] );
      edges_visited[current_edge][current_level] = true;
      next_edge = (next_level > 0 ? pos_mod(current_edge+1, nedges) 
                                  : pos_mod(current_edge-1, nedges) );
      next_global_ei = vert.in_bd_of[next_edge];
    } else {
      next_edge = (current_level > 0 ? pos_mod(current_edge-1, nedges) 
                                     : pos_mod(current_edge+1, nedges) );
      next_global_ei = vert.in_bd_of[next_edge];
      edges_visited[current_edge][current_level] = true;
      next_level = (next_global_ei > 0 
                      ? edge_pdperms[next_global_ei].inverse_map[current_level] 
                      : edge_pdperms[-next_global_ei].map[current_level] );
    }
    if (next_edge == start_edge && next_level == start_level) {
      boundary = -1;
      return;
    }
    if (next_level == 0) { //we've reached the boundary
      boundary = (direction > 0 ? current_edge : next_edge);
      return;
    }
    current_edge = next_edge;
    current_level = next_level;
    global_ei = next_global_ei;
  }
}

/*****************************************************************************
 * use the gluing to determine the number of vertices over the given vertex
 * ***************************************************************************/
int BranchedSurface::num_vertices_over_vertex(int vi, bool quiet) {
  const Vertex& vert = C->vertices[vi];
  int nedges = vert.in_bd_of.size();
  //this vector records whether we've gone over the edge
  //we go over an edge when we touch the right side (regardless of edge direction)
  std::vector<DSList<bool> > edges_visited(nedges);
  for (int i=0; i<nedges; ++i) {
    int ei = vert.in_bd_of[i];
    if (ei > 0) {
      edges_visited[i] = DSList<bool>( edge_pdperms[ei].smin(), edge_pdperms[ei].smax(), false );
    } else {
      edges_visited[i] = DSList<bool>( edge_pdperms[-ei].dmin(), edge_pdperms[-ei].dmax(), false);
    }
  }
  if (verbose > 2 && !quiet) std::cout << "Visiting vertex " << vi << "\n";
  int num_cycles = 0; //the number of cycles (including boundary cycles)
  int twice_num_extra_vertices_needed=0; //the number of boundary cycles which don't match up (and must be glued)
  while (true) {
    //find an edge that we haven't explored
    int start_edge=-1;
    int start_level=0;
    for (int i=0; i<nedges; ++i) {
      for (int j=edges_visited[i].min(); j<=edges_visited[i].max(); ++j) {
        if (j==0) continue;
        if (edges_visited[i][j] == false) {
          start_level = j; break;
        }
      }
      if (start_level != 0) {
        start_edge = i; break;
      }
    }
    if (start_edge == -1) break;
    
    int pos_dir_boundary; //this records the edge index that is the boundary
    follow_gluing_around_vertex(vert, start_edge, start_level, 1, edges_visited, pos_dir_boundary);
    
    if (pos_dir_boundary == -1) { //if there's no boundary, there's just a cycle
      ++num_cycles;
      continue;
    }
    
    int neg_dir_boundary;
    follow_gluing_around_vertex(vert, start_edge, start_level, -1, edges_visited, neg_dir_boundary);
    
    //this only works for a vertex of valence 4!
    if (vert.in_bd_of.size() != 4) std::cout << "Vertex has valence > 4?\n";
    if (pos_dir_boundary%2 != neg_dir_boundary%2) ++twice_num_extra_vertices_needed;
    ++num_cycles;
  }
  if (verbose > 2 && !quiet) std::cout << "I found " << num_cycles << " cycles and " << twice_num_extra_vertices_needed/2 << " extra vertices\n";
  if (twice_num_extra_vertices_needed%2 != 0) std::cout << "Extra vertex count weird?\n";
  
  return num_cycles - (twice_num_extra_vertices_needed/2);
}

/*****************************************************************************
 * compute chi, but only use vertices determined by the edges in edge_is_set
 * for the other vertices, it assumes the best case, so it produces 
 * an upper bound for the real chi
 *****************************************************************************/
int BranchedSurface::partially_defined_chi(const std::vector<bool>& edge_is_set) {
  int ncells = 0;
  for (int i=1; i<(int)cell_coefficients.size(); ++i) {
    ncells += cell_coefficients[i].first + cell_coefficients[i].second;
  }
  if (verbose>1) std::cout << "I found there were " << ncells << " cells\n";
  int nedges = 0;
  for (int i=1; i<(int)C->edges.size(); ++i) {
    nedges += edge_pdperms[i].max_size();
  }
  if (verbose>1) std::cout << "I found there were " << nedges << " edges\n";
  int nverts = 0;
  for (int i=1; i<(int)C->vertices.size(); ++i) {
    //determine if the vertex is computable
    const Vertex& vert = C->vertices[i];
    bool can_compute = true;
    int num_sheets = 0;
    for (int j=0; j<(int)vert.in_bd_of.size(); ++j) {
      int ei = abs(vert.in_bd_of[j]);
      if (!edge_is_set[ei]) {
        can_compute = false;
      }
      int sheets_over_edge = edge_pdperms[ei].max_size();
      num_sheets  = (num_sheets > sheets_over_edge ? num_sheets : sheets_over_edge);
    }
    if (can_compute) {
      nverts += num_vertices_over_vertex(i);
    } else {
      nverts += num_sheets;
    }
  }
  if (verbose>1) std::cout << "I found there were " << nverts << " vertices\n";
  return nverts - nedges + ncells;
  return 0;
}



/*****************************************************************************
 * optimize over all gluings for the given coefficients
 * this leaves the edge pdperms in the best gluing
 * ***************************************************************************/
int BranchedSurface::brute_minimal_gluing() {
  //first, we should get a good guess for what the gluings
  //should be, as this should let us trim the tree more
  int best_chi_found = guess_minimal_gluing();
  
  //save the best edge pdperms
  std::vector<PDPerm> best_gluing = edge_pdperms;
  
  //reset to the initial position
  init_edge_pdperms();
  
  //this vector records which edges we have decided
  std::vector<bool> edge_is_set(C->edges.size(), false);
  
  //the stack records which edges we have altered, in order
  std::vector<int> stack(0);
  
  //initialize the stack
  edge_is_set[1] = true;
  stack.push_back(1);
  int next_edge;
  
  if (verbose>1) std::cout << "Doing brute force gluing optimization\n";
  
  while (true) {
    if (verbose > 2) std::cout << "Current stack: " << stack << "\n";
    //compute the potential chi
    int chi_ub = partially_defined_chi(edge_is_set);
    if (verbose>2) std::cout << "Found chi upper bound of " << chi_ub << " compared to best chi of " << best_chi_found << "\n";
    if (chi_ub <= best_chi_found) goto BACKTRACK;
    //choose the next edge
    next_edge = 0;
    for (int i=1; i<(int)edge_is_set.size(); ++i) {
      if (edge_is_set[i] == false) {
        next_edge = i;
        break;
      }
    }
    if (verbose > 2) std::cout << "Found next edge " << next_edge << "\n";
    //if we can't choose the next edge, our chi calculation was a real one
    if (next_edge == 0) {
      if (chi_ub > best_chi_found) {
        if (verbose > 2) std::cout << "** new record\n";
        best_chi_found = chi_ub;
        best_gluing = edge_pdperms;
      }
      goto BACKTRACK;
    }
    //otherwise, set it up
    edge_pdperms[next_edge].reset(); //initialize it to the first pdperm
    edge_is_set[next_edge] = true;   //say we set it
    stack.push_back(next_edge);      //push it on
    continue;
    
    BACKTRACK:
    if (verbose > 2) std::cout << "Backtracking\n";
    bool completely_done = false;
    while (true) {
      if (stack.size() == 0) {
        completely_done = true; 
        break;
      }
      int current_ei = stack.back();
      if (verbose > 2) std::cout << "Finding next perm on edge " << current_ei << "\n";
      //advance to the next pdperm
      bool done_this_edge = edge_pdperms[current_ei].next();
      if (done_this_edge) {
        if (verbose > 2) std::cout << "Done this edge; backtracking more\n";
        edge_is_set[current_ei] = false;
        stack.pop_back();
      } else {
        break;
      }
    }
    if (completely_done) break;
  }
  edge_pdperms = best_gluing;
  return best_chi_found;
}

/*****************************************************************************
 * Try to get a good gluing by just hillclimbing
 *****************************************************************************/
int BranchedSurface::hillclimb_minimal_gluing() {
  return 0;
}

/*****************************************************************************
 * just guess a good gluing
 *****************************************************************************/
int BranchedSurface::guess_minimal_gluing(int seed) {
  init_edge_pdperms();
  return chi();
}


/******************************************************************************
 * print out a branched surface
 * ****************************************************************************/
void BranchedSurface::print(std::ostream& os) {
  os << "Branched surface on cellulation with " << C->vertices.size()-1 << " vertices, " 
                                                << C->edges.size()-1 << " edges, and " 
                                                << C->cells.size()-1 << " cells\n";
  os << "Cell coefficients:\n";
  for (int i=1; i<(int)C->cells.size(); ++i) {
    os << "(" << i << ": " << cell_coefficients[i] << "), ";
  }
  os << "\n";
  os << "Edge PDPerms:\n";
  for (int i=1; i<(int)C->edges.size(); ++i) {
    os << i << ": (" << C->edges[i] << ") " << edge_pdperms[i] << "\n";
  }
}












  
  
int main(int argc, char* argv[]) {

  std::vector<std::string> words(0);
  int genus, nboundaries;
  int current_arg = 1;
  int verbose=1;
  bool demo=false;
  bool winding_numbers=false;
  if (argc < 2 || (argc < 4 && argv[1][1] != 'd')) {
    std::cout << "usage: ./branched -v[n] [-w] [-d] <genus> <nbounaries> <loops>\n";
    std::cout << "\t-v[n]: verbose (of level n; default=2)\n";
    std::cout << "\t   -w: show complementary regions and winding numbers\n"; 
    std::cout << "\t   -d: demo (genus 3 example)\n";
    return 0;
  }
  while (current_arg < argc && argv[current_arg][0] == '-') {
    switch(argv[current_arg][1]) {
      case 'v':
        if (argv[current_arg][2] == '\0') {
          verbose = 2;
        } else {
          verbose = atoi(&argv[current_arg][2]);
        }
        break;
      case 'w':
        winding_numbers=true;
        break;
      case 'd':
        demo = true;
        break;
    }
    ++current_arg;
  }
  if (demo) {
    words.resize(3);
    words[0] = std::string("cdfABFaeCFc");
    words[1] = std::string("aDfAEBCDDedBE");
    words[2] = std::string("bbddbD");
    genus = 3;
    nboundaries = 0;
  } else {
    genus = atoi(argv[current_arg]);
    nboundaries = atoi(argv[current_arg+1]);
    for (int i=current_arg+2; i<argc; ++i) {
      words.push_back(std::string(argv[i]));
    }
  }
  Surface S(genus, nboundaries, verbose);
  if (verbose>1) {
    S.print(std::cout);
  }
  LoopArrangement LA(S, words, verbose);
  if (verbose>1) {
    std::cout << "After just initializing, " << LA.count_crossings() << " crossings:\n";
    LA.print(std::cout);
    LA.show();
  }

  LA.minimal_position();
  LA.find_crossing_data();
  std::cout << "Minimal position has " << LA.count_crossings() << " crossings.\n";
  if (verbose>1) {
    LA.print(std::cout);
  }
  LA.show();
  
  if (!winding_numbers) return 0;
  
  Cellulation C = LA.cellulation_from_loops();
  C.verbose = verbose;
  
  //compute the winding numbers
  C.compute_winding_numbers(LA);
  
  std::cout << "Displaying complementary regions and winding numbers.\n";
  if (verbose>1) { 
    C.print(std::cout);
  }
  LA.show(&C);
  
  return 0;
  
  int chi_ub = C.chi_upper_bound(LA);
  Rational scl_lb(-chi_ub, 2);
  std::cout << "Naive upper bound on chi: " << chi_ub << "\n";
  std::cout << "Implies lower bound on scl: " << scl_lb << "\n";
  LA.show(&C);
  
  BranchedSurface BS(&C, verbose);
  std::cout << "Initialized branched surface\n";
  BS.init_edge_pdperms();
  std::cout << "Initialized edge perms\n";
  BS.print(std::cout);
  
  std::cout << "Finding initial chi:\n";
  int chi = BS.chi();
  std::cout << "chi = " << chi << "\n";
  
  std::cout << "Finding a good guess chi:\n";
  chi = BS.guess_minimal_gluing();
  std::cout << "chi = " << chi << "\n";
  
  std::cout << "Finding chi via hillclimb\n";
  chi = BS.hillclimb_minimal_gluing();
  std::cout << "chi = " << chi << "\n";
  
  std::cout << "Finding chi with brute force\n";
  chi = BS.brute_minimal_gluing();
  std::cout << "chi = " << chi << "\n";
  
  std::cout << "Best gluing:\n";
  BS.print(std::cout);
  
  return 0;

}





