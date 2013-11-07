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

  //draw the cells; make them random gray levels
  srand(2);
  for (int i=1; i<(int)cells.size(); ++i) {
    if (cells[i].sign < 0 || cells[i].contains_boundary) continue;
    double rand_gray_level = (double)rand()/(double)RAND_MAX;
    rand_gray_level = rand_gray_level*0.3 + 0.5;
    int col = X.get_rgb_color(rand_gray_level, rand_gray_level, rand_gray_level);
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
    std::string edge_label = edge_label_s.str();
    X.draw_arrowed_labeled_line(temp1, temp2, black_color, 1, edge_label);
    if (edges[i].two_sided) {
      temp1 = Point2d<float>(edges[i].start_neg.x.get_d(), edges[i].start_neg.y.get_d());
      temp2 = Point2d<float>(edges[i].end_neg.x.get_d(), edges[i].end_neg.y.get_d());
      X.draw_arrowed_labeled_line(temp1, temp2, black_color, 1, edge_label);
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
    for (int i=0; i<(int)cells.size(); ++i) {
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
    if (verbose > 2) std::cout << "Total cells: " << cell_chi << "\n";
    int edge_chi = 0;
    for (int i=1; i<(int)edges.size(); ++i) {
      int edge_val = max( abs(cells[edges[i].in_bd_pos[0]].winding_number + off),
                          abs(cells[edges[i].in_bd_neg[0]].winding_number + off) );
      edge_chi += edge_val;
      if (verbose > 2) std::cout << "edge " << i << " contributes " << edge_val << "\n";
    }
    if (verbose > 2) std::cout << "Total edges: " << edge_chi << "\n";
    int vertex_chi = 0;
    for (int i=1; i<(int)vertices.size(); ++i) {
      int max_wn = 0;
      for (int j=0; j<(int)vertices[i].in_bd_of.size(); ++j) {
        int e = vertices[i].in_bd_of[j];
        int c = (e>0 ? edges[e].in_bd_neg[0] : edges[-e].in_bd_pos[0]);
        int wn = cells[c].winding_number + off;
        if (abs(wn) > max_wn) max_wn = abs(wn);
      }
      if (verbose > 2) std::cout << "vertex " << i << " contributes " << max_wn << "\n";
      vertex_chi += max_wn;
    }
    if (verbose > 2) std::cout << "Total vertices: " << vertex_chi << "\n";
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
  
  
  
  
  
  
  
  
int main(int argc, char* argv[]) {

  std::vector<std::string> words(0);
  int genus, nboundaries;
  int current_arg = 1;
  int verbose=1;
  bool demo=false;
  if (argc < 2 || (argc < 4 && argv[1][1] != 'd')) {
    std::cout << "usage: ./branched -v[n] <genus> <nbounaries> <loops>\n";
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
      case 'd':
        demo = true;
        break;
    }
    ++current_arg;
  }
  if (demo) {
    words.resize(3);
    words[0] = std::string("efEFabbca");
    words[1] = std::string("ababaccdcdeaea");
    words[2] = std::string("aBeDfcbaEd");
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
  S.print(std::cout);
  LoopArrangement LA(S, words, verbose);
  std::cout << "After just initializing, " << LA.count_crossings() << " crossings:\n";
  LA.print(std::cout);
  LA.show();

  LA.minimal_position();
  LA.find_crossing_data();
  std::cout << "After minimizing, " << LA.count_crossings() << " crossings:\n";
  LA.print(std::cout);
  LA.show();
  
  Cellulation C = LA.cellulation_from_loops();
  C.verbose = verbose;
  
  //compute the winding numbers
  C.compute_winding_numbers(LA);
  
  std::cout << "Cellulation from loops:\n";
  C.print(std::cout);
  
  int chi_ub = C.chi_upper_bound(LA);
  Rational scl_lb(-chi_ub, 2);
  std::cout << "Naive upper bound on chi: " << chi_ub << "\n";
  std::cout << "Implies lower bound on scl: " << scl_lb << "\n";
  
  LA.show(&C);
  
  return 0;

}





