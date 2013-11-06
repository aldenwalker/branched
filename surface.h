#ifndef __SURFACE_H__
#define __SURFACE_H__

#include <vector>
#include <string>
#include <iostream>
#include <map>

#include "point.h"
#include "rational.h"
#include "branched.h"

struct LoopArrangement;

struct Crossing;

struct Segment;

struct GenPosition; 

std::ostream& operator<<(std::ostream& os, std::vector<int>& L);



/****************************************************************************
 * Surface class
 * mostly records the cyclic order of things, and geodesicifying
 ****************************************************************************/
struct Surface {
  int genus;
  int nboundaries;
  int ngens;
  int verbose;
  std::vector<int> cyclic_order;
  std::map<int, int> cyclic_order_map;
  std::vector<int> relator;
  std::vector<int> relator_inverse;
  std::map<int, int> relator_map;
  std::map<int, int> relator_inverse_map;
  
  //these record the points of the polygon
  //note these are close rational approximations of roots of unity
  std::map<int, Point2d<Rational> > gen_edge_start; 
  std::map<int, Point2d<Rational> > gen_edge_end;
  
  //use the default string aBAbcDCd...fFgGhH...
  //The weird order is so that the relator is [a,b][c,d]...fgh...
  Surface(int g, int nb, int verbose=1); 
  
  //check if a "generator" is actually a boundary loop placeholder
  bool is_boundary_loop(int gen);
  
  //check if a generator is a boundary generator
  bool is_boundary_gen(int gen);
  
  //print out all the data
  void print(std::ostream& os);
  
  //apply a relator at a position in a word (cyclically)
  //(replace xyz with xy'z where y(y'^-1) = relator, and y has length len
  //if inverse==true, it'll use the relator inverse instead
  //this ASSUMES that the subword at pos of length len is a subword 
  //of relator/inverse !!!
  void apply_relator(std::vector<int>& w, int pos, int len, bool inv=false);
  
  //reduce a word into a combinatorial geodesic (a string)
  std::string geodesic(std::string& w);
  
  //reduce a word into a combinatorial geodesic
  //this is unique if you ask for it
  void make_geodesic(std::vector<int>& w, bool unique_geodesics=false);
  
  //returns +/-1 depending on whether the gen list x,y,z is 
  //positively cyclically ordered
  int cyclically_ordered(SignedInd x, SignedInd y, SignedInd z);
  
  //returns +/-1 depending on how the words diverge
  //the starting letters must be the same!
  //+1 if the cyclic order is (start, w1, w2)
  int cyclically_ordered(std::vector<int>& w1, int start1, int dir1,
                         std::vector<int>& w2, int start2, int dir2);
  
};


//this helps with arranging the geodesics
struct GenPosition {
  Surface* S;             //the surface we're working in
  std::vector<int>* word; //pointer to the actual word
  int w;                  //word index
  int i;                  //letter index
};


//this is needed to use Point2d<Rational>'s as keys
struct dictionary_order {
  bool operator()(const Point2d<Rational>& a, const Point2d<Rational>& b) {
    if (a.x == b.x) {
      return (a.y < b.y);
    }
    return (a.x < b.x);
  }
};

/****************************************************************************
 * A LoopArrangement is a collection of homotopy classes
 * this is what handles finding the minimal position, etc
 ****************************************************************************/
struct LoopArrangement {
  Surface* S;
  std::vector<std::string> W_words;
  std::vector<std::vector<int> > W;
  
  //these record the indices of the letters
  std::vector<std::vector<int> > positions_by_letter;
  std::vector<std::vector<GenPosition> > positions_by_gen;
  
  //this is a list of Segments, which record the crossing data, etc
  //this list will be 1-indexed!
  std::vector<Segment> segments;
  std::vector<std::vector<int> > segments_by_letter;

  //this is a list of Crossings, which record which segments are incident, etc
  //this list will also be 1-indexed to be the same as the segment list
  std::map< Point2d<Rational>, int, dictionary_order> crossings_by_coords;
  std::vector<Crossing> crossings;

  //should be print messages? 1=normal 2=some 3=lots
  int verbose;
  
  LoopArrangement(Surface& S, int verbose=1);
  LoopArrangement(Surface& S, std::vector<std::string>& W_words, int verbose=1);
  LoopArrangement(Surface& S, std::vector<std::vector<int> >& W, int verbose=1);
  void init_from_vectors(Surface& S, 
                        std::vector<std::vector<int> >& W,
                        bool unique_geodesics=false);
  
  //given the positions along the gen edges, get the position for each letter
  void generate_positions_by_letter();
  
  //check the cyclic order of three positions
  int cyclically_ordered_positions(int gen1, int pos1, 
                                   int gen2, int pos2,
                                   int gen3, int pos3);
  
  //check if two strands cross at the given location
  bool check_cross(int w1, int i1, int w2, int i2);
  
  //just do a basic count of the crossings
  //note higher valence vertices will be counted multiple times
  int count_crossings();
  
  //find the positions of all the segments (not the crossing data)
  void find_segment_coordinates();
  
  void find_segment_crossing_coordinates(Segment& s1, Segment& s2, bool& do_cross,
                                         Rational& s1_t, Rational& s2_t,
                                         Point2d<Rational>& cross_coords);
  
  //compute the algebraic intersection of a line segment wilth all the segments
  int algebraic_segment_intersection_number(const Point2d<Rational>& p1,
                                            const Point2d<Rational>& p2);
  
  //find all the detailed crossing information
  void find_crossing_data();
  
  
  //finds the locations along the boundary edges associated to every letter 
  //in the input words such that the intersection number is minimized
  //this makes sure the input is combinatorially geodesic
  //this may change the input loops! (to something the same in the group)
  void minimal_position();
  
  //this draws a loop arrangement of the surface
  void show(Cellulation* C=NULL);
  
  //prints out the data
  void print(std::ostream& os);
  
  //produces a branched surface (cellulation) corresponding to 
  //the complementary regions of the loops
  //The cellulation is topologically the surface, and the 
  //marked loops are part of the returned cellulation
  Cellulation cellulation_from_loops();
  
};


//this is the data associated with a segment (single pair of letters)
//these record the actual rational coordinates of the start and endpoints
//note that it'll start inside the generator which is -W[w][i1]
struct Segment {
  Surface* S;
  LoopArrangement* LA;
  int w, i1, i2;  //this is the index of the first and second letters
  Point2d<Rational> start;
  Point2d<Rational> end;
  std::vector<std::pair<Rational, int> > crossings;
};

std::ostream& operator<<(std::ostream& os, Segment& s);

//this allows us to package a crossing
//note crossings happen between letters, so each crossing
//gives the index of the letter just *before* the crossing for each word
struct Crossing {
  Surface* S;
  LoopArrangement* LA;
  Point2d<Rational> coords;
  std::vector<int> segments; //the list of incident (signed) segments, in cyclic order
};

std::ostream& operator<<(std::ostream& os, Crossing& c);
  







#endif