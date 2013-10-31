#ifndef __SURFACE_H__
#define __SURFACE_H__

#include <vector>
#include <string>
#include <iostream>
#include <map>

#include "branched.h"

struct LoopArrangement;

struct Surface {
  int genus;
  int nboundaries;
  int ngens;
  std::vector<int> cyclic_order;
  std::map<int, int> cyclic_order_map;
  std::vector<int> relator;
  std::vector<int> relator_inverse;
  std::map<int, int> relator_map;
  std::map<int, int> relator_inverse_map;
  
  //use the default string aBAbcDCd...fFgGhH...
  //The weird order is so that the relator is [a,b][c,d]...fgh...
  Surface(int g, int nb); 
  
  //print out all the data
  void print(std::ostream& os);
  
  //apply a relator at a position in a word (cyclically)
  //(replace xyz with xy'z where y(y'^-1) = relator, and y has length len
  //if inverse==true, it'll use the relator inverse instead
  //this ASSUMES that the subword at pos of length len is a subword of relator/inverse !!!
  void apply_relator(std::vector<int>& w, int pos, int len, bool inv=false);
  
  //reduce a word into a combinatorial geodesic (a string)
  std::string geodesic(std::string& w);
  
  //reduce a word into a combinatorial geodesic
  //this is NOT unique obviously
  std::vector<int> geodesic(std::vector<int>& w);
  
  //returns +/-1 depending on whether the gen list x,y,z is positively cyclically ordered
  int cyclically_ordered(SignedInd x, SignedInd y, SignedInd z);
  
  //returns +/-1 depending on how the words diverge
  //the starting letters must be the same!
  //+1 if the cyclic order is (start, w1, w2)
  int cyclically_ordered(std::vector<int>& w1, int start1, int dir1,
                         std::vector<int>& w2, int start2, int dir2);
  
  //finds the locations along the boundary edges associated to every letter 
  //in the input words such that the intersection number is minimized
  //this makes sure the input is combinatorially geodesic
  //this may change the input loops! (to something the same in the group)
  LoopArrangement minimal_intersection_position(std::vector<std::string>& W_words);
  
  //this draws a loop arrangement of the surface
  void draw_loop_arrangement(LoopArrangement& LA);
  
  //produces a branched surface (cellulation) corresponding to 
  //the complementary regions of the loops
  //The cellulation is topologically the surface, and the 
  //marked loops are part of the returned cellulation
  void cellulation_from_loops(std::vector<std::string> W);
};

//this helps with arranging the geodesics
struct GenPosition {
  std::vector<int>* word; //pointer to the actual word
  int w;                  //word index
  int i;                  //letter index
  Surface* S;             //the surface we're working in
};

//this allows us to package a loop arrangement
struct LoopArrangement {
  Surface* S;
  std::vector<std::string> W_words;
  std::vector<std::vector<int> > W;
  std::vector<std::vector<int> > positions;
  std::vector<int> gen_counts;
};

#endif