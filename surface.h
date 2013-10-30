#ifndef __SURFACE_H__
#define __SURFACE_H__

#include <vector>
#include <string>
#include <iostream>
#include <map>

#include "branched.h"

struct Surface {
  int genus;
  int nboundaries;
  int ngens;
  std::string cyclic_order;
  std::map<int, int> cyclic_order_map;
  std::string relator;
  std::string relator_inverse;
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
  void apply_relator(std::string& w, int pos, int len, bool inv=false);
  
  //reduce a word into a combinatorial geodesic
  //obviously, this is NOT unique if it's a closed surface
  std::string geodesic(std::string w);
  
  //returns +/-1 depending on whether the gen list x,y,z is positively cyclically ordered
  int cyclic_order(SignedInd x, SignedInd y, SignedInd z);
  
  //finds the locations along the boundary edges associated to every letter 
  //in the input words such that the intersection number is minimized
  //this makes sure the input is combinatorially geodesic
  //this may change the input loops! (to something the same in the group)
  void minimal_intersection_position(std::vector<std::string>& W, 
                                     std::vector<std::vector<int> >& positions,
                                     std::vector<int>& gen_counts);
  
  //produces a branched surface (cellulation) corresponding to 
  //the complementary regions of the loops
  //The cellulation is topologically the surface, and the 
  //marked loops are part of the returned cellulation
  void cellulation_from_loops(std::vector<std::string> W);
};

//this helps with arranging the geodesics
struct GenPosition {
  int w;    //word index
  int i;    //letter index
  int sign; //+1 means positive (lower case)
  std::string* word; //pointer to the actual word
  Surface* S; //the surface we're working in
};


#endif