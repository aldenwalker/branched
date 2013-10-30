#ifndef __SURFACE_H__
#define __SURFACE_H__

#include <vector>
#include <string>
#include <iostream>
#include <map>

struct Surface {
  int genus;
  int nboundaries;
  int ngens;
  std::string cyclic_order;
  std::map<int, int> cyclic_order_map;
  std::string relator;
  
  //use the default string aBAbcDCd...fFgGhH...
  //The weird order is so that the relator is [a,b][c,d]...fgh...
  Surface(int g, int nb); 
  
  //print out all the data
  void print(std::ostream& os);
  
  //reduce a word into a combinatorial geodesic
  //obviously, this is NOT unique if it's a closed surface
  std::string geodesic(std::string w);
  
  //finds the locations along the boundary edges associated to every letter 
  //in the input words such that the intersection number is minimized
  //this makes sure the input is combinatorially geodesic
  std::vector<std::vector<int> > minimal_intersection_position(std::vector<std::string> W);
  
  //produces a branched surface (cellulation) corresponding to 
  //the complementary regions of the loops
  //The cellulation is topologically the surface, and the 
  //marked loops are part of the returned cellulation
  void cellulation_from_loops(std::vector<std::string> W);
};

#endif