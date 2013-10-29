#ifndef __SURFACE_H__
#define __SURFACE_H__

#include <string>

struct Surface {
  int genus;
  int nboundaries;
  std::string cyclic_order;
  
  //Use a cyclic order string; the letter gives the position of the 
  //*destination* of the action of that letter.
  Surface(std::string cyclic_order);
  
  //use the default string aBAbcDCd...fFgGhH...
  //The weird order is so that the boundary is [a,b][c,d]... + 
  Surface(int g, int nb); 
  

#endif