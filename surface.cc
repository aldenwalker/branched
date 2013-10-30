#include <sstream>

#include "surface.h"


char alpha_ind_to_letter(SignedInd ind) {
  if (ind > 0) {
    return (char)((ind-1)+97);
  } else {
    return (char)(((-ind)-1)+65);
  }
}

SignedInd letter_to_alpha_ind(char let) {
  int sgn = (let >= 97 ? 1 : -1);
  int val = (let >= 97 ? let-97 : let-65);
  return sgn*(val+1);
}

Surface::Surface(int g, int nb) {
  genus = g;
  nboundaries = nb;
  ngens = 2*g + nb;
  std::stringstream CO;
  cyclic_order_map.clear();
  std::stringstream R;
  relator_map.clear();
  relator_inverse_map.clear();
  //go around the normal closed surface part
  for (int i=0; i<genus; ++i) {
    CO << alpha_ind_to_letter( 2*i+1 )
       << alpha_ind_to_letter( -((2*i+1)+1) )
       << alpha_ind_to_letter( -(2*i+1) )
       << alpha_ind_to_letter( (2*i+1)+1 );
    cyclic_order_map[2*i+1] = 4*i;
    cyclic_order_map[-((2*i+1)+1)] = 4*i+1;
    cyclic_order_map[-(2*i+1)] = 4*i+2;
    cyclic_order_map[(2*i+1)+1] = 4*i+3;
    R << alpha_ind_to_letter( 2*i+1 )
      << alpha_ind_to_letter( (2*i+1)+1 )
      << alpha_ind_to_letter( -(2*i+1) )
      << alpha_ind_to_letter( -((2*i+1)+1) );
    relator_map[2*i+1] = 4*i;
    relator_map[(2*i+1)+1] = 4*i+1;
    relator_map[-(2*i+1)] = 4*i+2;
    relator_map[-((2*i+1)+1)] = 4*i+3;
  }
  //go over the boundary components
  for (int i=0; i<nboundaries; ++i) {
    CO << alpha_ind_to_letter(2*genus+i+1) << alpha_ind_to_letter(-(2*genus+i+1));
    cyclic_order_map[2*genus+i+1] = 4*genus + 2*i;
    cyclic_order_map[-(2*genus+i+1)] = 4*genus + 2*i + 1;
    R << alpha_ind_to_letter(2*genus+i+1);
    relator_map[2*genus+i+1] = 4*genus + i;
  }
  //create the relator inverse map
  for (int i=1; i<=2*genus; ++i) {
    relator_inverse_map[i] = (4*genus+nboundaries)-1-relator_map[-i];
    relator_inverse_map[-i] = (4*genus+nboundaries)-1-relator_map[i];
  }
  for (int i=1; i<=nboundaries; ++i) {
    relator_inverse_map[-(2*genus+i)] = (4*genus+nboundaries)-1-relator_map[(2*genus+i)];
  }
  relator = R.str();
  cyclic_order = CO.str();
  
}
  
  
void Surface::print(std::ostream& os) {
  os << "Surface of genus " << genus << " with " << nboundaries << " boundary components\n";
  os << "Cyclic order action on polygon: " << cyclic_order << "\n";
  os << "Relator: " << relator << "\n";
  os << "Gens: ";
  for (int i=1; i<=ngens; ++i) {
    os << alpha_ind_to_letter(i) << " ";
  }
  os << "\n";
  os << "Cyclic order map: ";
  for (std::map<int,int>::iterator it=cyclic_order_map.begin(); 
       it != cyclic_order_map.end(); 
       it++) {
    os << alpha_ind_to_letter(it->first) << ":" << it->second << ", ";
  }
  os << "\n";
  os << "Relator map: ";
  for (std::map<int,int>::iterator it=relator_map.begin(); 
       it != relator_map.end(); 
       it++) {
    os << alpha_ind_to_letter(it->first) << ":" << it->second << ", ";
  }
  os << "\n";
  os << "Relator inverse map: ";
  for (std::map<int,int>::iterator it=relator_inverse_map.begin(); 
       it != relator_inverse_map.end(); 
       it++) {
    os << alpha_ind_to_letter(it->first) << ":" << it->second << ", ";
  }
  os << "\n";
}
  
  