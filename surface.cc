
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
  cyclic_order = std::string("");
  cyclic_order_map.clear();
  relator = std::string("");
  for (int i=0; i<genus; ++i {
    cyclic_order = cyclic_order + std::string(alpha_ind_to_letter( 2*i+1 ))
                                + std::string(alpha_ind_to_letter( -((2*i+1)+1) ))
                                + std::string(alpha_ind_to_letter( -(2*i+1) ))
                                + std::string(alpha_ind_to_letter( (2*i+1)+1 ));
    cyclic_order_map[2*i+1] = 4*i;
    cyclic_order_map[-((2*i+1)+1)] = 4*i+1;
    cyclic_order_map[-(2*i+1)] = 4*i+2;
    cyclic_order_map[(2*i+1)+1] = 4*i+3;
    relator = relator + std::string(alpha_ind_to_letter( 2*i+1 ))
                      + std::string(alpha_ind_to_letter( -(2*i+1) ))
                      + std::string(alpha_ind_to_letter( (2*i+1)+1 ))
                      + std::string(alpha_ind_to_letter( -((2*i+1)+1) ));
  }
  for (int i=1; i<=nboundaries; ++i) {
    cyclic_order = cyclic_order + std::string(alpha_ind_to_letter(2*genus+i))
                                + std::string(alpha_ind_to_letter(-(2*genus+i)));
    cyclic_order_map[2*genus+i] = 4*genus + 2*i;
    cyclic_order_map[-(2*genus+i)] = 4*genus + 2*i + 1;
    relator = relator + std::string(alpha_ind_to_letter(2*genus+i));
  }
}
  
  
void Surface::print(std::ostream& os) {
  os << "Surface of genus " << genus << " with " << nboundaries << " boundary components\n";
  os << "Cyclic order action on polygon: " << cyclic_order << "\n";
  os << "Relator: " << relator << "\n";
  os << "Gens: ";
  for (int i=1; i<=ngens; ++i) {
    os << alpha_ind_to_letter(i) << " ";
  }
  os << "Cyclic order map: ";
  for (int i=1; i<=ngens; ++i) {
    os << alpha_ind_to_letter(i) << ": " << cyclic_order_map(i) << ", ";
    os << alpha_ind_to_letter(-i) << ": " << cyclic_order_map(-i) << ", ";
  }
  os << "\n";
}
  
  