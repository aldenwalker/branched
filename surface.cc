#include <sstream>
#include <algorithm>

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

char swap_case(char c) {
  if (std::isupper(c)) return std::tolower(c);
  else if (std::islower(c)) return std::toupper(c);
  else return 0;
}

std::string inverse(std::string& w) {
  std::string ans = w;
  std::reverse(ans.begin(), ans.end());
  std::transform(ans.begin(), ans.end(), ans.begin(), swap_case);
  return ans;
}

void reduce(std::string& w) {
  int i=0; 
  while (i < (int)w.size()-1) {
    if (w[i] == swap_case(w[i+1])) {
      w.replace(i,2, std::string(""));
      if (i>0) --i;
    } else {
      ++i;
    }
  }
}

void cyclically_reduce(std::string& w) {
  reduce(w);
  int wL = w.size();
  int i=0;
  while (i < wL/2 && w[i] == swap_case(w[wL-1-i])) {
    ++i;
  }
  w = w.erase(wL-1-(i-1), i);
  w = w.erase(0,i);
}
  
std::string cyclic_subword(std::string& w, int pos, int len) {
  int wL = w.size();
  if (pos+len > wL) {
    return w.substr(pos, wL-pos) + cyclic_subword(w, 0, len-(wL-pos));
  } else {
    return w.substr(pos, len);
  }
}

int cyclic_word_agreement_length(std::string w1, int pos1, std::string w2, int pos2) {
  int w1L = w1.size();
  int w2L = w2.size();
  int ans = 0;
  int w1pos = pos1;
  int w2pos = pos2;
  while (w1[w1pos%w1L] == w2[w2pos%w2L]) {
    ++ans;
    ++w1pos;
    ++w2pos;
    if (ans == w1L || ans == w2L) return ans;
  }
  return ans;
}


Surface::Surface(int g, int nb) {
  genus = g;
  nboundaries = nb;
  ngens = 2*g + nb;
  std::stringstream CO;
  std::stringstream R;
  //go around the normal closed surface part
  for (int i=0; i<genus; ++i) {
    CO << alpha_ind_to_letter( 2*i+1 )
       << alpha_ind_to_letter( -((2*i+1)+1) )
       << alpha_ind_to_letter( -(2*i+1) )
       << alpha_ind_to_letter( (2*i+1)+1 );
    R << alpha_ind_to_letter( 2*i+1 )
      << alpha_ind_to_letter( (2*i+1)+1 )
      << alpha_ind_to_letter( -(2*i+1) )
      << alpha_ind_to_letter( -((2*i+1)+1) );
  }
  //go over the boundary components
  for (int i=0; i<nboundaries; ++i) {
    CO << alpha_ind_to_letter(2*genus+i+1) << alpha_ind_to_letter(-(2*genus+i+1));
    R << alpha_ind_to_letter(2*genus+i+1);
  }
  relator = R.str();
  cyclic_order = CO.str();
  relator_inverse = inverse(relator);
  
  //create the maps for fast index access
  cyclic_order_map.clear();
  relator_map.clear();
  relator_inverse_map.clear();
  for (int i=0; i<(int)cyclic_order.size(); ++i) {
    cyclic_order_map[ letter_to_alpha_ind( cyclic_order[i])] = i;
  }
  for (int i=0; i<(int)relator.size(); ++i) {
    relator_map[ letter_to_alpha_ind( relator[i])] = i;
    relator_inverse_map[ letter_to_alpha_ind( relator_inverse[i])] = i;
  }
}
  
  
void Surface::print(std::ostream& os) {
  os << "Surface of genus " << genus << " with " << nboundaries << " boundary components\n";
  os << "Cyclic order action on polygon: " << cyclic_order << "\n";
  os << "Relator: " << relator << "\n";
  os << "Relator inverse: " << relator_inverse << "\n";
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


void Surface::apply_relator(std::string& w, int pos, int len, bool inv) {
  int wL = w.size();
  //get the position in the relator of the letter at pos+len
  int rel_pos = ( inv  
                  ? relator_inverse_map[letter_to_alpha_ind(w[(pos+len-1)%wL])]
                  : relator_map[letter_to_alpha_ind(w[(pos+len-1)%wL])] );
  rel_pos = (rel_pos+1)%relator.size();
  //std::cout << "Replacing subword of length " << len << " at pos " << pos << "\n";
  //std::cout << "Got relator position " << rel_pos << "\n";
  //get the cyclic subword of relator starting at rel_pos of 
  //length relator_length - len
  std::string new_relator_chunk = cyclic_subword( (inv ? relator_inverse : relator), 
                                                  rel_pos, 
                                                  relator.size()-len);
  std::cout << "Got the new relator chunk " << new_relator_chunk << "\n";
  new_relator_chunk = inverse(new_relator_chunk);
  std::cout << "Took the inverse " << new_relator_chunk << "\n";
  //if the w subword wraps around, just remove it and tack 
  //on the new chunk at the end, otherwise, replace in the middle
  if (pos + len > wL) {
    w = cyclic_subword(w, (pos+len)%wL, wL-len);
    w = w + new_relator_chunk;
    w = cyclic_subword(w, pos, w.size());
  } else {
    w.replace(pos, len, new_relator_chunk);
  }
}



std::string Surface::geodesic(std::string w) {
  //we scan the cyclic word w looking for segments in 
  //common with the relator or relator inverse
  //if all these segments have length <= 1/2 the relator length, 
  //then it's geodesic
  cyclically_reduce(w);
  std::cout << "Geodesicifying " << w << "\n";
  int too_long_length = relator.size()/2 + 1;  //if a segment is >= this length, we can shorten it
  while (true) {
    //scan to find a spot which is not geodesic
    bool did_something = false;
    for (int i=0; i<(int)w.size(); ++i) {
      //find the longest subword in common with the relator or relator inverse
      //starting at this position
      int r_position = relator_map[letter_to_alpha_ind(w[i])];
      int R_position = relator_inverse_map[letter_to_alpha_ind(w[i])];
      int r_agree_length = cyclic_word_agreement_length(w, i, relator, r_position);
      int R_agree_length = cyclic_word_agreement_length(w, i, relator_inverse, R_position);
      if (r_agree_length >= too_long_length) {
        std::cout << "Found agreement of length " << r_agree_length << " at pos " << i << "\n";
        apply_relator(w, i, r_agree_length);
        std::cout << "After replacing: " << w << "\n";
        did_something = true;
        break;
      } else if (R_agree_length >= too_long_length) {
        std::cout << "Found inverse agreement of length " << R_agree_length << " at pos " << i << "\n";
        apply_relator(w, i, R_agree_length, true);
        std::cout << "After replacing: " << w << "\n";
        did_something = true;
        break;
      }
    }
    if (!did_something) break;
  }
  return w;
}

int Surface::cyclic_order(SignedInd x, SignedInd y, SignedInd z) {
  int xpos = cyclic_order_map[x];
  int ypos = cyclic_order_map[y];
  int zpos = cyclic_order_map[z];
  while (xpos < zpos) xpos += cyclic_order.size();
  while (ypos < zpos) ypos += cyclic_order.size();
  return ( xpos < ypos ? 1 : -1 );;
}
  


bool sort_at_gen_positions(const GenPosition& gp1, const GenPosition& gp2) {
  //scan left and right to determine the cyclic order
  //this does NOT take into account anything with the relator
  //it's just for getting a good idea for nonstupid initial placement
  
  //the words can go in different directions and everything.
  //we'll change the signs so that we can think about both words 
  //going in the same direction
  std::string* w1 = gp1.word;
  int i1 = gp1.i;
  std::string* w2 = gp2.word;
  int i2 = gp2.i;
  int w1s = gp1.sign;
  int w2s = gp2.sign;
  int w1L = w1->size();
  int w2L = w2->size();
  //sanity check
  if (w1s * letter_to_alpha_ind((*w1)[i1]) != w2s * letter_to_alpha_ind((*w2)[i2])) {
    std::cout << "You're trying to sort based at two different letters\n";
    return false;
  }
  int forward_agreement = 0;
  while (true) {
    int w1_index = pos_mod(i1 + w1s*forward_agreement, w1L);
    int w2_index = pos_mod(i2 + w2s*forward_agreement, w2L);
    if (w1s*letter_to_alpha_ind( (*w1)[w1_index] ) !=
        w2s*letter_to_alpha_ind( (*w2)[w2_index] )) {
      break;
    }
    ++forward_agreement;
  }
  int w1_forward_diverge_ind = pos_mod(i1 + w1s*forward_agreement, w1L);
  int w1_forward_diverge_gen = w1s*letter_to_alpha_ind( (*w1)[w1_forward_diverge_ind] );
  int w2_forward_diverge_ind = pos_mod(i2 + w2s*forward_agreement, w2L);
  int w2_forward_diverge_gen = w2s*letter_to_alpha_ind( (*w2)[w2_forward_diverge_ind] );
  int forward_prev_index = pos_mod(i1 + w1s*(forward_agreement-1), w1L);
  int forward_back_gen = w1s*letter_to_alpha_ind( (*w1)[forward_prev_index] );
  int CO_forward = cyclic_order( w1_forward_diverge_gen,
                                 w2_forward_diverge_gen,
                                 forward_back_gen );
  
  
  int backward_agreement = 0; 
  while (true) {
    int w1_index = pos_mod(i1 + w1s*backward_agreement, w1L);
    int w2_index = pos_mod(i2 + w2s*backward_agreement, w2L);
    if (w1s*letter_to_alpha_ind( (*w1)[w1_index] ) !=
        w2s*letter_to_alpha_ind( (*w2)[w2_index] )) {
      break;
    }
    ++backward_agreement;
  }
  
}

void Surface::minimal_intersection_position(std::vector<std::string>& W, 
                                            std::vector<std::vector<int> >& positions,
                                            std::vector<int>& gen_counts) {
  //reduce the words to geodesics
  int nwords = W.size();
  for (int i=0; i<nwords; ++i) {
    W[i] = geodesic(W[i]);
  }
  
  //make the positions the right size
  positions.resize(nwords);
  for (int i=0; i<nwords; ++i) {
    positions[i].resize(W[i].size(), -1);
  }
  
  //count the gens and find where they are
  std::vector<std::vector<GenPosition> > gen_positions(ngens+1, std::vector<GenPosition>());
  GenPosition temp_gen_pos;
  temp_gen_pos.S = this;
  for (int i=0; i<nwords; ++i) {
    temp_gen_pos.w = i;
    temp_gen_pos.word = &W[i];
    for (int j=0; j<(int)W[i].size(); ++j) {
      int gen = letter_to_alpha_ind(W[i][j]);
      temp_gen_pos.i = j;
      temp_gen_pos.sign = ( gen < 0 ? -1 : 1 );
      gen_positions[ abs(gen) ].push_back(temp_gen_pos);
    }
  } 
  
  //sort the gen positions as a good first guess
  
  
  //now we go through and find every single crossing and check 
  //each of 4 directions for whether it can be reduced with a bigon move.
  //every time we do anything we must start over, because it's hard to know 
  //the effect of our moves
  
}


  
  