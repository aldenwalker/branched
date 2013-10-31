#include <sstream>
#include <algorithm>

#include "surface.h"


std::ostream& operator<<(std::ostream& os, std::vector<int>& L) {
  if (L.size() == 0) return os;
  os << "[" << L[0];
  for (int i=1; i<(int)L.size(); ++i) {
    os << "," << L[i];
  }
  os << "]";
  return os;
}

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

std::string word_from_vector(std::vector<int>& w) {
  std::string ans("");
  for (int i=0; i<(int)w.size(); ++i) {
    ans += alpha_ind_to_letter(w[i]);
  }
  return ans;
}

std::vector<int> vector_from_word(std::string& w) {
  std::vector<int> ans(w.size());
  for (int i=0; i<(int)w.size(); ++i) {
    ans[i] = letter_to_alpha_ind(w[i]);
  }
  return ans;
}

char swap_case(char c) {
  if (std::isupper(c)) return std::tolower(c);
  else if (std::islower(c)) return std::toupper(c);
  else return 0;
}

int negate(int x) {
  return -x;
}

std::string inverse(std::string& w) {
  std::string ans = w;
  std::reverse(ans.begin(), ans.end());
  std::transform(ans.begin(), ans.end(), ans.begin(), swap_case);
  return ans;
}

void invert(std::string& w) {
  std::reverse(w.begin(), w.end());
  std::transform(w.begin(), w.end(), w.begin(), swap_case);
}

std::vector<int> inverse(std::vector<int>& w) {
  std::vector<int> ans = w;
  std::reverse(ans.begin(), ans.end());
  std::transform(ans.begin(), ans.end(), ans.begin(), negate);
  return ans;
}

void invert(std::vector<int>& w) {
  std::reverse(w.begin(), w.end());
  std::transform(w.begin(), w.end(), w.begin(), swap_case);
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

void reduce(std::vector<int>& w) {
int i=0; 
  while (i < (int)w.size()-1) {
    if (w[i] == -w[i+1]) {
      w.erase(w.begin()+i, w.begin()+i+2);
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

void cyclically_reduce(std::vector<int>& w) {
  reduce(w);
  int wL = w.size();
  int i=0;
  while (i < wL/2 && w[i] == -w[wL-1-i]) {
    ++i;
  }
  w.erase(w.begin()+(wL-1-(i-1)), w.end());
  w.erase(w.begin(),w.begin()+i);
}

  
std::string cyclic_subword(std::string& w, int pos, int len) {
  int wL = w.size();
  if (pos+len > wL) {
    return w.substr(pos, wL-pos) + cyclic_subword(w, 0, len-(wL-pos));
  } else {
    return w.substr(pos, len);
  }
}

std::vector<int> cyclic_subword(std::vector<int>& w, int pos, int len) {
  int wL = w.size();
  std::vector<int> ans;
  if (pos+len > wL) {
    ans.insert(ans.end(), w.begin()+pos, w.end());
    ans.insert(ans.end(), w.begin(), w.begin()+(len-(wL-pos)));
  } else {
    ans.insert(ans.end(), w.begin()+pos, w.begin()+pos+len);
  }
  return ans;
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

int cyclic_word_agreement_length(std::vector<int>& w1, int pos1, 
                                 std::vector<int>& w2, int pos2) {
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

int cyclic_word_agreement_length(std::vector<int>& w1, int pos1, int dir1, 
                                 std::vector<int>& w2, int pos2, int dir2) {
  int w1L = w1.size();
  int w2L = w2.size();
  int ans = 0;
  int w1pos = pos1;
  int w2pos = pos2;
  while (dir1*w1[pos_mod(w1pos,w1L)] == dir2*w2[pos_mod(w2pos,w2L)]) {
    ++ans;
    w1pos += dir1;
    w2pos += dir2;
    if (ans == w1L || ans == w2L) return ans;
  }
  return ans;
}

Surface::Surface(int g, int nb) {
  genus = g;
  nboundaries = nb;
  ngens = 2*g + nb;
  cyclic_order.resize(0);;
  relator.resize(0);;
  //go around the normal closed surface part
  for (int i=0; i<genus; ++i) {
    cyclic_order.push_back( 2*i+1 );
    cyclic_order.push_back( -((2*i+1)+1) );
    cyclic_order.push_back( -(2*i+1) );
    cyclic_order.push_back( (2*i+1)+1 );
    relator.push_back( 2*i+1 );
    relator.push_back( (2*i+1)+1 );
    relator.push_back( -(2*i+1) );
    relator.push_back( -((2*i+1)+1) );
  }
  //go over the boundary components
  for (int i=0; i<nboundaries; ++i) {
    cyclic_order.push_back( 2*genus+i+1 );
    cyclic_order.push_back( -(2*genus+i+1) );
    relator.push_back( 2*genus+i+1 );
  }
  relator_inverse = inverse(relator);
  
  //create the maps for fast index access
  cyclic_order_map.clear();
  relator_map.clear();
  relator_inverse_map.clear();
  for (int i=0; i<(int)cyclic_order.size(); ++i) {
    cyclic_order_map[ cyclic_order[i] ] = i;
  }
  for (int i=0; i<(int)relator.size(); ++i) {
    relator_map[ relator[i] ] = i;
    relator_inverse_map[ relator_inverse[i] ] = i;
  }
}
  
  
void Surface::print(std::ostream& os) {
  os << "Surface of genus " << genus << " with " << nboundaries << " boundary components\n";
  os << "Cyclic order action on polygon: " << cyclic_order << "\n";
  os << "Relator: " << relator << "(" << word_from_vector(relator) << ")\n";
  os << "Relator inverse: " << relator_inverse << "(" << word_from_vector(relator_inverse) << ")\n";
  os << "Gens: ";
  for (int i=1; i<=ngens; ++i) {
    os << alpha_ind_to_letter(i) << " ";
  }
  os << "\n";
  os << "Cyclic order map: ";
  for (std::map<int,int>::iterator it=cyclic_order_map.begin(); 
       it != cyclic_order_map.end(); 
       it++) {
    os << it->first << "(" << alpha_ind_to_letter(it->first) << "):" << it->second << ", ";
  }
  os << "\n";
  os << "Relator map: ";
  for (std::map<int,int>::iterator it=relator_map.begin(); 
       it != relator_map.end(); 
       it++) {
    os << it->first << "(" << alpha_ind_to_letter(it->first) << "):" << it->second << ", ";
  }
  os << "\n";
  os << "Relator inverse map: ";
  for (std::map<int,int>::iterator it=relator_inverse_map.begin(); 
       it != relator_inverse_map.end(); 
       it++) {
    os << it->first << "(" << alpha_ind_to_letter(it->first) << "):" << it->second << ", ";
  }
  os << "\n";
}


void Surface::apply_relator(std::vector<int>& w, int pos, int len, bool inv) {
  int wL = w.size();
  //get the position in the relator of the letter at pos+len
  int rel_pos = ( inv  
                  ? relator_inverse_map[ w[(pos+len-1)%wL] ]
                  : relator_map[ w[(pos+len-1)%wL] ] );
  rel_pos = (rel_pos+1)%relator.size();
  //std::cout << "Replacing subword of length " << len << " at pos " << pos << "\n";
  //std::cout << "Got relator position " << rel_pos << "\n";
  //get the cyclic subword of relator starting at rel_pos of 
  //length relator_length - len
  std::vector<int> new_relator_chunk = cyclic_subword( (inv ? relator_inverse : relator), 
                                                        rel_pos, 
                                                        relator.size()-len);
  std::cout << "Got the new relator chunk " << new_relator_chunk << "\n";
  new_relator_chunk = inverse(new_relator_chunk);
  std::cout << "Took the inverse " << new_relator_chunk << "\n";
  //if the w subword wraps around, just remove it and tack 
  //on the new chunk at the end, otherwise, replace in the middle
  if (pos + len > wL) {
    w = cyclic_subword(w, (pos+len)%wL, wL-len);
    w.insert(w.end(), new_relator_chunk.begin(), new_relator_chunk.end());
    w = cyclic_subword(w, wL-pos, w.size());
  } else {
    w.erase(w.begin()+pos, w.begin()+pos+len);
    w.insert(w.begin()+pos, new_relator_chunk.begin(), new_relator_chunk.end());
  }
}



std::string Surface::geodesic(std::string& w) {
  std::vector<int> wv = vector_from_word(w);
  std::vector<int> wv_g = geodesic(wv);
  return word_from_vector(wv_g);
}

std::vector<int> Surface::geodesic(std::vector<int>& w) {
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
      int r_position = relator_map[w[i]];
      int R_position = relator_inverse_map[w[i]];
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

int Surface::cyclically_ordered(SignedInd x, SignedInd y, SignedInd z) {
  int xpos = cyclic_order_map[x];
  int ypos = cyclic_order_map[y];
  int zpos = cyclic_order_map[z];
  while (xpos < zpos) xpos += cyclic_order.size();
  while (ypos < zpos) ypos += cyclic_order.size();
  return ( xpos < ypos ? 1 : -1 );;
}
  
int Surface::cyclically_ordered(std::vector<int>& w1, int start1, int dir1,
                                std::vector<int>& w2, int start2, int dir2) {
  if (dir1*w1[start1] != dir2*w2[start2]) {
    std::cout << "Error; first letters don't match\n";
    return 0;
  }
  int w1L = w1.size();
  int w2L = w2.size();
  int agreement = cyclic_word_agreement_length(w1, start1, dir1, w2, start2, dir2);
  int w1_ind = pos_mod(start1 + dir1*agreement, w1L);
  int w2_ind = pos_mod(start2 + dir2*agreement, w2L);
  int backward_gen = dir1*w1[pos_mod(w1_ind-1, w1L)];
  return cyclically_ordered(backward_gen, dir1*w1[w1_ind], dir2*w2[w2_ind]);
}


bool sort_at_gen_positions(GenPosition& gp1, GenPosition& gp2) {
  //scan left and right to determine the cyclic order
  //this does NOT take into account anything with the relator
  //it's just for getting a good idea for nonstupid initial placement
  
  //the words can go in different directions and everything.
  //we'll change the signs so that we can think about both words 
  //going in the same direction
  std::vector<int>& w1 = *gp1.word;
  int i1 = gp1.i;
  std::vector<int>& w2 = *gp2.word;
  int i2 = gp2.i;
  int w1s = sgn(w1[i1]);
  int w2s = sgn(w2[i2]);
  //sanity check
  if (w1s * w1[i1] != w2s * w2[i2]) {
    std::cout << "You're trying to sort based at two different letters\n";
    return false;
  }
  int CO_forward = gp1.S->cyclically_ordered(w1, i1, w1s, w2, i2, w2s);
  int CO_backward = gp1.S->cyclically_ordered(w2, i2, -w2s, w1, i1, -w1s);
  if (CO_forward == CO_backward) {
    //(start, w1, w2), forward and (start, w2, w1) backward have the 
    //same sign; this means that the words are unlinked, so they just 
    //pull apart.  If it's positively ordered, the w1 comes before w2
    return (CO_forward > 0 ? true : false);
  }
  //the orders don't agree, so the words must cross; to determine the order, 
  //we'll put the crossing in the middle
  int forward_length = cyclic_word_agreement_length(w1, i1, w1s, w2, i2, w2s);
  int backward_length = cyclic_word_agreement_length(w2, i1, -w2s, w1, i1, -w1s);
  if (forward_length > backward_length || forward_length == backward_length) {
    //so we'll use the order from the forward direction
    return (CO_forward > 0 ? true : false);
  } else {
    return (CO_backward > 0 ? true : false);
  }
}

void Surface::minimal_intersection_position(std::vector<std::string>& W_words, 
                                            std::vector<std::vector<int> >& positions,
                                            std::vector<int>& gen_counts) {
  //turn the words into gen lists
  std::vector<std::vector<int> > W(W_words.size());
  for (int i=0; i<(int)W.size(); ++i) {
    W[i] = vector_from_word(W_words[i]);
  }
  
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
    temp_gen_pos.word = &W[i];
    for (int j=0; j<(int)W[i].size(); ++j) {
      int gen = letter_to_alpha_ind(W[i][j]);
      temp_gen_pos.i = j;
      gen_positions[ abs(gen) ].push_back(temp_gen_pos);
    }
  } 
  
  //sort the gen positions as a good first guess
  
  
  //now we go through and find every single crossing and check 
  //each of 4 directions for whether it can be reduced with a bigon move.
  //every time we do anything we must start over, because it's hard to know 
  //the effect of our moves
  
}


  
  