#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cmath>

#include "surface.h"
#include "graphics.h"

/****************************************************************************
 * Everything is a vector, so it is handy to write it out
 ****************************************************************************/
std::ostream& operator<<(std::ostream& os, std::vector<int>& L) {
  if (L.size() == 0) return os;
  os << "[" << L[0];
  for (int i=1; i<(int)L.size(); ++i) {
    os << "," << L[i];
  }
  os << "]";
  return os;
}

/*****************************************************************************
 * These functions translate between character generators a,b,c... 
 * and integers 1,2,3,..
 *****************************************************************************/
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

/****************************************************************************
 * These are some handy string functions
 * They are duplicated for int vectors, which represent all our words
 ****************************************************************************/
char swap_case(char c) {
  if (std::isupper(c)) return std::tolower(c);
  else if (std::islower(c)) return std::toupper(c);
  else return 0;
}

int negate(int x) {
  return -x;
}

int max(int a, int b) {
  return (a<b ? b : a);
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
  
/*****************************************************************************
 * These functions record how long the cyclic words agree at the given spot
 *****************************************************************************/
int cyclic_word_agreement_length(std::string w1, int pos1, 
                                 std::string w2, int pos2) {
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
  int max_len = max(w1L, w2L);
  while (w1[w1pos%w1L] == w2[w2pos%w2L]) {
    ++ans;
    ++w1pos;
    ++w2pos;
    if (ans == max_len) return ans;
  }
  return ans;
}

/*****************************************************************************
 * This function is similar, but it also takes a direction, which says 
 * which way to go.  Note that gens get negated as we go backwards
 *****************************************************************************/
int cyclic_word_agreement_length(std::vector<int>& w1, int pos1, int dir1, 
                                 std::vector<int>& w2, int pos2, int dir2) {
  int w1L = w1.size();
  int w2L = w2.size();
  int ans = 0;
  int w1pos = pos1;
  int w2pos = pos2;
  int max_len = max(w1L, w2L);
  while (dir1*w1[pos_mod(w1pos,w1L)] == dir2*w2[pos_mod(w2pos,w2L)]) {
    ++ans;
    w1pos += dir1;
    w2pos += dir2;
    if (ans == max_len) return ans;
  }
  return ans;
}

/*****************************************************************************
 * Construct a surface with the default polygon so the relator is [a,b][c,d]...
 *****************************************************************************/
Surface::Surface(int g, int nb, int verbose) {
  genus = g;
  nboundaries = nb;
  ngens = 2*g + nb;
  this->verbose = verbose;
  cyclic_order.resize(0);
  relator.resize(0);
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
    cyclic_order.push_back( 2*genus+nb + i + 1 ); //this is the boundary loop
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
  
  //create the polygon rational positions
  Point2d<double> prev_point(1,0);
  Point2d<double> current_point;
  double PI = 3.1415926535;
  double current_angle = 0;
  double angle_step = (2*PI)/cyclic_order.size();
  for (int i=0; i<(int)cyclic_order.size(); ++i) {
    current_angle += angle_step;
    current_point = Point2d<double>(cos(current_angle), sin(current_angle));
    int gen = cyclic_order[i];
    gen_edge_start[gen] = Point2d<Rational>(approx_rat(prev_point.x, 0.01),
                                            approx_rat(prev_point.y, 0.01));
    gen_edge_end[gen] = Point2d<Rational>(approx_rat(current_point.x, 0.01),
                                          approx_rat(current_point.y, 0.01));
    prev_point = current_point;
  }
}

/*****************************************************************************
 * Check if a gen is actually just a boundary loop placeholder
 * ***************************************************************************/
bool Surface::is_boundary_loop(int gen) {
  return (gen > 2*genus+nboundaries);
}
  
  
/*****************************************************************************
 * Check if a gen is a boundary generator
 * ***************************************************************************/
bool Surface::is_boundary_gen(int gen) {
  return (abs(gen) > 2*genus);
}

/*****************************************************************************
 * Print out a surface
 *****************************************************************************/
void Surface::print(std::ostream& os) {
  os << "Surface of genus " << genus << " with " 
                            << nboundaries << " boundary components\n";
  os << "Cyclic order action on polygon: " << cyclic_order << "\n";
  os << "Relator: " << relator << "(" << word_from_vector(relator) << ")\n";
  os << "Relator inverse: " << relator_inverse << "(" << word_from_vector(relator_inverse) << ")\n";
  os << "Gens: ";
  for (int i=1; i<=ngens; ++i) {
    os << i << "(" << alpha_ind_to_letter(i) << ") ";
  }
  os << "\n";
  os << "Boundary loops: ";
  for (int i=1; i<=nboundaries; ++i) {
    os << 2*genus + nboundaries + i << "(" << alpha_ind_to_letter(2*genus+nboundaries+i) << ") ";
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

/*****************************************************************************
 * Given a spot and length, this replaces the string using the relator, 
 * usually to make it shorter.  This uses the fact that each relator contains 
 * every letter at most once, and it assumes that the input is sensical 
 * (actually a substring of the relator)
 *****************************************************************************/
void Surface::apply_relator(std::vector<int>& w, int pos, int len, bool inv) {
  int wL = w.size();
  //get the position in the relator of the letter at pos+len
  int rel_pos = ( inv  
                  ? relator_inverse_map[ w[(pos+len-1)%wL] ]
                  : relator_map[ w[(pos+len-1)%wL] ] );
  rel_pos = (rel_pos+1)%relator.size();
  //get the cyclic subword of relator starting at rel_pos of 
  //length relator_length - len
  std::vector<int> new_relator_chunk = cyclic_subword( (inv ? relator_inverse : relator), 
                                                        rel_pos, 
                                                        relator.size()-len);
  new_relator_chunk = inverse(new_relator_chunk);
  if (verbose > 1) {
    std::cout << "Applying the relator to " << w << "\n";
    std::cout << "Replacing subword of length " << len << " at pos " << pos << "\n";
    std::cout << "Got relator position " << rel_pos << "\n";
    std::cout << "Got the bit to insert: " << new_relator_chunk << "\n";
  }
  //if the w subword wraps around, just remove it and tack 
  //on the new chunk at the end, otherwise, replace in the middle
  if (pos + len > wL) {
    w = cyclic_subword(w, (pos+len)%wL, wL-len);
    w.insert(w.end(), new_relator_chunk.begin(), new_relator_chunk.end());
    //w = cyclic_subword(w, wL-pos, w.size());
  } else {
    w.erase(w.begin()+pos, w.begin()+pos+len);
    w.insert(w.begin()+pos, new_relator_chunk.begin(), new_relator_chunk.end());
  }
}


/****************************************************************************
 * return a new string which is a geodesic rep of w
 ****************************************************************************/
std::string Surface::geodesic(std::string& w) {
  std::vector<int> wv = vector_from_word(w);
  make_geodesic(wv);
  return word_from_vector(wv);
}

/*****************************************************************************
 * turn a string (vector) into a geodesic.  If unique_geodesics=true, then 
 * it'll make an arbitrary choice about any ambiguity with an even-length
 * relator
 *****************************************************************************/
void Surface::make_geodesic(std::vector<int>& w, bool unique_geodesics) {
  //we scan the cyclic word w looking for segments in 
  //common with the relator or relator inverse
  //if all these segments have length <= 1/2 the relator length, 
  //then it's geodesic
  cyclically_reduce(w);
  if (verbose > 1) std::cout << "Geodesicifying " << w << "\n";
  
  //get the lengths we can shorten and stuff
  bool do_unique_reduction = (relator.size()%2 == 0 && unique_geodesics);
  int half_length = relator.size()/2;
  int too_long_length = relator.size()/2 + 1;
  std::vector<bool> places_to_shorten;
  std::vector<bool> places_to_shorten_inverse;
  if (do_unique_reduction) {
    places_to_shorten.resize(relator.size());
    places_to_shorten_inverse.resize(relator.size());
    if (verbose > 1) std::cout << "Doing unique reduction\n";
    for (int i=0; i<(int)relator.size(); ++i) {
      int im1 = (i==0 ? relator.size()-1 : i-1);
      places_to_shorten[i] = (relator[i] > -relator[im1]);
      places_to_shorten_inverse[i] = (relator_inverse[i] > -relator_inverse[im1]);
    }
    if (verbose>2) {
      std::cout << "Places to shorten relator: \n";
      for (int i=0; i<(int)relator.size(); ++i) {
        std::cout << i << ":" << places_to_shorten[i] << " ";
      }
      std::cout << "\n";
    }
  }
  
  while (true) {
    //scan to find a spot which is not geodesic
    bool did_something = false;
    for (int i=0; i<(int)w.size(); ++i) {
      //find the longest subword in common with the relator or relator inverse
      //starting at this position
      int r_position = relator_map[w[i]];
      int R_position = relator_inverse_map[w[i]];
      int r_agree_length = cyclic_word_agreement_length(w, 
                                                        i, 
                                                        relator, 
                                                        r_position);
      int R_agree_length = cyclic_word_agreement_length(w, 
                                                        i, 
                                                        relator_inverse, 
                                                        R_position);
      if (r_agree_length >= too_long_length
          || (do_unique_reduction 
              && r_agree_length == half_length 
              && places_to_shorten[r_position])) {
        if (verbose > 2) std::cout << "Found agreement of length " << r_agree_length << " at pos " << i << "\n";
        apply_relator(w, i, r_agree_length);
        if (verbose > 2) std::cout << "After replacing: " << w << "\n";
        did_something = true;
        break;
      } else if (R_agree_length >= too_long_length
                 || (do_unique_reduction 
                      && R_agree_length == half_length 
                      && places_to_shorten_inverse[R_position])) {
        if (verbose > 2) std::cout << "Found inverse agreement of length " << R_agree_length << " at pos " << i << "\n";
        apply_relator(w, i, R_agree_length, true);
        if (verbose > 2) std::cout << "After replacing: " << w << "\n";
        did_something = true;
        break;
      }
    }
    if (!did_something) break;
  }
}

/*****************************************************************************
 * Just check the cyclic order
 *****************************************************************************/
int Surface::cyclically_ordered(SignedInd x, SignedInd y, SignedInd z) {
  int xpos = cyclic_order_map[x];
  int ypos = cyclic_order_map[y];
  int zpos = cyclic_order_map[z];
  while (xpos < zpos) xpos += cyclic_order.size();
  while (ypos < zpos) ypos += cyclic_order.size();
  return ( xpos < ypos ? 1 : -1 );;
}
  
/*****************************************************************************
 * determine the cyclic order of two words which agree at the starting 
 * locations, so the cyclic order is well-defined.  The words can 
 * be read in either direction determined by dir1 and dir2
 *****************************************************************************/
int Surface::cyclically_ordered(std::vector<int>& w1, int start1, int dir1,
                                std::vector<int>& w2, int start2, int dir2) {
  if (dir1*w1[start1] != dir2*w2[start2]) {
    std::cout << "Error; first letters don't match\n";
    return 0;
  }
  int w1L = w1.size();
  int w2L = w2.size();
  int agreement = cyclic_word_agreement_length(w1, start1, dir1, w2, start2, dir2);
  if (verbose > 2) 
  std::cout << "I found that the agreement  of (" << w1 << "," << start1 << "," << dir1 
            << ") (" << w2 << "," << start2 << "," << dir2 << ") is " << agreement << "\n";
  if (agreement == max(w1L, w2L)) return 0;
  int w1_ind = pos_mod(start1 + dir1*agreement, w1L);
  int w2_ind = pos_mod(start2 + dir2*agreement, w2L);
  int backward_gen = -dir1*w1[pos_mod(w1_ind-dir1, w1L)];
  if (verbose > 2) std::cout << "Returning cyclic order of " << backward_gen << ", " << dir1*w1[w1_ind] << ", " << dir2*w2[w2_ind] <<"\n";
  return cyclically_ordered(backward_gen, dir1*w1[w1_ind], dir2*w2[w2_ind]);
}


/*****************************************************************************
 * construct an empty loop arrangement
 *****************************************************************************/
LoopArrangement::LoopArrangement(Surface& S, int verbose) {
  W_words.resize(0);
  W.resize(0);
  this->S = &S;
  positions_by_letter.resize(0);
  positions_by_gen.resize(0);
  this->verbose = verbose;
}

/*****************************************************************************
 * Construct a loop arrangement from a list of words
 *****************************************************************************/
LoopArrangement::LoopArrangement(Surface& S, 
                                 std::vector<std::string>& W_words, 
                                 int verbose) {
  std::vector<std::vector<int> > W_in(W_words.size());
  this->verbose = verbose;
  for (int i=0; i<(int)W_in.size(); ++i) {
    W_in[i] = vector_from_word(W_words[i]);
  }
  init_from_vectors(S, W_in);
}  
  
/*****************************************************************************
 * construct a list arrangement from a list of vector words 
 *****************************************************************************/
LoopArrangement::LoopArrangement(Surface& S, 
                                 std::vector<std::vector<int> >& W_in, 
                                 int verbose) {
  this->verbose = verbose;
  init_from_vectors(S, W_in);
}
  
/*****************************************************************************
 * this is what does the actual initialization
 *****************************************************************************/
void LoopArrangement::init_from_vectors(Surface& S, 
                                        std::vector<std::vector<int> >& W_in,
                                        bool unique_geodesics)  {
  this->S = &S;
  W = W_in;
  W_words.resize(W.size());
  
  //make the input geodesic
  for (int i=0; i<(int)W.size(); ++i) {
    S.make_geodesic(W[i], unique_geodesics);
    W_words[i] = word_from_vector(W[i]);
  }
  
  //count the gens and find where they are
  positions_by_gen = std::vector<std::vector<GenPosition> >(S.ngens+1, std::vector<GenPosition>());
  GenPosition temp_gen_pos;
  temp_gen_pos.S = &S;
  for (int i=0; i<(int)W.size(); ++i) {
    temp_gen_pos.word = &(W[i]);
    temp_gen_pos.w = i;
    for (int j=0; j<(int)W[i].size(); ++j) {
      int gen = W[i][j];
      temp_gen_pos.i = j;
      positions_by_gen[ abs(gen) ].push_back(temp_gen_pos);
    }
  } 
  
  generate_positions_by_letter();
  
  segments.resize(0);
  
}

/*****************************************************************************
 * the positions of the loop as they pass by the generators are listed by 
 * generator; this scans through and produces a new list indexed by letter
 *****************************************************************************/
void LoopArrangement::generate_positions_by_letter() {
  //make the positions vector the right size
  positions_by_letter.resize(W.size());
  for (int i=0; i<(int)W.size(); ++i) {
    positions_by_letter[i].resize(W[i].size(), -1);
  }
  
  //copy the positions over
  for (int i=1; i<=S->ngens; ++i) {
    for (int j=0; j<(int)positions_by_gen[i].size(); ++j) {
      positions_by_letter[positions_by_gen[i][j].w][positions_by_gen[i][j].i] = j;
    }
  }
}



/*****************************************************************************
 * Given two positions on two loops, say which one should be 
 * "less" in the order of the edge; this is what positions all the 
 * loops to eliminate bigons
 *****************************************************************************/
bool sort_at_gen_positions(const GenPosition& gp1, const GenPosition& gp2) {
  //scan left and right to determine the cyclic order
  //this does NOT take into account anything with the relator
  //however, it should give the minimal position IF the geodesics are unique
  
  //the words can go in different directions and everything.
  //we'll change the signs so that we can think about both words 
  //going in the same direction
  std::vector<int>& w1 = *gp1.word;
  int i1 = gp1.i;
  std::vector<int>& w2 = *gp2.word;
  int i2 = gp2.i;
  int w1s = sgn(w1[i1]);
  int w2s = sgn(w2[i2]);
  int verbose = gp1.S->verbose;
  //sanity check
  if (w1s * w1[i1] != w2s * w2[i2]) {
    std::cout << "You're trying to sort based at two different letters\n";
    return false;
  }
  int CO_forward = gp1.S->cyclically_ordered(w1, i1, w1s, w2, i2, w2s);
  if (verbose > 2) 
    std::cout << "I found that the cyclic order on (" << w1 << "," << i1 << "," << w1s 
            << ") (" << w2 << "," << i2 << "," << w2s << ") is " << CO_forward << "\n";
  if (CO_forward == 0) { 
    //this means the words are (cyclically) the same word
    //so we can sort them by saying that the word of lower 
    //index is lower in the order.  Or, if they are the same word, then 
    //the position of lower index is lower in the order 
    if (verbose > 2) std::cout << "Found a cyclic duplicate\n";
    
    if (gp1.w != gp2.w) { 
      if (verbose > 2) std::cout << "Words are different, so I'm returning " << (w1s == 1 ? (gp1.w < gp2.w) : !(gp1.w < gp2.w)) << "\n";
      return (w1s == 1 ? (gp1.w < gp2.w) : !(gp1.w < gp2.w));
    }
    if (verbose > 2) std::cout << "Words are the same, so returning " << (w1s == 1 ? (i1 < i2) : !(i1 < i2)) << "\n";
    return (w1s == 1 ? (i1 < i2) : !(i1 < i2));
  }
    
  int CO_backward = gp1.S->cyclically_ordered(w2, i2, -w2s, w1, i1, -w1s);
  if (verbose > 2) 
    std::cout << "I found that the cyclic order on (" << w2 << "," << i2 << "," << -w2s 
          << ") (" << w1 << "," << i1 << "," << -w1s << ") is " << CO_backward << "\n";
  if (CO_forward == CO_backward) {
    //(start, w1, w2), forward and (start, w2, w1) backward have the 
    //same sign; this means that the words are unlinked, so they just 
    //pull apart.  If it's positively ordered, the w1 comes before w2
    if (verbose > 2) std::cout << "These pull apart: returning " << (CO_forward > 0) << "\n";
    return (CO_forward > 0 ? true : false);
  }
  //the orders don't agree, so the words must cross; to determine the order, 
  //we'll put the crossing in the middle
  int forward_length = cyclic_word_agreement_length(w1, i1, w1s, w2, i2, w2s);
  int backward_length = cyclic_word_agreement_length(w2, i2, -w2s, w1, i1, -w1s);
  if (verbose > 2) std::cout << "Forward length " << forward_length << " and backward length " << backward_length << "\n";
  if (forward_length > backward_length || forward_length == backward_length) {
    //so we'll use the order from the forward direction
    if (verbose > 2) std::cout << "They cross, and forward agreement is longer: " << (CO_forward > 0) << "\n";
    return (CO_backward > 0 ? true : false);
  } else {
    if (verbose > 2) std::cout << "They cross, and backward agreement is longer: " << (CO_backward > 0) << "\n";
    return (CO_forward > 0 ? true : false);
  }
}


/*****************************************************************************
 * Since there are multiple positions per generator, checking whether 
 * three loop-segment-ends are positively cyclically ordered is trickier
 *****************************************************************************/
int LoopArrangement::cyclically_ordered_positions(int gen1, int pos1, 
                                                  int gen2, int pos2,
                                                  int gen3, int pos3) {
  if (gen2 == gen3) {
    if (gen1 == gen2) {
      int preorder = ((pos1 < pos2 and pos2 < pos3) ? 1 : -1);
      return sgn(gen1)*preorder;
    }
    int preorder = ((pos2 < pos3) ? 1 : -1);
    return sgn(gen2)*preorder;
  }
  if (gen1 == gen2) {
    int preorder = ((pos1 < pos2) ? 1 : -1);
    return sgn(gen1)*preorder;
  }
  if (gen1 == gen3) {
    int preorder = ((pos1 > pos3) ? 1 : -1);
    return sgn(gen1)*preorder;
  }
  return S->cyclically_ordered(gen1, gen2, gen3);
}


/*****************************************************************************
 * Check whether two loop segments cross
 *****************************************************************************/
bool LoopArrangement::check_cross(int w1, int i1, int w2, int i2) {
  int w1L = W[w1].size();
  int w2L = W[w2].size();
  int i1p1 = (i1+1)%w1L;
  int i2p1 = (i2+1)%w2L;
  int ord1 = cyclically_ordered_positions(-W[w1][i1], positions_by_letter[w1][i1],
                                          -W[w2][i2], positions_by_letter[w2][i2],
                                           W[w2][i2p1], positions_by_letter[w2][i2p1]);
  int ord2 = cyclically_ordered_positions(W[w1][i1p1], positions_by_letter[w1][i1p1],
                                          -W[w2][i2], positions_by_letter[w2][i2],
                                          W[w2][i2p1], positions_by_letter[w2][i2p1]);
  if (verbose > 3) {
    std::cout << "Checking whether " << W[w1] << "," << i1 << " and " << W[w2] << "," << i2 << " cross\n";
    std::cout << "Checking cyclic order of: (" << -W[w1][i1] << "," << positions_by_letter[w1][i1] 
                                          << "),(" << -W[w2][i2] << "," << positions_by_letter[w2][i2]
                                          << "),(" << W[w2][i2p1] << "," << positions_by_letter[w2][i2p1] << ")\n";
    std::cout << "Got cyclic order " << ord1 << "\n";
    std::cout << "Checking cyclic order of: (" << W[w1][i1p1] << "," << positions_by_letter[w1][i1p1] 
                                          << "),(" << -W[w2][i2] << "," << positions_by_letter[w2][i2]
                                          << "),(" << W[w2][i2p1] << "," << positions_by_letter[w2][i2p1] << ")\n";
    std::cout << "Got cyclic order " << ord2 << "\n";
  }
  return (ord1 != ord2);
}


/*****************************************************************************
 * find all crossings brute force
 *****************************************************************************/
int LoopArrangement::count_crossings() {
  //for every pair of segments, check if they intersect
  int num_crossings = 0;
  for (int i=0; i<(int)W.size(); ++i) {
    for (int j=0; j<(int)W[i].size(); ++j) {
      for (int k=i; k<(int)W.size(); ++k) {
        for (int m=(k==i ? j+1 : 0); m<(int)W[k].size(); ++m) {
          if (k==i && m==j) continue;
          if (check_cross(i, j, k, m)) ++num_crossings;
        }
      }
    }
  }
  return num_crossings;
}   

/*****************************************************************************
 * Compute the positions of all the segments
 *****************************************************************************/
void LoopArrangement::find_segment_coordinates() {
  //create the step vectors for every generator
  std::map<int, Point2d<Rational> > gen_edge_step;
  for (int i=0; i<(int)S->cyclic_order.size(); ++i) {
    int gen = S->cyclic_order[i];
    if (S->is_boundary_loop(abs(gen))) continue;
    gen_edge_step[gen] = S->gen_edge_end[gen] - S->gen_edge_start[gen];
    gen_edge_step[gen] = Rational(1, positions_by_gen[abs(gen)].size()+1) * gen_edge_step[gen];
  }
  
  segments = std::vector<Segment>(1);
  Segment temp_seg;
  temp_seg.S = S;
  temp_seg.LA = this;
  segments_by_letter = std::vector<std::vector<int> >(W.size());
  for (int i=0; i<(int)W.size(); ++i) {
    segments_by_letter[i].resize(W[i].size());
    temp_seg.w = i;
    for (int j=0; j<(int)W[i].size(); ++j) {
      segments_by_letter[i][j] = segments.size();
      temp_seg.i1 = j;
      temp_seg.i2 = (j+1)%int(W[i].size()); 
      int leaving_gen = -W[i][temp_seg.i1];
      int target_gen = W[i][temp_seg.i2];
      int j1_pos_ind = positions_by_letter[i][temp_seg.i1];
      int j2_pos_ind = positions_by_letter[i][temp_seg.i2];
      temp_seg.start = (  leaving_gen > 0 
                ? S->gen_edge_start[leaving_gen] + Rational(j1_pos_ind+1,1)*gen_edge_step[leaving_gen]
                : S->gen_edge_end[leaving_gen] - Rational(j1_pos_ind+1,1)*gen_edge_step[leaving_gen] );
      temp_seg.end = (  target_gen > 0 
                ? S->gen_edge_start[target_gen] + Rational(j2_pos_ind+1,1)*gen_edge_step[target_gen]
                : S->gen_edge_end[target_gen] - Rational(j2_pos_ind+1,1)*gen_edge_step[target_gen] );
      
      segments.push_back(temp_seg);
  
      //std::cout << "Letters " << j1 << " and " << j2 << " leaving gen: " << leaving_gen << 
      //             " target gen: " << target_gen << "\n";
      //std::cout << "j1_pos_ind: " << j1_pos_ind << "\n";
      //std::cout << "gen_edge_start: " << gen_edge_start[leaving_gen] << "\n";
      //std::cout << "gen_edge_end: " << gen_edge_end[leaving_gen] << "\n";
      //std::cout << "gen_edge_step: " << gen_edge_step[leaving_gen] << "\n";
      //std::cout << "j2_pos_ind: " << j2_pos_ind << "\n";
      //std::cout << "gen_edge_start: " << gen_edge_start[target_gen] << "\n";
      //std::cout << "gen_edge_end: " << gen_edge_end[target_gen] << "\n";
      //std::cout << "gen_edge_step: " << gen_edge_step[target_gen] << "\n";
    }
  }
}

/*****************************************************************************
 * Given two segments, find the rational coordinates of the location where
 * they intersect, and the rational multipliers of the whole segments
 * which get to that location
 *****************************************************************************/
void LoopArrangement::find_segment_crossing_coordinates(Segment& s1, Segment& s2, bool& do_cross,
                                                        Rational& s1_t, Rational& s2_t,
                                                        Point2d<Rational>& cross_coords) {
  Point2d<Rational> delta_vec_1 = s1.end - s1.start;
  Point2d<Rational> delta_vec_2 = s2.end - s2.start;
  
  if (verbose > 2) {
    std::cout << "Finding the crossing coordinates for segments:\n";
    std::cout << s1 << "\n";
    std::cout << s2 << "\n";
    std::cout << "Delta vector 1: " << delta_vec_1 << "\n";
    std::cout << "Delta vector 2: " << delta_vec_2 << "\n";
  }
  
  //solve the vector equation
  //s1.start + s1_t*delta_vec_1 = s2.start + s2_t*delta_vec_2
  Point2d<Rational> RHS = s2.start - s1.start;
  Rational det = delta_vec_1.x*(-delta_vec_2.y) - (-delta_vec_2.x)*delta_vec_1.y;
  if (verbose > 2) {
    std::cout << "Det: " << det << "\n";
  }
  if (det == 0) {
    do_cross = false;
    return;
  }
  Rational a11 = (-delta_vec_2.y);
  Rational a12 = -(-delta_vec_2.x);
  Rational a21 = -delta_vec_1.y;
  Rational a22 = delta_vec_1.x;
  Point2d<Rational> ans(a11*RHS.x + a12*RHS.y, a21*RHS.x + a22*RHS.y);
  ans =  det.inv() * ans;
  
  if (verbose > 2) {
    std::cout << "Found t-values " << ans.x << " and " << ans.y << "\n";
  }
  
  s1_t = ans.x;
  s2_t = ans.y;
  if (s1_t > 1 || s1_t < 0 || s2_t > 1 || s2_t < 0) {
    do_cross = false;
    return;
  }
  cross_coords = s1.start + s1_t*delta_vec_1;
  do_cross = true;
  if (verbose > 2) {
    std::cout << "Found crossing coordinates " << cross_coords << "\n";
  }
}

/*****************************************************************************
 * find the algebraic intersection number of an interval with the segments
 * ***************************************************************************/
int LoopArrangement::algebraic_segment_intersection_number(const Point2d<Rational>& p1,
                                                           const Point2d<Rational>& p2) {
  Point2d<Rational> diff = p2 - p1;
  Rational t1, t2;
  Point2d<Rational> cross_coords;
  bool do_cross;
  int ans = 0;
  for (int i=1; i<(int)segments.size(); ++i) {
    //std::cout << "Intersecting " << p1 << "->" << p2 << " and " << segments[i].start << "->" << segments[i].end << "\n";
    segment_intersection(p1, p2, segments[i].start, segments[i].end, do_cross, t1, t2, cross_coords);
    if (!do_cross) {
      //std::cout << "They don't cross\n"; 
      continue;
    }
    Point2d<Rational> diff2 = segments[i].end - segments[i].start;
    ans += (cross_z(diff, diff2) < 0 ? 1 : -1);
    //if (cross_z(diff, diff2) > 0) {
    //  std::cout << "They cross positively; winding number now " << ans << "\n";
    //} else {
    //  std::cout << "They cross negatively; winding number now " << ans << "\n";
    //}
  }
  return ans;
}


/*****************************************************************************
 * Find the data on the crossings and segments
 * This finds the segment coordinates, and then finds all 
 * the crossing locations and incidences between segments and crossings etc
 ******************************************************************************/
void LoopArrangement::find_crossing_data() {
  find_segment_coordinates();
  crossings = std::vector<Crossing>(1);
  crossings_by_coords.clear();
  for (int i=1; i<(int)segments.size(); ++i) {
    for (int j=i+1; j<(int)segments.size(); ++j) {
      bool do_cross;
      Point2d<Rational> cross_coords;
      Rational s1_t;
      Rational s2_t;
      find_segment_crossing_coordinates(segments[i], segments[j], 
                                        do_cross, s1_t, s2_t, cross_coords);
      if (do_cross != check_cross(segments[i].w, segments[i].i1,
                                  segments[j].w, segments[j].i1)) {
        std::cout << "Crossing error\n";
      }
      if (!do_cross) continue;
      //check if this crossing already exists
      std::map<Point2d<Rational>,int, dictionary_order>::iterator it = crossings_by_coords.find(cross_coords);
      int crossing_index;
      if (it == crossings_by_coords.end()) {
        Crossing temp_cross;
        temp_cross.S = S;
        temp_cross.LA = this;
        temp_cross.coords = cross_coords;
        temp_cross.segments.resize(0);
        crossings.push_back(temp_cross);
        crossing_index = crossings.size()-1;
        crossings_by_coords[cross_coords] = crossing_index;
      } else {
        crossing_index = it->second;
      }
      //add the segments to the crossings
      //these will get sorted later
      crossings[crossing_index].segments.push_back(i);
      crossings[crossing_index].segments.push_back(-i);
      crossings[crossing_index].segments.push_back(j);
      crossings[crossing_index].segments.push_back(-j);
      
      //add the crossing to the segments
      segments[i].crossings.push_back( std::make_pair(s1_t, crossing_index) );
      segments[j].crossings.push_back( std::make_pair(s2_t, crossing_index) );
    }
  }
  //now all the segments and crossings know about each other
  //but they are not sorted
  //sort the crossing appearances on the segments
  for (int i=1; i<(int)segments.size(); ++i) {
    std::sort(segments[i].crossings.begin(), segments[i].crossings.end());
  }
  //for each crossing, sort the incident segments
  for (int i=1; i<(int)crossings.size(); ++i) {
    std::vector<std::pair<double, int> > segs_with_angles(crossings[i].segments.size());
    for (int j=0; j<(int)crossings[i].segments.size(); ++j) {
      Point2d<Rational> other_end;
      int seg_ind = crossings[i].segments[j];
      if (seg_ind < 0) {
        other_end = segments[abs(seg_ind)].start;
      } else {
        other_end = segments[seg_ind].end;
      }
      Point2d<Rational> delta_v = other_end - crossings[i].coords;
      segs_with_angles[j].first = atan2(delta_v.y.get_d(), delta_v.x.get_d());
      segs_with_angles[j].second = seg_ind;
    }
    std::sort(segs_with_angles.begin(), segs_with_angles.end());
    for (int j=0; j<(int)crossings[i].segments.size(); ++j) {
      crossings[i].segments[j] = segs_with_angles[j].second;
    }
  }
      
  
}
        

/****************************************************************************
 * basically just by sorting, put the loops into minimal position
 ****************************************************************************/
void LoopArrangement::minimal_position() {
  
  //Unnecessary:
  //find_all_crossings();
  //std::cout << "Showing unsorted positions\n";
  //std::cout << "There are " << crossings.size() << " crossings\n";
  //show();
  //////////////////////
  
  //reget the geodesics and everything, except force uniqueness
  init_from_vectors(*S, W, true);
  
  if (verbose > 2) {
    std::cout << "re-inited with unique geodesics:\n";
    print(std::cout);
  }
  
  //sort the gen positions; now this will be minimal position
  for (int i=1; i<=S->ngens; ++i) {
    std::sort(positions_by_gen[i].begin(), 
              positions_by_gen[i].end(), 
              sort_at_gen_positions);
  }
  
  if (verbose > 1) {
    std::cout << "Sorted the gen positions\n";
  }
  
  //regenerate the letter positions
  generate_positions_by_letter();
  
  //unnecessary
  //find_all_crossings();
  //std::cout << "Showing sorted positions\n";
  //std::cout << "There are " << crossings.size() << " crossings\n";
  //show();
  ///////////////////////////

}


/*****************************************************************************
 * Produce a cellulation which gives the complementary regions of 
 * the loop arrangement
 * ***************************************************************************/
Cellulation LoopArrangement::cellulation_from_loops() {
  
  //make sure that the crossing data has been computed
  find_crossing_data();
  
  Cellulation C;
  C.vertices.resize(1);
  C.edges.resize(1);
  C.cells.resize(1);
  
  if (verbose > 1) {
    std::cout << "Creating cellulation from loops\n";
  }
  
  //every crossing becomes a vertex
  for (int i=1; i<(int)crossings.size(); ++i) {
    Vertex temp_vert;
    temp_vert.in_bd_of = crossings[i].segments;
    C.vertices.push_back(temp_vert);
  }

  if (verbose > 2) {
    std::cout << "Created " << C.vertices.size()-1 << " vertices from crossings\n";
  }
  
  //make an edge for every chunk of every segment
  std::vector<std::vector<int> > edges_along_segments(segments.size()+1);
  for (int i=1; i<(int)segments.size(); ++i) {
    edges_along_segments[i].resize(0);
    Edge temp_edge;
    temp_edge.in_bd_pos.resize(0);
    temp_edge.in_bd_neg.resize(0);
    temp_edge.two_sided = false; //all of these are not two sided
    temp_edge.boundary_loop = false;
    for (int j=0; j<(int)segments[i].crossings.size()+1; ++j) {
      edges_along_segments[i].push_back(C.edges.size());
      C.edges.push_back(temp_edge);
    }
  }
  
  //now scan through the crossings and fill in the details
  for (int i=1; i<(int)segments.size(); ++i) {
    
    //we don't know the vertices for the start and end, 
    //but we do know the locations
    C.edges[edges_along_segments[i][0]].start_pos = segments[i].start;
    C.edges[edges_along_segments[i][edges_along_segments[i].size()-1]].end_pos = segments[i].end;
    
    for (int j=0; j<(int)segments[i].crossings.size(); ++j) {
      int this_vert = segments[i].crossings[j].second;
      C.edges[edges_along_segments[i][j]].end = this_vert;
      C.edges[edges_along_segments[i][j]].end_pos = crossings[this_vert].coords;
      C.edges[edges_along_segments[i][j+1]].start = this_vert;
      C.edges[edges_along_segments[i][j+1]].start_pos = crossings[this_vert].coords;
      //fill in the correct edges -- note we need to look in the crossing data
      //because we are editing the vertices all the time
      for (int k=0; k<(int)crossings[this_vert].segments.size(); ++k) {
        if (crossings[this_vert].segments[k] == i) {
          C.vertices[this_vert].in_bd_of[k] = edges_along_segments[i][j+1];
        } else if (crossings[this_vert].segments[k] == -i) {
          C.vertices[this_vert].in_bd_of[k] = -edges_along_segments[i][j];
        }
      }
    }
  }
  
  if (verbose>2) {
    std::cout << "Added " << C.edges.size()-1 << " edges from segments\n";
  }
  
  //now we need to add edges and vertices all along the polygon boundary
  
  //make the big vertex at the center of the polygon edges
  //the edges outgoing go in the reversed order of the normal relator
  //e.g. DCdcBAba
  int main_vertex_ind = C.vertices.size();
  Vertex temp_main_vert;
  temp_main_vert.in_bd_of.resize(S->relator.size());
  C.vertices.push_back(temp_main_vert);
  if (verbose > 2) {
    std::cout << "The main vertex index is " << main_vertex_ind << "\n";
  }
  
  //make a vertex for every boundary
  std::map<int, int> vertices_for_boundary_gens;
  for (int i=2*S->genus+1; i<=2*S->genus+S->nboundaries; ++i) {
    Vertex temp_vert;
    temp_vert.in_bd_of.resize(0);
    vertices_for_boundary_gens[i] = C.vertices.size();
    C.vertices.push_back(temp_vert);
  }
  
  //add a bunch of edge placeholders for every generator
  std::vector<std::vector<int> > edges_along_gens(S->ngens+1);
  for (int i=1; i<=S->ngens; ++i) {
    edges_along_gens[i].resize(0);
    Edge temp_edge;
    temp_edge.in_bd_pos.resize(0);
    temp_edge.in_bd_neg.resize(0);
    temp_edge.boundary_loop = false;
    for (int j=0; j<(int)positions_by_gen[i].size()+1; ++j) {
      edges_along_gens[i].push_back(C.edges.size());
      C.edges.push_back(temp_edge);
    }
  }
  
  //now go back through and deal with the vertices
  for (int i=1; i<=S->ngens; ++i) {
    //for the first edge in the gen, add it to the initial vertex
    int first_edge = edges_along_gens[i][0];
    C.edges[first_edge].start = main_vertex_ind;
    temp_main_vert.in_bd_of[S->relator.size() - (S->relator_map[i]+1)] = first_edge;
    C.edges[first_edge].two_sided = true;
    C.edges[first_edge].start_pos = S->gen_edge_start[i];
    C.edges[first_edge].start_neg = S->gen_edge_end[-i];
    
    //for the final edge, add it also to the final vertex
    //(this might be the same edge)
    //if it's a boundary generator, attach it instead to the vertex for that boundary
    int last_edge = edges_along_gens[i][edges_along_gens[i].size()-1];
    if (S->is_boundary_gen(i)) {
      C.edges[last_edge].end = vertices_for_boundary_gens[i];
      C.edges[last_edge].two_sided = true;
      C.edges[last_edge].end_pos = S->gen_edge_end[i];
      C.edges[last_edge].end_neg = S->gen_edge_start[-i];
      C.vertices[vertices_for_boundary_gens[i]].in_bd_of.push_back(-last_edge);
    } else {
      C.edges[last_edge].end = main_vertex_ind;
      temp_main_vert.in_bd_of[S->relator.size() - (S->relator_map[-i]+1)] = -last_edge;
      C.edges[last_edge].two_sided = true;
      C.edges[last_edge].end_pos = S->gen_edge_end[i];
      C.edges[last_edge].end_neg = S->gen_edge_start[-i];
    }
    
    //fill in the vertices
    for (int j=0; j<(int)edges_along_gens[i].size()-1; ++j) {
      int word = positions_by_gen[i][j].w;
      int letter = positions_by_gen[i][j].i;
      int sign = sgn(W[word][letter]);
      int prev_letter = pos_mod(letter-1, W[word].size());
      int out_segment = segments_by_letter[word][letter];
      int in_segment = segments_by_letter[word][prev_letter];
      int out_edge = edges_along_segments[out_segment][0];
      int in_edge = edges_along_segments[in_segment][edges_along_segments[in_segment].size()-1];     
      Point2d<Rational> out_segment_start = segments[out_segment].start;
      Point2d<Rational> in_segment_end = segments[in_segment].end;
      
      //connect all the edges
      int this_vert = C.vertices.size();
      C.edges[edges_along_gens[i][j]].two_sided = true;
      C.edges[edges_along_gens[i][j]].end = this_vert;
      C.edges[edges_along_gens[i][j]].end_pos = (sign > 0 ? in_segment_end : out_segment_start);
      C.edges[edges_along_gens[i][j]].end_neg = (sign > 0 ? out_segment_start : in_segment_end);
      C.edges[edges_along_gens[i][j+1]].two_sided = true;
      C.edges[edges_along_gens[i][j+1]].start = this_vert;
      C.edges[edges_along_gens[i][j+1]].start_pos = (sign > 0 ? in_segment_end : out_segment_start);
      C.edges[edges_along_gens[i][j+1]].start_neg = (sign > 0 ? out_segment_start : in_segment_end);
      C.edges[out_edge].start = this_vert;
      C.edges[in_edge].end = this_vert;
      
      //create the new vertex
      Vertex temp_vert;
      temp_vert.in_bd_of.resize(4);
      temp_vert.in_bd_of[0] = (sign > 0 ? -in_edge : out_edge);
      temp_vert.in_bd_of[1] = -edges_along_gens[i][j];
      temp_vert.in_bd_of[2] = (sign > 0 ? out_edge : -in_edge);
      temp_vert.in_bd_of[3] = edges_along_gens[i][j+1];
      C.vertices.push_back(temp_vert);
    }
  }
  
  //copy the temp main vert boundary so it's correct
  C.vertices[main_vertex_ind].in_bd_of = temp_main_vert.in_bd_of;

  //we just need to add one more edge for each boundary
  for (int i=2*S->genus+1; i<=2*S->genus+S->nboundaries; ++i) {
    int this_boundary_vert = vertices_for_boundary_gens[i];
    int this_edge = C.edges.size();
    Edge temp_edge;
    temp_edge.start = temp_edge.end = this_boundary_vert;
    C.vertices[this_boundary_vert].in_bd_of.push_back(-this_edge);
    C.vertices[this_boundary_vert].in_bd_of.push_back(this_edge);
    temp_edge.two_sided = false;
    temp_edge.boundary_loop = true;
    temp_edge.start_pos = S->gen_edge_end[i];
    temp_edge.end_pos = S->gen_edge_start[-i];
    C.edges.push_back(temp_edge);
  }
  
  
  //the edges along the loops form the loops of the cellulation
  C.loops.resize(0);
  for (int i=0; i<(int)W.size(); ++i) {
    std::vector<int> temp_loop(0);
    for (int j=0; j<(int)W[i].size(); ++j) {
      int seg_ind = segments_by_letter[i][j];
      temp_loop.insert(temp_loop.end(), edges_along_segments[seg_ind].begin(),
                                        edges_along_segments[seg_ind].end());
    }
    C.loops.push_back(temp_loop);
  }
  
  //now we need to go through and find the cells
  Cell temp_cell;
  temp_cell.computed_winding_number = false;
  std::vector<std::pair<bool,bool> > done_edge(C.edges.size());
  for (int i=1; i<(int)C.edges.size(); ++i) done_edge[i] = std::make_pair(false, false);
  while (true) {
    int start_edge = 0;
    for (int i=1; i<(int)done_edge.size(); ++i) {
      if (!done_edge[i].first) {
        start_edge = i;
        done_edge[i].first = true;
        break;
      } else if (!done_edge[i].second) {
        start_edge = -i;
        done_edge[i].second = true;
        break;
      }
    }
    if (start_edge == 0) break;
    temp_cell.sign = 1;
    temp_cell.bd = C.follow_edge(start_edge);
    temp_cell.contains_boundary = false;
    temp_cell.coords = Point2d<Rational>(0,0);
    for (int i=0; i<(int)temp_cell.bd.size(); ++i) {
      int e = temp_cell.bd[i];
      if (e>0) {
        C.edges[e].in_bd_pos.push_back(C.cells.size());
        done_edge[e].first = true;
        temp_cell.coords = temp_cell.coords + C.edges[e].start_pos;
      } else {
        C.edges[-e].in_bd_neg.push_back(C.cells.size());
        done_edge[-e].second = true;
        if (C.edges[-e].two_sided == true) {
          temp_cell.coords = temp_cell.coords + C.edges[-e].end_neg;
        } else {
          temp_cell.coords = temp_cell.coords + C.edges[-e].end_pos;
        } 
      }
      if (C.edges[abs(temp_cell.bd[i])].boundary_loop) {
        temp_cell.contains_boundary = true;
      }
    }
    temp_cell.coords = Rational(1,temp_cell.bd.size())*temp_cell.coords;
    C.cells.push_back(temp_cell);
  }
  
  return C;
}
  









/*****************************************************************************
 * print out a segment
 *****************************************************************************/
std::ostream& operator<<(std::ostream& os, Segment& s) {
  os << "Segment " << s.w << "," << s.i1 << "; start: " << s.start << " end: " << s.end;
  os << "; Crossings: ";
  os << std::setprecision(3);
  for (int i=0; i<(int)s.crossings.size(); ++i) {
    os << s.crossings[i].first.get_d() << "," << s.crossings[i].second << "; ";
  }
  os << std::setprecision(6);
  return os;
}
  
/*****************************************************************************
 * print out a crossing
 *****************************************************************************/
std::ostream& operator<<(std::ostream& os, Crossing& c) {
  os << "Crossing at " << c.coords << " segments: " << c.segments;
  return os;
}


/*****************************************************************************
 * Draw a loop arrangement to an X11 window
 *****************************************************************************/
void LoopArrangement::show(Cellulation* C) {
  
  //start an 800x900 graphics windows with range [-1,1]x[-1.1,1]
  Point2d<float> translate(410,480);
  XGraphics X(820, 890, (float)400, translate);
  
  //make sure we've computed the endpoints of each segment
  if (segments.size() < 2) {
    find_segment_coordinates();
  }
  
  //if we have a cellulation, draw it now so everything goes over it
  if (C != NULL) {
    (*C).draw_to_xgraphics(X);
  }
  
  //get the colors for each word
  const char* color_list_arr[] = {"red", "blue", "limegreen", "gold", "cyan", "fuchsia", "orange"};
  std::vector<std::string> color_list(color_list_arr, color_list_arr+7);
  std::vector<int> word_colors(W.size());
  for (int i=0; i<(int)W.size(); ++i) {
    word_colors[i] = X.get_color(color_list[i%color_list.size()].c_str());
  }
  
  //draw the polygon
  int black_color = X.get_color("black");
  for (int i=0; i<(int)S->cyclic_order.size(); ++i) {
    int gen = S->cyclic_order[i];
    if (verbose>1) {
      std::cout << "Polygon edge " << gen << " has coordinates " 
               << S->gen_edge_start[gen] << " and " << S->gen_edge_end[gen] << "\n";
    }
    Point2d<float> ge_start(S->gen_edge_start[gen].x.get_d(),
                            S->gen_edge_start[gen].y.get_d());
    Point2d<float> ge_end(S->gen_edge_end[gen].x.get_d(),
                          S->gen_edge_end[gen].y.get_d());
    if (S->is_boundary_loop(gen)) {
      X.draw_line(ge_start, ge_end, black_color, 2);
    } else {
      X.draw_line(ge_start, ge_end, black_color, 1);
    }
    //and the label
    Point2d<float> label_spot = (float)(1.05*0.5)*(ge_start + ge_end);
    if (S->is_boundary_loop(gen)) {
      // I guess do nothing?
    } else {
      std::string label(1,alpha_ind_to_letter(gen));
      X.draw_text_centered(label_spot, label, black_color);
    } 
  }
  
  //print out the key
  Point2d<float> box_pos(-0.90, -0.9);
  Point2d<float> word_pos(-0.88, -0.9);
  Point2d<float> step(0, -0.03);
  for (int i=0; i<(int)W.size(); ++i) {
    X.draw_box_radius(box_pos, 0.01, word_colors[i]);
    X.draw_text(word_pos, W_words[i], black_color);
    box_pos = box_pos + step;
    word_pos = word_pos + step;
  }
  
  //draw the loops
  for (int i=1; i<(int)segments.size(); ++i) {
    Point2d<float> start(segments[i].start.x.get_d(), segments[i].start.y.get_d());
    Point2d<float> end(segments[i].end.x.get_d(), segments[i].end.y.get_d());
    int col = word_colors[segments[i].w];
    X.draw_arrowed_labeled_line(start, end, col, 2, std::string(""));
    //X.draw_line(start, end, col, 2);
    //std::stringstream edge_label_s;
    //edge_label_s << segments[i].i1;
    //std::string edge_label = edge_label_s.str();
    //X.draw_text_centered(start + (float)0.5*(end-start), edge_label, black_color);
  }
  
  //draw the crossings, if we have them
  for (int i=1; i<(int)crossings.size(); ++i) {
    Point2d<float> center(crossings[i].coords.x.get_d(),
                          crossings[i].coords.y.get_d());
    X.draw_circle(center, 4, black_color);
  }
  
  std::string key_press;
  X.flush();
  key_press = X.wait_for_key();
}
  
/*****************************************************************************
 * Print the data of a loop arrangement to the screen
 *****************************************************************************/
void LoopArrangement::print(std::ostream& os) {
  os << "Loop arrangement of " << W.size() << " loops\n";
  os << "Verbose: " << verbose << "\n";
  for (int i=0; i<(int)W.size(); ++i) {
    os << i << ": " << W[i] << "(" << W_words[i] << ")\n";
  }
  std::cout << "Positions by gen:\n";
  for (int i=1; i<=S->ngens; ++i) {
    std::cout << i << ": ";
    for (int j=0; j<(int)positions_by_gen[i].size(); ++j) {
      std::cout << "(" << positions_by_gen[i][j].w << "," << positions_by_gen[i][j].i << ") ";
    }
    std::cout << "\n";
  }
  if (segments.size() > 1) {
    os << "Segments:\n";
    for (int i=1; i<(int)segments.size(); ++i) {
      os << i << ": " << segments[i] << "\n";
    }
  }
  if (crossings.size() > 1) {
    os << "Crossings:\n";
    for (int i=1; i<(int)crossings.size(); ++i) {
      os << i << ": " << crossings[i] << "\n";
    }
  }
    
}
  
  
  
  
  
  
  
  
  
  