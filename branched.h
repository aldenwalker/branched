#ifndef __BRANCHED_H__
#define __BRANCHED_H__

#include <vector>
#include <string>
#include <iostream>


typedef SignedInd int;

/****************************************************************************
 * Simplex                                                                  
 ****************************************************************************/

//the boundary is listed in standard form b[i] = \hat{i}, 
//i.e. the boundary simplex on which i is omitted
class Simplex {
  private:
    int dim;
    std::vector<SignedInd> bd;
    std::vector<SignedInd> in_bd_of;
  public:
    Simplex();
    Simplex(int n);
    Simplex(std::vector<SignedInd>& bd);
    
    int dim();
    SignedInd& boundary(int i);
    SignedInd& operator[](int i);
    SignedInd& in_boundary_of(int i);
    std::string ihat_string(int i); //returns the appropriate string
    
    friend std::ostream& operator<<(std::ostream& os, Simplex& S);
};
    

/*****************************************************************************
 * Triangulation                                                             
 *****************************************************************************/
//all lists in here are ONE-BASED, so that SignedInds are easier
class Triangulation {
  private:
    std::vector<Simplex> triangles;
    std::vector<Simplex> edges;
    std::vector<Simplex> vertices;

  public:
    Triangulation();
    void read_file(std::string filename);
    void write_file(std::string filename);
    void print(ostream& os);
    Simplex& triangle(int i);
    Simplex& edge(int i);
    Simplex& vertex(int i);
  
    Triangulation branched_surface_from_vector(std::vector<int>& weights);
    Triangulation resolve_branched_surface();
};

#endif