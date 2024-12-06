#ifndef CELL_H
#define CELL_H

#include "particle.h"

class cell {
 public:
  cell();
//  int    number;                 // number of this cell
//  int    domain;                 // domain number it belongs to
    
//  Scalar x,y,z;                  // position of the left cell boundary in wavelengths
  Scalar dens[2];                // densities for each species in units n_c
  int             np[2];         // # of electrons [0] and ions [1]
  int             npart;         // # particles

  struct particle *first;        // pointer to the first particle
  struct particle *last;         // pointer to the last particle
  struct particle *insert;       // pointer to particle, in front of which new particles
                                 // have to be inserted
  cell& operator = (int c) { };
  cell&	operator = (cell& source) ;
  void delete_particles();
  void insert_particles(struct particle* part);
};


#endif


