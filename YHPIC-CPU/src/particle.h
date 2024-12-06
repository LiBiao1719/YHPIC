#ifndef PARTICLE_H
#define PARTICLE_H
#include "otypen.h"

struct particle {

//  int             number;        // number of this particle
  int             species;       // particle species, 0=electron, 1=ion
  Scalar x,y,z;         // position and shift within one timestep
  Scalar igamma;                 // inverse gamma factor 
  Scalar ux, uy, uz;             // gamma * velocity

  struct particle *prev;         // pointer to the previous particle in this cell
  struct particle *next;         // pointer to the next particle in this cell   
};

#endif
