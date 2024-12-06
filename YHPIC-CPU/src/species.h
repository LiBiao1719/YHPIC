#ifndef SPECIES_H
#define SPECIES_H

struct species {

  int    fix;                    // fixed species? 0->no, 1->yes
  int    ze;                     // charge of the micro particle in units of e
  int     m;                      // mass of the micro particle in units of m_e
  Scalar zm;                     // specific charge, z/m
  Scalar n;                      // particle density in units of n_c
  Scalar zn;                     // contribution of the particle to the charge density
                                 // in units of n_c ( = z * n )
};

#endif
