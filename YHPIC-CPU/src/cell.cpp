#include "cell.h"
#include "matrix.h"

cell::cell()
{
  first=NULL;
  last =NULL;
  insert=NULL;
}

cell&  cell::operator = (cell& source)
{
  delete_particles();
  np[0] = source.np[0];
  np[1] = source.np[1];
  dens[0] = source.dens[0];
  dens[1] = source.dens[1];
  npart   = source.npart;
  first   = source.first;
  last    = source.last;
  insert  = source.insert;

  source.np[0] = 0;
  source.np[1] = 0;
//  source.dens[0] = 0.0;
//  source.dens[1] = 0.0;
  source.npart   = 0;
  source.first   = NULL;
  source.last    = NULL;
  source.insert  = NULL;
}

void cell::delete_particles()
{
  struct particle *part,*old;

  if ( first != NULL ) {
    part = first;               // delete particles
    do 
      { npart --;
	np[part->species] --;
	old  = part;
	part = part->next;
 	delete old;
      }
    while( part != NULL );
    first = NULL;
    last  = NULL;
  }
}

void cell::insert_particles(struct particle* part)
{
  int    i;                         // recieve from next domain
  struct particle *insert_pointer;

  insert_pointer = first;

  if (insert_pointer!=NULL){         // insert always in front of insert pointer
    part->next       = insert_pointer;
    part->prev       = insert_pointer->prev;
    part->next->prev = part;
    if (part->prev==NULL) first = part;
    else part->prev->next = part;
  }
  else{                              // insert always on bottom
    part->next    = NULL;
    part->prev    = last; 
    if (part->prev!=NULL) part->prev->next = part;
    else first = part;
    last    = part;
  }

  npart ++;                     // update cell's particle bookkeeping
  np[part->species] ++;  
}
