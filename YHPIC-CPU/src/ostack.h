#ifndef OSTACK_H
#define OSTACK_H

#include "particle.h"
#include "cellinfo.h"

class CellInfo;

typedef struct stack_member_struct {
  struct stack_member_struct *next;
  struct particle            *part;
  int i;
  int j;
  int k;
} stack_member;

class Stack {

 public:
                   Stack( );
  void      put_on_stack( struct particle *part,int i,int j,int k );
  void remove_from_stack( stack_member *member );
  void   insert_particle( CellInfo& cells,struct particle *part,int i,int j,int k ); 

  stack_member *zero;
  stack_member *hole;
};

#endif


