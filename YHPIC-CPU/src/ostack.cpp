#include "ostack.h"

//////////////////////////////////////////////////////////////////////////////////////////

Stack::Stack( )
{
  zero = new stack_member;
  hole = new stack_member;

  zero->next = hole;
}

//////////////////////////////////////////////////////////////////////////////////////////

void Stack::put_on_stack( struct particle *part,int i,int j,int k )
{

  stack_member *new_member;

  new_member = new stack_member;
  if (!new_member) { printf( "allocation error\n" ); exit(1); }

  new_member->part     = part;
  new_member->i        = i;
  new_member->j        = j;
  new_member->k        = k;
  new_member->next = zero->next;
  zero->next       = new_member;

}

//////////////////////////////////////////////////////////////////////////////////////////

void Stack::remove_from_stack( stack_member *member )
{
  zero->next = member->next;
  delete member;
}

//////////////////////////////////////////////////////////////////////////////////////////

void Stack::insert_particle( CellInfo& cells,struct particle *part,int i0,int j0,int k0 ) 
{
  int i1,j1,k1;

  i1 = (int)FLOOR(part->x);
  j1 = (int)FLOOR(part->y);
  k1 = (int)FLOOR(part->z);

  if (part->prev!=NULL) part->prev->next = part->next; // extract particle from old chain 
  else cells.CellsM[i0][j0][k0].first = part->next;
  if (part->next!=NULL) part->next->prev = part->prev;
  else cells.CellsM[i0][j0][k0].last = part->prev;

  if (cells.CellsM[i0][j0][k0].insert == part) cells.CellsM[i0][j0][k0].insert = part->next;
  // new insert mark!!

  part->prev           = NULL;                    // insert particle into the new chain
  part->next           = cells.CellsM[i1][j1][k1].first;  
  cells.CellsM[i1][j1][k1].first      = part;
  if (part->next==NULL) cells.CellsM[i1][j1][k1].last = part;
  else part->next->prev = part;

  cells.CellsM[i0][j0][k0].np[ part->species ] --;
  cells.CellsM[i0][j0][k0].npart               --;
  cells.CellsM[i1][j1][k1].np[ part->species ] ++;
  cells.CellsM[i1][j1][k1].npart               ++;

}




//////////////////////////////////////////////////////////////////////////////////////////
//eof
