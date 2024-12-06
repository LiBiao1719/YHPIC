//////////////////////////////////////////////////////////////////////////////////////////
//
// tools for reading ascii files
// e.g. 'namelist' input
//
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef READFILE_H
#define READFILE_H

#define MAX_LINE_LENGTH 1000
#define MAX_COL 50

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class readfile {
 public:

  readfile();

  void           openinput( char* );
  void          closeinput( void );
  int             setinput( char* );
  char*           getinput( char* );
  char*             setget( char*, char* );
  int        read_one_line( void );
  void      write_one_line( void );

  int       file_structure( int *col_min, int *col_max, int *rows );
  int            read_line( int *narg, char **arg );
  int             read_col( double *data, int col, int max_rows );

  void           copy_file( char *input, char *output );
  int        compare_files( char *name1, char *nsma2 );

  char**           cmatrix( long nrh, long nch );
  void        free_cmatrix( char **m );
  char*            cvector( long nch );
  void        free_cvector( char *m );

 private:

  int  already_open;
  FILE *fd;
  char *buffer;
  char *result;
};

#endif

