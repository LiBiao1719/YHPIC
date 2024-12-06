#ifndef MATRIX_INCLUDED
#define MATRIX_INCLUDED
#include <stdio.h>
#include <assert.h>
#include "vector.h"

#ifndef FALSE
#define FALSE (1==0)
#endif
#ifndef TRUE
#define TRUE (0==0)
#endif

// 
// vector/plane/matrix with distributed file access functionality
//

//#define INBOUNDS(a,b,i) ((a[i]>=b[i]) && (a[i]<=b[i+1]) && (a[i+1]>=b[i]) && (a[i+1]<=b[i+1]) && (a[i]<=a[i+1]) && (b[i]<=b[i+1]))

//class ParParams;
template<class T> class Matrix;

//------------------------------------------ Vektor ------------------------------------------------
template < class T >
class Vector {
public:
                Vector() {data=0;size[0]=0;}
    int         Init(int* sz) {
                    data+=size[0];
                    if(data) delete[] data;
                    assert(sz[0]<=sz[1]);
                    size[0]=sz[0];size[1]=sz[1];
                    data=new T[size[1]-size[0]+1];
                    if(data==0) return FALSE;
                    data-=size[0];
                    for(int i=size[0];i<=size[1];i++) data[i]=0;
                    return TRUE;
    }
                Vector(int* sz) {
                    data=0;
                    size[0]=0;
                    Init(sz);
                }
    inline T&   operator[] (int i) const {
#ifndef NDEBUG
                    if(!((data + size[0]) && (i >= size[0]) && (i <= size[1]))) {
                        printf("Matrix access at ?,?,%d (%d/%d)\n",i,size[0],size[1]);
                        assert((data + size[0]) && (i >= size[0]) && (i <= size[1]));
                    }
#endif
                    return data[i];
                }
    T           Min() {
                    T min=data[size[1]],minc;
                    for(int i=size[0];i<size[1];i++) if(min>(minc=data[i])) min=minc;
                    return min;
                }
    T           Max() {
                    T max=data[size[1]],maxc;
                    for(int i=size[0];i<size[1];i++) if(max<(maxc=data[i])) max=maxc;
                    return max;
                }
    int         Size(int* bigsize) {return(2*sizeof(int)+sizeof(T)*(bigsize[1]-bigsize[0]+1));}

               ~Vector() {
                    data+=size[0];
                    if(data) delete[] data;
                }
//protected:
    T*          data;
    int         size[2];
    int         kl, kr;
};

//------------------------------------------ Plane ------------------------------------------------
template < class T >
class Plane {
public:
                Plane() {
                    data=0;
                    size[0]=0;
                }
    int         Init(int* sz) {
                    data+=size[0];
                    if(data) delete[] data;
                    assert(sz[0]<=sz[1] && sz[2]<=sz[3]);
                    size[0]=sz[0];size[1]=sz[1];
                    size[2]=sz[2];size[3]=sz[3];
                    data=new Vector<T>[size[1]-size[0]+1];
                    if (data==0) return FALSE;
                    data-=size[0];
                    for(int i=size[0];i<=size[1];i++) data[i].Init(&(size[2]));
                    return TRUE;
                }
                Plane(int* sz) {
                    data=0;
                    size[0]=0;
                    Init(sz);
                }
    inline Vector<T>& operator[] (int i) const{
#ifndef NDEBUG
                    if(!((data + size[0]) && (i >= size[0]) && (i <= size[1]))) {
                        printf("Matrix access at ?,%d,? (%d/%d)\n",i,size[0],size[1]);
                        assert((data + size[0]) && (i >= size[0]) && (i <= size[1]));
                    }
#endif
                    return data[i];
                }
    T           Min() {
                    T min=data[size[1]].Min(),minc;
                    for(int i=size[0];i<size[1];i++) if(min>(minc=data[i].Min())) min=minc;
                    return min;
                }
    T           Max() {
                    T max=data[size[1]].Max(),maxc;
                    for(int i=size[0];i<size[1];i++) if(max<(maxc=data[i].Max())) max=maxc;
                    return max;
                }
    int         Size(int* bigsize) {return(4*sizeof(int)+data[size[0]].Size(&(bigsize[2]))*(bigsize[1]-bigsize[0]+1));}
               ~Plane() {
                    data+=size[0];
                    if(data) delete[] data;
                }
//protected:
    int         size[4];
    Vector<T>*  data;
};

//--------------------------------------- Matrix ---------------------------------------------------------------
template < class T >
class Matrix {

public:
                Matrix() {
                    data=0;
                    size[0]=0;
                }
    int         Init(int* sz) {
                    data+=size[0];
                    if(data) delete[] data;
                    assert(sz[0]<=sz[1] && sz[2]<=sz[3] && sz[4]<=sz[5]);
                    int i;
                    for(i=0;i<6;i++) size[i]=sz[i];
                    data=new Plane<T>[size[1]-size[0]+1];
                    if(data==0) return FALSE;
                    data-=size[0];
                    for(i=size[0];i<=size[1];i++) data[i].Init(&(size[2]));
                    return TRUE;
                }
                Matrix(int* sz) {
                    data=0;
                    size[0]=0;
                    Init(sz);
                }
    inline Plane<T>& operator[] (int i) const {
#ifndef NDEBUG
                    if(!((data + size[0]) && (i >= size[0]) && (i <= size[1]))) {
                        printf("Matrix access at %d,?,? (%d/%d)\n",i,size[0],size[1]);
                        assert((data + size[0]) && (i >= size[0]) && (i <= size[1]));
                    }
#endif
                    return data[i];
                }
    

    T           Min() {
                    T min=data[size[1]].Min(),minc;
                    for(int i=size[0];i<size[1];i++) if(min>(minc=data[i].Min())) min=minc;
                    return min;
                }
    T           Max() {
                    T max=data[size[1]].Max(),maxc;
                    for(int i=size[0];i<size[1];i++) if(max<(maxc=data[i].Max())) max=maxc;
                    return max;
                }
    void        Dump(char* str) {
                    int i,j,k;
                    for(i=size[0];i<=size[1];i++) for(j=size[2];j<=size[3];j++) for(k=size[4];k<=size[5];k++) 
                        printf("%s[%3d][%3d][%3d]=%le\n",str,i,j,k,data[i][j][k]);
                }
               ~Matrix() {
                    data+=size[0];
                    if(data) delete[] data;
                }
//protected:
    int         size[6];
    Plane<T>*   data;
};


inline Vector3 DIV_n(Matrix<Vector3> &A,int i,int j,int k)
{
    return 
           Vector3((A[i][j][k].e1()-A[i-1][j][k].e1())/SCALE,
                   (A[i][j][k].e2()-A[i][j-1][k].e2())/SCALE,
                   (A[i][j][k].e3()-A[i][j][k-1].e3())
                   );
}

inline Vector3 DIV_p(Matrix<Vector3> &A,int i,int j,int k)
{
    return 
           Vector3((A[i+1][j][k].e1()-A[i][j][k].e1())/SCALE,
                   (A[i][j+1][k].e2()-A[i][j][k].e2())/SCALE,
                   (A[i][j][k+1].e3()-A[i][j][k].e3())
                  );
}

inline Vector3 ROATE_p(Matrix<Vector3> &A,int i,int j,int k)
{
    return 
    Vector3( (A[i][j+1][k].e3()-A[i][j][k].e3()-SCALE*A[i][j][k+1].e2()+SCALE*A[i][j][k].e2() )/SCALE,
             (SCALE*A[i][j][k+1].e1()-SCALE*A[i][j][k].e1()-A[i+1][j][k].e3()+A[i][j][k].e3() )/SCALE,
             (A[i+1][j][k].e2()-A[i][j][k].e2()-A[i][j+1][k].e1()+A[i][j][k].e1() )/SCALE
           );
}

inline Vector3 ROATE_n(Matrix<Vector3> &A,int i,int j,int k)
{
    return 
    Vector3( (A[i][j][k].e3()-A[i][j-1][k].e3()-SCALE*A[i][j][k].e2()+SCALE*A[i][j][k-1].e2() )/SCALE,
             (SCALE*A[i][j][k].e1()-SCALE*A[i][j][k-1].e1()-A[i][j][k].e3()+A[i-1][j][k].e3() )/SCALE,
             (A[i][j][k].e2()-A[i-1][j][k].e2()-A[i][j][k].e1()+A[i][j-1][k].e1() )/SCALE
           );
}

inline Vector3 AVERAGE(Matrix<Vector3> &A,int i,int j,int k)
{
    return
    Vector3();
}

inline Vector3 AVERAGE(Matrix<Scalar> &A,int i,int j,int k)
{
    return
    Vector3();
}

//plane
inline Vector3 DIV_n(Plane<Vector3> &A,int i,int j,int k)
{
    return 
           Vector3( 0.0,
                   (A[j][k].e2()-A[j-1][k].e2()),
                   (A[j][k].e3()-A[j][k-1].e3())
                   );
}

inline Vector3 DIV_p(Plane<Vector3> &A,int i,int j,int k)
{
    return 
           Vector3( 0.0,
                   (A[j+1][k].e2()-A[j][k].e2()),
                   (A[j][k+1].e3()-A[j][k].e3())
                  );
}

inline Vector3 ROATE_p(Plane<Vector3> &A,int i,int j,int k)
{
    return 
    Vector3( (A[j+1][k].e3()-A[j][k].e3()-A[j][k+1].e2()+A[j][k].e2() ),
             (A[j][k+1].e1()-A[j][k].e1()-0.0 ),
             (0.0                        -A[j+1][k].e1()+A[j][k].e1() )
           );
}

inline Vector3 ROATE_n(Plane<Vector3> &A,int i,int j,int k)
{
    return 
    Vector3( (A[j][k].e3()-A[j-1][k].e3()-A[j][k].e2()+A[j][k-1].e2() ),
             (A[j][k].e1()-A[j][k-1].e1()-0.0 ),
             (0.0                        -A[j][k].e1()+A[j-1][k].e1() )
           );
}
#endif


