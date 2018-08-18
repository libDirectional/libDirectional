
#ifndef BINGHAM_UTIL_H
#define BINGHAM_UTIL_H

#if defined(_MSC_VER) || defined(__MINGW32__)
#include <malloc.h>     // alloca
#else
#include <alloca.h>     // alloca
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define test_alloc(X) do{ if ((void *)(X) == NULL){ fprintf(stderr, "Out of memory in %s, (%s, line %d).\n", __FUNCTION__, __FILE__, __LINE__); exit(1); }} while (0)
#define safe_calloc(x, n, type) do{ x = (type*)calloc(n, sizeof(type)); test_alloc(x); } while (0)
#define safe_malloc(x, n, type) do{ x = (type*)malloc((n)*sizeof(type)); test_alloc(x); } while (0)
#define safe_realloc(x, n, type) do{ x = (type*)realloc(x,(n)*sizeof(type)); test_alloc(x); } while(0)

double tetrahedron_volume(double x[], double y[], double z[], double w[], int n);   /* calculate the volume of a tetrahedron */

int count(int x[], int n);                                            /* count the non-zero elements of x */
int find(int *k, int x[], int n);                                     /* computes a dense array of the indices of x's non-zero elements */
int findinv(int *k, int x[], int n);                                  /* computes a sparse array of the indices of x's non-zero elements */
double arr_max(double x[], int n);                                    /* computes the max of x */
double arr_min(double x[], int n);                                    /* computes the min of x */
double norm(double x[], int n);                                       /* computes the norm of x */
double dist(double x[], double y[], int n);                           /* computes the norm of x-y */
double dist2(double x[], double y[], int n);                          /* computes the norm^2 of x-y */
double dot(double x[], double y[], int n);                            /* computes the dot product of x and y */
void add(double z[], double x[], double y[], int n);                  /* adds two vectors, z = x+y */
void sub(double z[], double x[], double y[], int n);                  /* subtracts two vectors, z = x-y */
void mult(double y[], double x[], double c, int n);                   /* multiplies a vector by a scalar, y = c*x */
void avg(double z[], double x[], double y[], int n);                  /* averages two vectors, z = (x+y)/2 */

double **new_matrix2(int n, int m);                                         /* create a new n-by-m 2d matrix of doubles */
int **new_matrix2i(int n, int m);                                           /* create a new n-by-m 2d matrix of ints */
void free_matrix2(double **X);                                                  /* free a 2d matrix of doubles */
void free_matrix2i(int **X);                                                    /* free a 2d matrix of ints */
void matrix_copy(double **Y, double **X, int n, int m);                         /* matrix copy, Y = X */
void reorder_rows(double **Y, double **X, int *idx, int n, int m);              /* reorder the rows of X, Y = X(idx,:) */

typedef struct ilist {
  int x;
  int len;
  struct ilist *next;
} ilist_t;

ilist_t *ilist_add(ilist_t *x, int a);                  /* add an element to a list */
int ilist_contains(ilist_t *x, int a);                  /* check if a list contains an element */
int ilist_find(ilist_t *x, int a);                      /* find the index of an element in a list (or -1 if not found) */
void ilist_free(ilist_t *x);                            /* free a list */


typedef struct {
  int index;
  ilist_t *neighbors;
  int *edges;
} vertex_t;

typedef struct {
  int i;
  int j;
} edge_t;

typedef struct {
  int i;
  int j;
  int k;
} face_t;

typedef struct {
  int nv;
  int ne;
  vertex_t *vertices;
  edge_t *edges;
} graph_t;

void graph_free(graph_t *g);                                /* free a graph */
int graph_find_edge(graph_t *g, int i, int j);              /* find the index of an edge in a graph */
void graph_smooth(double **dst, double **src, graph_t *g, int d, double w);  /* smooth the edges of a graph */

typedef struct {
  int nv;
  int ne;
  int nf;
  int *vertices;
  edge_t *edges;
  face_t *faces;
  int *nvn;                /* # vertex neighbors */
  int *nen;                /* # edge neighbors */
  int **vertex_neighbors;  /* vertex -> {vertices} */
  int **vertex_edges;      /* vertex -> {edges} */
  int **edge_neighbors;    /* edge -> {vertices} */
  int **edge_faces;        /* edge -> {faces} */
  /* internal vars */
  int _vcap;
  int _ecap;
  int _fcap;
  int *_vncap;
  int *_encap;
} meshgraph_t;

meshgraph_t *meshgraph_new(int vcap, int dcap);
int meshgraph_find_edge(meshgraph_t *g, int i, int j);
int meshgraph_find_face(meshgraph_t *g, int i, int j, int k);
int meshgraph_add_edge(meshgraph_t *g, int i, int j);
int meshgraph_add_face(meshgraph_t *g, int i, int j, int k);

typedef struct {
  unsigned char r;
  unsigned char g;
  unsigned char b;
} color_t;

extern const color_t colormap[256];

#ifdef __cplusplus
}
#endif

#endif
