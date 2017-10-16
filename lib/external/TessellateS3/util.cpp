#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"

#include <malloc.h>
#ifdef __GNUC__
    #include <alloca.h>
#endif

const color_t colormap[256] =
  {{0, 0, 131},
   {0, 0, 135},
   {0, 0, 139},
   {0, 0, 143},
   {0, 0, 147},
   {0, 0, 151},
   {0, 0, 155},
   {0, 0, 159},
   {0, 0, 163},
   {0, 0, 167},
   {0, 0, 171},
   {0, 0, 175},
   {0, 0, 179},
   {0, 0, 183},
   {0, 0, 187},
   {0, 0, 191},
   {0, 0, 195},
   {0, 0, 199},
   {0, 0, 203},
   {0, 0, 207},
   {0, 0, 211},
   {0, 0, 215},
   {0, 0, 219},
   {0, 0, 223},
   {0, 0, 227},
   {0, 0, 231},
   {0, 0, 235},
   {0, 0, 239},
   {0, 0, 243},
   {0, 0, 247},
   {0, 0, 251},
   {0, 0, 255},
   {0, 4, 255},
   {0, 8, 255},
   {0, 12, 255},
   {0, 16, 255},
   {0, 20, 255},
   {0, 24, 255},
   {0, 28, 255},
   {0, 32, 255},
   {0, 36, 255},
   {0, 40, 255},
   {0, 44, 255},
   {0, 48, 255},
   {0, 52, 255},
   {0, 56, 255},
   {0, 60, 255},
   {0, 64, 255},
   {0, 68, 255},
   {0, 72, 255},
   {0, 76, 255},
   {0, 80, 255},
   {0, 84, 255},
   {0, 88, 255},
   {0, 92, 255},
   {0, 96, 255},
   {0, 100, 255},
   {0, 104, 255},
   {0, 108, 255},
   {0, 112, 255},
   {0, 116, 255},
   {0, 120, 255},
   {0, 124, 255},
   {0, 128, 255},
   {0, 131, 255},
   {0, 135, 255},
   {0, 139, 255},
   {0, 143, 255},
   {0, 147, 255},
   {0, 151, 255},
   {0, 155, 255},
   {0, 159, 255},
   {0, 163, 255},
   {0, 167, 255},
   {0, 171, 255},
   {0, 175, 255},
   {0, 179, 255},
   {0, 183, 255},
   {0, 187, 255},
   {0, 191, 255},
   {0, 195, 255},
   {0, 199, 255},
   {0, 203, 255},
   {0, 207, 255},
   {0, 211, 255},
   {0, 215, 255},
   {0, 219, 255},
   {0, 223, 255},
   {0, 227, 255},
   {0, 231, 255},
   {0, 235, 255},
   {0, 239, 255},
   {0, 243, 255},
   {0, 247, 255},
   {0, 251, 255},
   {0, 255, 255},
   {4, 255, 251},
   {8, 255, 247},
   {12, 255, 243},
   {16, 255, 239},
   {20, 255, 235},
   {24, 255, 231},
   {28, 255, 227},
   {32, 255, 223},
   {36, 255, 219},
   {40, 255, 215},
   {44, 255, 211},
   {48, 255, 207},
   {52, 255, 203},
   {56, 255, 199},
   {60, 255, 195},
   {64, 255, 191},
   {68, 255, 187},
   {72, 255, 183},
   {76, 255, 179},
   {80, 255, 175},
   {84, 255, 171},
   {88, 255, 167},
   {92, 255, 163},
   {96, 255, 159},
   {100, 255, 155},
   {104, 255, 151},
   {108, 255, 147},
   {112, 255, 143},
   {116, 255, 139},
   {120, 255, 135},
   {124, 255, 131},
   {128, 255, 128},
   {131, 255, 124},
   {135, 255, 120},
   {139, 255, 116},
   {143, 255, 112},
   {147, 255, 108},
   {151, 255, 104},
   {155, 255, 100},
   {159, 255, 96},
   {163, 255, 92},
   {167, 255, 88},
   {171, 255, 84},
   {175, 255, 80},
   {179, 255, 76},
   {183, 255, 72},
   {187, 255, 68},
   {191, 255, 64},
   {195, 255, 60},
   {199, 255, 56},
   {203, 255, 52},
   {207, 255, 48},
   {211, 255, 44},
   {215, 255, 40},
   {219, 255, 36},
   {223, 255, 32},
   {227, 255, 28},
   {231, 255, 24},
   {235, 255, 20},
   {239, 255, 16},
   {243, 255, 12},
   {247, 255, 8},
   {251, 255, 4},
   {255, 255, 0},
   {255, 251, 0},
   {255, 247, 0},
   {255, 243, 0},
   {255, 239, 0},
   {255, 235, 0},
   {255, 231, 0},
   {255, 227, 0},
   {255, 223, 0},
   {255, 219, 0},
   {255, 215, 0},
   {255, 211, 0},
   {255, 207, 0},
   {255, 203, 0},
   {255, 199, 0},
   {255, 195, 0},
   {255, 191, 0},
   {255, 187, 0},
   {255, 183, 0},
   {255, 179, 0},
   {255, 175, 0},
   {255, 171, 0},
   {255, 167, 0},
   {255, 163, 0},
   {255, 159, 0},
   {255, 155, 0},
   {255, 151, 0},
   {255, 147, 0},
   {255, 143, 0},
   {255, 139, 0},
   {255, 135, 0},
   {255, 131, 0},
   {255, 128, 0},
   {255, 124, 0},
   {255, 120, 0},
   {255, 116, 0},
   {255, 112, 0},
   {255, 108, 0},
   {255, 104, 0},
   {255, 100, 0},
   {255, 96, 0},
   {255, 92, 0},
   {255, 88, 0},
   {255, 84, 0},
   {255, 80, 0},
   {255, 76, 0},
   {255, 72, 0},
   {255, 68, 0},
   {255, 64, 0},
   {255, 60, 0},
   {255, 56, 0},
   {255, 52, 0},
   {255, 48, 0},
   {255, 44, 0},
   {255, 40, 0},
   {255, 36, 0},
   {255, 32, 0},
   {255, 28, 0},
   {255, 24, 0},
   {255, 20, 0},
   {255, 16, 0},
   {255, 12, 0},
   {255, 8, 0},
   {255, 4, 0},
   {255, 0, 0},
   {251, 0, 0},
   {247, 0, 0},
   {243, 0, 0},
   {239, 0, 0},
   {235, 0, 0},
   {231, 0, 0},
   {227, 0, 0},
   {223, 0, 0},
   {219, 0, 0},
   {215, 0, 0},
   {211, 0, 0},
   {207, 0, 0},
   {203, 0, 0},
   {199, 0, 0},
   {195, 0, 0},
   {191, 0, 0},
   {187, 0, 0},
   {183, 0, 0},
   {179, 0, 0},
   {175, 0, 0},
   {171, 0, 0},
   {167, 0, 0},
   {163, 0, 0},
   {159, 0, 0},
   {155, 0, 0},
   {151, 0, 0},
   {147, 0, 0},
   {143, 0, 0},
   {139, 0, 0},
   {135, 0, 0},
   {131, 0, 0},
   {128, 0, 0}};

// count the non-zero elements of x
int count(int x[], int n)
{
  int i;
  int cnt = 0;
  for (i = 0; i < n; i++)
    if (x[i] != 0)
      cnt++;

  return cnt;
}

// returns a dense array of the indices of x's non-zero elements
int find(int *k, int x[], int n)
{
  int i;
  int cnt = 0;
  for (i = 0; i < n; i++)
    if (x[i] != 0)
      k[cnt++] = i;
  return cnt;
}

// returns a sparse array of the indices of x's non-zero elements
int findinv(int *k, int x[], int n)
{
  int i;
  int cnt = 0;
  for (i = 0; i < n; i++)
    if (x[i] != 0)
      k[i] = cnt++;
  return cnt;
}

// computes the max of x
double arr_max(double x[], int n)
{
  int i;

  double y = x[0];
  for (i = 1; i < n; i++)
    if (x[i] > y)
      y = x[i];

  return y;
}

// computes the min of x
double arr_min(double x[], int n)
{
  int i;

  double y = x[0];
  for (i = 1; i < n; i++)
    if (x[i] < y)
      y = x[i];

  return y;
}

// computes the norm of x
double norm(double x[], int n)
{
  double d = 0.0;
  int i;

  for (i = 0; i < n; i++)
    d += x[i]*x[i];

  return sqrt(d);
}

// computes the norm of x-y
double dist(double x[], double y[], int n)
{
  double d = 0.0;
  int i;

  for (i = 0; i < n; i++)
    d += (x[i]-y[i])*(x[i]-y[i]);

  return sqrt(d);
}

// computes the norm^2 of x-y
double dist2(double x[], double y[], int n)
{
  double d = 0.0;
  int i;

  for (i = 0; i < n; i++)
    d += (x[i]-y[i])*(x[i]-y[i]);

  return d;
}

// computes the dot product of z and y
double dot(double x[], double y[], int n)
{
  int i;
  double z = 0.0;
  for (i = 0; i < n; i++)
    z += x[i]*y[i];
  return z;
}

// adds two vectors, z = x+y
void add(double z[], double x[], double y[], int n)
{
  int i;
  for (i = 0; i < n; i++)
    z[i] = x[i] + y[i];
}

// subtracts two vectors, z = x-y
void sub(double z[], double x[], double y[], int n)
{
  int i;
  for (i = 0; i < n; i++)
    z[i] = x[i] - y[i];
}

// multiplies a vector by a scalar, y = c*x
void mult(double y[], double x[], double c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    y[i] = c*x[i];
}

// averages two vectors, z = (x+y)/2
void avg(double z[], double x[], double y[], int n)
{
  add(z, x, y, n);
  mult(z, z, .5, n);
}

// add an element to the front of a list
ilist_t *ilist_add(ilist_t *x, int a)
{
  ilist_t *head;
  safe_malloc(head, 1, ilist_t);
  head->x = a;
  head->next = x;
  head->len = (x ? 1 + x->len : 1);

  return head;
}

// check if a list contains an element
int ilist_contains(ilist_t *x, int a)
{
  if (!x)
    return 0;

  ilist_t *tmp;
  for (tmp = x; tmp; tmp = tmp->next)
    if (tmp->x == a)
      return 1;
  return 0;
}

// find the index of an element in a list (or -1 if not found)
int ilist_find(ilist_t *x, int a)
{
  int i = 0;
  ilist_t *tmp;
  for (tmp = x; tmp; tmp = tmp->next) {
    if (tmp->x == a)
      return i;
    i++;
  }

  return -1;
}

// free a list
void ilist_free(ilist_t *x)
{
  ilist_t *tmp, *tmp2;
  tmp = x;
  while (tmp) {
    tmp2 = tmp->next;
    free(tmp);
    tmp = tmp2;
  }  
}

// create a new n-by-m 2d matrix of doubles
double **new_matrix2(int n, int m)
{
  if (n*m == 0) return NULL;
  int i;
  double *raw, **X;
  safe_calloc(raw, n*m, double);
  safe_malloc(X, n, double*);

  for (i = 0; i < n; i++)
    X[i] = raw + m*i;

  return X;
}

// create a new n-by-m 2d matrix of ints
int **new_matrix2i(int n, int m)
{
  if (n*m == 0) return NULL;
  int i, *raw, **X;
  safe_calloc(raw, n*m, int);
  safe_malloc(X, n, int*);

  for (i = 0; i < n; i++)
    X[i] = raw + m*i;

  return X;
}

// free a 2d matrix of doubles
void free_matrix2(double **X)
{
  if (X == NULL) return;
  free(X[0]);
  free(X);
}

// free a 2d matrix of ints
void free_matrix2i(int **X)
{
  if (X == NULL) return;
  free(X[0]);
  free(X);
}

// calculate the volume of a tetrahedron
double tetrahedron_volume(double x1[], double x2[], double x3[], double x4[], int n)
{
  double U = dist2(x1, x2, n);
  double V = dist2(x1, x3, n);
  double W = dist2(x2, x3, n);
  double u = dist2(x3, x4, n);
  double v = dist2(x2, x4, n);
  double w = dist2(x1, x4, n);

  double a = v+w-U;
  double b = w+u-V;
  double c = u+v-W;

  return sqrt( (4*u*v*w - u*a*a - v*b*b - w*c*c + a*b*c) ) / 12.0 ;
}

// matrix copy, Y = X 
void matrix_copy(double **Y, double **X, int n, int m)
{
  memcpy(Y[0], X[0], n*m*sizeof(double));
}

// reorder the rows of X, Y = X(idx,:)
void reorder_rows(double **Y, double **X, int *idx, int n, int m)
{
  int i;
  double **Y2 = (X==Y ? new_matrix2(n,m) : Y);
  for (i = 0; i < n; i++)
    memcpy(Y2[i], X[idx[i]], m*sizeof(double));
  if (X==Y) {
    matrix_copy(Y, Y2, n, m);
    free_matrix2(Y2);
  }
}

// free a graph
void graph_free(graph_t *g)
{
  int i;
  free(g->edges);
  for (i = 0; i < g->nv; i++) {
    free(g->vertices[i].edges);
    ilist_free(g->vertices[i].neighbors);
  }
  free(g);
}

// find the index of an edge in a graph
int graph_find_edge(graph_t *g, int i, int j)
{
  int k = ilist_find(g->vertices[i].neighbors, j);

  if (k < 0)
    return -1;

  return g->vertices[i].edges[k];
}

/*
 * Create a new meshgraph with initial vertex capacity 'vcap' and degree capacity 'dcap'.
 */
meshgraph_t *meshgraph_new(int vcap, int dcap)
{
  int i;
  meshgraph_t *g;
  safe_calloc(g, 1, meshgraph_t);

  safe_malloc(g->vertices, vcap, int);
  safe_malloc(g->edges, vcap, edge_t);
  safe_malloc(g->faces, vcap, face_t);

  g->_vcap = g->_ecap = g->_fcap = vcap;
  safe_malloc(g->_vncap, vcap, int);
  safe_malloc(g->_encap, vcap, int);

  safe_malloc(g->vertex_neighbors, vcap, int *);
  safe_malloc(g->vertex_edges, vcap, int *);
  safe_calloc(g->nvn, vcap, int);
  for (i = 0; i < vcap; i++) {
    safe_malloc(g->vertex_neighbors[i], dcap, int);
    safe_malloc(g->vertex_edges[i], dcap, int);
    g->_vncap[i] = dcap;
  }

  safe_malloc(g->edge_neighbors, vcap, int *);
  safe_malloc(g->edge_faces, vcap, int *);
  safe_calloc(g->nen, vcap, int);
  for (i = 0; i < vcap; i++) {
    safe_malloc(g->edge_neighbors[i], dcap, int);
    safe_malloc(g->edge_faces[i], dcap, int);
    g->_encap[i] = dcap;
  }

  return g;
}

void meshgraph_free(meshgraph_t *g)
{
  int i;

  free(g->vertices);
  free(g->edges);
  free(g->faces);
  free(g->nvn);
  free(g->nen);
  free(g->_vncap);
  free(g->_encap);

  for (i = 0; i < g->nv; i++) {
    free(g->vertex_neighbors[i]);
    free(g->vertex_edges[i]);
  }
  free(g->vertex_neighbors);
  free(g->vertex_edges);

  for (i = 0; i < g->ne; i++) {
    free(g->edge_neighbors[i]);
    free(g->edge_faces[i]);
  }
  free(g->edge_neighbors);
  free(g->edge_faces);

  free(g);
}

int meshgraph_find_edge(meshgraph_t *g, int i, int j)
{
  int n;
  for (n = 0; n < g->nvn[i]; n++)
    if (g->vertex_neighbors[i][n] == j)
      return g->vertex_edges[i][n];

  return -1;
}

int meshgraph_find_face(meshgraph_t *g, int i, int j, int k)
{
  int e = meshgraph_find_edge(g, i, j);
  if (e < 0)
    return -1;

  int n;
  for (n = 0; n < g->nen[e]; n++)
    if (g->edge_neighbors[e][n] == k)
      return g->edge_faces[e][n];

  return -1;
}

static inline int meshgraph_add_vertex_neighbor(meshgraph_t *g, int i, int vertex, int edge)
{
  int n = g->nvn[i];
  if (n == g->_vncap[i]) {
    g->_vncap[i] *= 2;
    safe_realloc(g->vertex_neighbors[i], g->_vncap[i], int);
    safe_realloc(g->vertex_edges[i], g->_vncap[i], int);
  }
  g->vertex_neighbors[i][n] = vertex;
  g->vertex_edges[i][n] = edge;
  g->nvn[i]++;

  return n;
}

int meshgraph_add_edge(meshgraph_t *g, int i, int j)
{
  int edge = meshgraph_find_edge(g, i, j);
  if (edge >= 0)
    return edge;

  // add the edge
  if (g->ne == g->_ecap) {
    int old_ecap = g->_ecap;
    g->_ecap *= 2;
    safe_realloc(g->edges, g->_ecap, edge_t);
    safe_realloc(g->edge_neighbors, g->_ecap, int *);
    safe_realloc(g->edge_faces, g->_ecap, int *);
    safe_realloc(g->nen, g->_ecap, int);
    safe_realloc(g->_encap, g->_ecap, int);

    int e;
    for (e = old_ecap; e < g->_ecap; e++) {
      g->nen[e] = 0;
      int dcap = g->_encap[0];
      safe_malloc(g->edge_neighbors[e], dcap, int);
      safe_malloc(g->edge_faces[e], dcap, int);
      g->_encap[e] = dcap;
    }
  }

  edge = g->ne;
  g->edges[edge].i = i;
  g->edges[edge].j = j;
  g->ne++;

  // add the vertex neighbors
  meshgraph_add_vertex_neighbor(g, i, j, edge);
  meshgraph_add_vertex_neighbor(g, j, i, edge);

  return edge;
}

static inline int meshgraph_add_edge_neighbor(meshgraph_t *g, int i, int vertex, int face)
{
  int n = g->nen[i];

  if (n == g->_encap[i]) {
    g->_encap[i] *= 2;
    safe_realloc(g->edge_neighbors[i], g->_encap[i], int);
    safe_realloc(g->edge_faces[i], g->_encap[i], int);
  }

  g->edge_neighbors[i][n] = vertex;
  g->edge_faces[i][n] = face;
  g->nen[i]++;

  return n;
}

int meshgraph_add_face(meshgraph_t *g, int i, int j, int k)
{
  int face = meshgraph_find_face(g, i, j, k);
  if (face >= 0)
    return face;

  // add the edges
  int edge_ij = meshgraph_add_edge(g, i, j);
  int edge_ik = meshgraph_add_edge(g, i, k);
  int edge_jk = meshgraph_add_edge(g, j, k);

  // add the face
  if (g->nf == g->_fcap) {
    g->_fcap *= 2;
    safe_realloc(g->faces, g->_fcap, face_t);
  }
  face = g->nf;
  g->faces[face].i = i;
  g->faces[face].j = j;
  g->faces[face].k = k;
  g->nf++;

  // add the edge neighbors
  meshgraph_add_edge_neighbor(g, edge_ij, k, face);
  meshgraph_add_edge_neighbor(g, edge_ik, j, face);
  meshgraph_add_edge_neighbor(g, edge_jk, i, face);

  return face;
}


