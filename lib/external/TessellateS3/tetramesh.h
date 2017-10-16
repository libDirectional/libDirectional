#ifndef BINGHAM_TETRAMESH_H
#define BINGHAM_TETRAMESH_H

#include "util.h"


/** Tools for subdividing and smoothing a tetrahedral mesh **/

typedef struct {
  int nv;              /* number of vertices */
  int nt;              /* number of tetrahedra */
  int d;               /* dimension of each vertex */
  double **vertices;   /* vertices */
  int **tetrahedra;    /* indices into vertices */
} tetramesh_t;

typedef struct {
  int num_vertices;
  int num_edges;
  int num_tetrahedra;
  double min_edge_len;
  double max_edge_len;
  double avg_edge_len;
  double std_edge_len;
  double min_skewness;
  double max_skewness;
  double avg_skewness;
  double std_skewness;
  double min_volume;
  double max_volume;
  double avg_volume;
  double std_volume;
} tetramesh_stats_t;


/* Convert a tetramesh to a set of centroids and volumes. */
void tetramesh_centroids(double **centroids, double *volumes, tetramesh_t *mesh);

/* Subdivide each tetrahedron in a mesh into 8 smaller tetrahedra. */
void tetramesh_subdivide(tetramesh_t *dst, tetramesh_t *src);

/* Smooth the tetrahedral mesh */
void tetramesh_smooth(tetramesh_t *dst, tetramesh_t *src, double w);
void tetramesh_smooth2(tetramesh_t *dst, tetramesh_t *src, double w);
void tetramesh_smooth_edges(tetramesh_t *dst, tetramesh_t *src, double w);

/* Create (allocate) the contents of a mesh. */
void tetramesh_new(tetramesh_t *mesh, int nv, int nt, int d);

/* Copy the contents of one mesh into another mesh. */
void tetramesh_copy(tetramesh_t *dst, tetramesh_t *src);
void tetramesh_copy_vertices(tetramesh_t *dst, tetramesh_t *src);
void tetramesh_copy_tetrahedra(tetramesh_t *dst, tetramesh_t *src);

/* Clone a mesh */
tetramesh_t *tetramesh_clone(tetramesh_t *src);

/* Free the contents of a tetrahedral mesh. */
void tetramesh_free(tetramesh_t *mesh);

/* Get the meshgraph of a tetramesh. */
meshgraph_t *tetramesh_meshgraph(tetramesh_t *mesh);

/* Get the graph of a tetramesh. */
graph_t *tetramesh_graph(tetramesh_t *mesh);

/* Save a tetrahedral mesh to PLY file. */
void tetramesh_save_PLY(tetramesh_t *mesh, meshgraph_t *graph, char *filename);
void tetramesh_save_PLY_colors(tetramesh_t *mesh, meshgraph_t *graph, char *filename, int *colors);


/* Compute stats of a mesh. */
tetramesh_stats_t tetramesh_stats(tetramesh_t *T);

/* Print the stats of a mesh. */
void tetramesh_print_stats(tetramesh_stats_t);

#endif
