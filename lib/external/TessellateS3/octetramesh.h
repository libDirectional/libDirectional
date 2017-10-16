#ifndef BINGHAM_OCTETRAMESH_H
#define BINGHAM_OCTETRAMESH_H


#include "tetramesh.h"
#include "util.h"


/** Tools for subdividing and smoothing a tetrahedral-octahedral (alternating cubic) mesh **/


typedef struct {
  int nv;              // number of vertices
  int nt;              // number of tetrahedra
  int no;              // number of octahedra
  int d;               // dimension of each vertex
  double **vertices;   // vertices
  int **tetrahedra;    // indices into vertices
  int **octahedra;     // indices into vertices
} octetramesh_t;

typedef struct {
  int num_vertices;
  int num_edges;
  int num_tetrahedra;
  int num_octahedra;
  double min_edge_len;
  double max_edge_len;
  double avg_edge_len;
  double std_edge_len;
  double min_tetra_skewness;
  double max_tetra_skewness;
  double avg_tetra_skewness;
  double std_tetra_skewness;
  double min_tetra_volume;
  double max_tetra_volume;
  double avg_tetra_volume;
  double std_tetra_volume;
} octetramesh_stats_t;



// Convert an octa-tetrahedral mesh to a tetrahedral mesh.
tetramesh_t *octetramesh_to_tetramesh(octetramesh_t *mesh);

// Build the graph of a mesh.
graph_t *octetramesh_graph(octetramesh_t *mesh);

/** Subdivide each octahedron into 6 octahedra and 8 tetrahedra,
    and subdivide each tetrahedron into 4 tetrahedra and an octahedron.*/
void octetramesh_subdivide(octetramesh_t *dst, octetramesh_t *src);

/** Subdivide each specificied octahedron into 6 octahedra and 8 tetrahedra,
    and subdivide each specified tetrahedron into 4 tetrahedra and an octahedron.*/
void octetramesh_subdivide_select(octetramesh_t *dst, octetramesh_t *src, int tetmask[], int octmask[]);

/** Multi-resolution subdivision based on approximating
    a scalar function f to a given resolution.*/
void octetramesh_subdivide_mres(octetramesh_t *dst, octetramesh_t *src,
				double(*f)(double *, void *), void *fdata, double resolution);

// Create (allocate) the contents of a mesh.
void octetramesh_new(octetramesh_t *mesh, int nv, int nt, int no, int d);

// Copy the contents of one mesh into another mesh.
void octetramesh_copy(octetramesh_t *dst, octetramesh_t *src);
void octetramesh_copy_vertices(octetramesh_t *dst, octetramesh_t *src);
void octetramesh_copy_tetrahedra(octetramesh_t *dst, octetramesh_t *src);
void octetramesh_copy_octahedra(octetramesh_t *dst, octetramesh_t *src);

// Clone a mesh
octetramesh_t *octetramesh_clone(octetramesh_t *src);

// Free the contents of a octetrahedral mesh.
void octetramesh_free(octetramesh_t *mesh);

// Save a octetrahedral mesh to PLY file.
void octetramesh_save_PLY(octetramesh_t *mesh, char *filename);

// Compute stats of a mesh.
octetramesh_stats_t octetramesh_stats(octetramesh_t *T);

// Print the stats of a mesh.
void octetramesh_print_stats(octetramesh_stats_t);

#endif
